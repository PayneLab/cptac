#   Copyright 2018 Samuel Payne sam_payne@byu.edu
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import bs4
import getpass
import logging
import os
import requests
import threading
import warnings
import webbrowser

import cptac
from cptac.tools.file_tools import *
from cptac.exceptions import CptacDevError, DownloadFailedError, InvalidParameterError, NoInternetError, PdcDownloadError
from queue import Queue
from werkzeug import Request, Response
from werkzeug.serving import make_server

# making box token global here allows for it to be used repeatedly as the _stream function gets called over and over.
# eventually authentication will be deprecated, otherwise we would make the authentication through box cleaner
__BOX_TOKEN__ = None

# TODO: Move Study ids map to another file?
STUDY_IDS_MAP = {
    "pdcbrca": {
        "acetylome": "PDC000239", # Prospective Breast BI Acetylome
        "phosphoproteome": "PDC000121", # Prospective BRCA Phosphoproteome S039-2
        "proteome": "PDC000120", # Prospective BRCA Proteome S039-1
    },
    "pdcccrcc": {
        "phosphoproteome": "PDC000128", # CPTAC CCRCC Discovery Study - Phosphoproteme S044-2
        "proteome": "PDC000127", # CPTAC CCRCC Discovery Study - Proteome S044-1
    },
    "pdccoad": {
        "phosphoproteome": "PDC000117", # Prospective COAD Phosphoproteome S037-3
        "proteome": "PDC000116", # Prospective COAD Proteome S037-2
    },
    "pdcgbm": {
        "acetylome": "PDC000245", # CPTAC GBM Discovery Study - Acetylome
        "phosphoproteome": "PDC000205", # CPTAC GBM Discovery Study - Phosphoproteome
        "proteome": "PDC000204", # CPTAC GBM Discovery Study - Proteome
    },
    "pdchnscc": {
        "phosphoproteome": "PDC000222", # CPTAC HNSCC Discovery Study - Phosphoproteome
        "proteome": "PDC000221", # CPTAC HNSCC Discovery Study - Proteome
    },
    "pdclscc": {
        "acetylome": "PDC000233", # CPTAC LSCC Discovery Study - Acetylome
        "phosphoproteome": "PDC000232", # CPTAC LSCC Discovery Study - Phosphoproteome
        "proteome": "PDC000234", # CPTAC LSCC Discovery Study - Proteome
        "ubiquitylome": "PDC000237", # CPTAC LSCC Discovery Study - Ubiquitylome
    },
    "pdcluad": {
        "acetylome": "PDC000224", # CPTAC LUAD Discovery Study - Acetylome
        "phosphoproteome": "PDC000149", # CPTAC LUAD Discovery Study - Phosphoproteome
        "proteome": "PDC000153", # CPTAC LUAD Discovery Study - Proteome
    },
    "pdcov": {
        "phosphoproteome": "PDC000119", # Prospective OV Phosphoproteome S038-3
        "proteome": "PDC000118", # Prospective OV Proteome S038-2
    },
    "pdcpdac": {
        "proteome": "PDC000270", # CPTAC PDAC Discovery Study - Proteome
        "phosphoproteome": "PDC000271", # CPTAC PDAC Discovery Study - Phosphoproteome
    },
    "pdcucec": {
        "acetylome": "PDC000226", # CPTAC UCEC Discovery Study - Acetylome
        "phosphoproteome": "PDC000126", # UCEC Discovery - Phosphoproteome S043-2
        "proteome": "PDC000125", # UCEC Discovery - Proteome S043-1
    },
}

# Some websites don't like requests from sources without a user agent. Let's preempt that issue.
USER_AGENT = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0)'
HEADERS = {'User-Agent': USER_AGENT}

def download(sources, cancers='all', version="latest", redownload=False):
    """Download data files for the specified cancers, sources, and datatypes. Defaults to downloading latest version on server.

    Parameters:
    sources (dict): Keys are source names (str), values are the datatypes (list of str)
    cancers (list of str): The cancers for which the sources/datatypes will be downloaded
    version (str, optional): Which version of the data files to download. Defaults to latest on server.
    redownload (bool, optional): Whether to redownload the data files, even if that version of the data is already downloaded. Default False.

    Returns:
    bool: Indicates whether download was successful.
    """

    # check if sources parameter is valid
    _validate_sources(sources)
    
    # check if cancers parameter is valid
    if cancers != 'all':
        _validate_cancers(cancers)

    # iterate through cancers and sources and download corresonding data files
    if cancers == 'all': cancers = cptac.get_cancer_options()
    success = True
    for cancer in cancers:
        for source, datatypes in sources.items():
            if type(datatypes) is not list:
                datatypes = list([datatypes])
            if not _stream(cancer, source, datatypes, version=version, redownload=redownload):
                success = False

    return success


def _stream(cancer, source, datatypes, version, redownload):
   
    dataset = source + "_" + cancer

    # Get our dataset path
    dataset_path = get_dataset_path(dataset)

    # Update the index
    update_index(dataset=dataset)

    # Load the index
    index = get_index(dataset)

    # Validate the version number, including parsing if it's "latest"
    version_number = validate_version(version, dataset, use_context="download")

    # Construct the path to the directory for this version
    version_path = os.path.join(dataset_path, f"{dataset}_v{version_number}")

    # Get the index for the desired version
    # If datatypes are specified, filter out the undesired datatypes
    version_index = index.get(version_number)
    if datatypes != "all":
        version_index = get_filtered_version_index(version_index=version_index, datatypes=datatypes, source=dataset)

    # Get list of files to download.
    files_to_download = _gather_files(version_path=version_path, version_index=version_index, redownload=redownload)
   
    # Retrurn true if no new files to download
    if files_to_download is None: 
        return True

    # Else Download the files
    password = None
    total_files = len(files_to_download)

    # verify box token or login to box and get a new token
    if source != 'awg':
        _authenticate()

    for data_file in files_to_download:

        file_index = version_index.get(data_file)
        server_hash = file_index.get("hash")
        file_url = file_index.get("url")

        file_path = os.path.join(version_path, data_file)
        file_number = files_to_download.index(data_file) + 1

        downloaded_path = download_file(file_url, file_path, server_hash,source=source, password=password, file_message=f"{dataset} v{version} data files", file_number=file_number, total_files=total_files)

        while downloaded_path == "wrong_password":
            if password is None:
                password = getpass.getpass(prompt=f'Password for {dataset} dataset: ') # We manually specify the prompt parameter so it shows up in Jupyter Notebooks
            else:
                password = getpass.getpass(prompt="Wrong password. Try again: ")
            print("\033[F", end='\r') # Use an ANSI escape sequence to move cursor back up to the beginning of the last line, so in the next line we can clear the password prompt
            print("\033[K", end='\r') # Use an ANSI escape sequence to print a blank line, to clear the password prompt

            downloaded_path = download_file(file_url, file_path, source, server_hash, password=password, file_message=f"{dataset} v{version} data files", file_number=file_number, total_files=total_files)

    return True

def _gather_files(version_path, version_index, redownload):
    if os.path.isdir(version_path): # See if they've downloaded this version before. 
        if redownload:
            files_to_download = list(version_index.keys())
        else:
            files_to_download = []
            for data_file in version_index.keys():
                # Compare the server and local hashes, to make sure there was no data corruption
                file_path = os.path.join(version_path, data_file)
                if os.path.isfile(file_path):
                    file_index = version_index.get(data_file)
                    server_hash = file_index.get("hash")
                    local_hash = hash_file(file_path) 
                    if local_hash == server_hash:
                        continue
                files_to_download.append(data_file)

            if len(files_to_download) == 0:
                return None
    else:
        os.mkdir(version_path)
        files_to_download = list(version_index.keys())

    return files_to_download

def _authenticate():
    global __BOX_TOKEN__
    if __BOX_TOKEN__ is None:
        __BOX_TOKEN__ = get_box_token()

def _validate_sources(sources):
    if type(sources) is not dict:
        raise InvalidParameterError("Sources must be a dict of form {'source':['datatypes']}. 'all' is a valid source and datatype.")
    
    valid_sources = cptac.get_source_options()
    for s in sources:
        if s not in valid_sources:
            raise InvalidParameterError(f"{s} is not a valid source! Call cptac.list_datasets() for valid options.")
    

def _validate_cancers(cancers):
    if type(cancers) is str and cancers == 'all':
        pass

def update_index(dataset):
    """Check if the index of the given dataset is up to date with server version, and update it if needed.

    Parameters:
    dataset (str): The name of the dataset to check the index of.

    Returns:
    bool: Indicates if we were able to check the index and update if needed (i.e. we had internet)
    """
    # Get the path to our dataset
    dataset_path = get_dataset_path(dataset)

    # Define our file names we'll need
    index_urls_file = "index_urls.tsv"
    index_hash_file = "index_hash.txt"
    index_file = "index.txt"

    # Get, from the server, what the md5 hash of our index file should be
    index_urls_path = os.path.join(dataset_path, index_urls_file)
    urls_dict = parse_tsv_dict(index_urls_path)
    index_hash_url = urls_dict.get(index_hash_file)

    checking_msg = f"Checking that {dataset} index is up-to-date..."
    print(checking_msg, end='\r')
    try:
        server_index_hash = download_text(index_hash_url)
    finally:
        print(" " * len(checking_msg), end='\r') # Erase the checking message, even if there was an internet error

    index_path = os.path.join(dataset_path, index_file)

    if os.path.isfile(index_path):
        local_index_hash = hash_file(index_path)
        if local_index_hash == server_index_hash:
            return True

    index_url = urls_dict.get(index_file)
    download_file(url=index_url, path=index_path, server_hash=server_index_hash, file_message=f"{dataset} index")

    if os.path.isfile(index_path):
        local_index_hash = hash_file(index_path)
        if local_index_hash == server_index_hash:
            return True

    # If we get here, something apparently went wrong with the download.
    raise NoInternetError("Insufficient internet. Check your internet connection.")

def download_text(url):
    """Download text from a direct download url for a text file.

    Parameters:
    url (str): The direct download url for the text.

    Returns:
    str: The downloaded text.
    """
    
    try:
        response = requests.get(url, headers=HEADERS, allow_redirects=True)
        response.raise_for_status() # Raises a requests HTTPError if the response code was unsuccessful
    except requests.RequestException: # Parent class for all exceptions in the requests module
        raise NoInternetError("Insufficient internet. Check your internet connection.") from None 

    text = response.text.strip()
    return text

def download_file(url, path, server_hash, source=None, password=None, file_message=None, file_number=None, total_files=None): 
    """Download a file from a given url to the specified location.

    Parameters:
    url (str): The direct download url for the file.
    path (str): The path to the file (not just the directory) to save the file to on the local machine.
    source(str): The source the file is coming from to help determine authentication needs
    server_hash (str): The hash for the file, to check it against. If check fails, try download one more time, then throw an exception.
    password (str, optional): If the file is password protected, the password for it. Unneeded otherwise.
    file_message (str, optional): Identifing message about the file, to be printed while it's downloading. Default None will cause the full file name to be printed.
    file_number (int, optional): Which file this is in a batch of files, if you want to print a "File 1/15", "File 2/15", etc. sort of message. Must also pass total_files parameter.
    total_files (int, optional): The total number of files in the download batch, if you're printing that. Must also pass file_number parameter.

    Returns:
    str: The path the file was downloaded to.
    """
    # We provide the option of displaying a message indicating which file this is in a batch of files we're currently downloading
    batch_status = ''
    if (file_number is not None) and (total_files is not None):
        batch_status = f" ({file_number}/{total_files})"

    if file_message is None:
        file_message = path.split(os.sep)[-1]

    download_msg = f"Downloading {file_message}{batch_status}..."
    print(download_msg, end='\r')

    for i in range(2):
        try:
            # check if the required file is from a pancan data source or not (except pdc)
            if source in ["bcm", "broad", "harmonized", "mssm", "umich", "washu"]: # We are using Box OAuth2
                download_url = f"https://api.box.com/2.0/files/{url}/content/" # url is actually file ID
                headers = dict(HEADERS)
                global __BOX_TOKEN__
                headers["Authorization"] = f"Bearer {__BOX_TOKEN__}"
                response = requests.get(download_url, headers=headers)
            
            elif password is None: # No password or OAuth2 (awg files)
                response = requests.get(url, headers=HEADERS, allow_redirects=True)

            else: # The file is password protected (awgconf files)
                with requests.Session() as session: # Use a session object to save cookies
                    # Construct the urls for our GET and POST requests
                    get_url = url
                    post_url = get_url.replace("https://byu.box.com/shared", "https://byu.app.box.com/public")

                    # Send initial GET request and parse the request token out of the response
                    get_response = session.get(get_url, headers=HEADERS) 
                    soup = bs4.BeautifulSoup(get_response.text, "html.parser")
                    token_tag = soup.find(id="request_token")
                    token = token_tag.get("value")

                    # Send a POST request, with the password and token, to get the data
                    payload = {
                        'password': password,
                        'request_token': token}
                    response = session.post(post_url, headers=HEADERS, data=payload)

            response.raise_for_status() # Raises a requests.HTTPError if the response code was unsuccessful
        except requests.RequestException as e: # Parent class for all exceptions in the requests module
            raise Exception(e) #from None
            
        local_hash = hash_bytes(response.content)
        if local_hash == server_hash: # Only replace the old file if the new one downloaded successfully.
            with open(path, 'wb') as dest:
                dest.write(response.content)
            print(" " * len(download_msg), end='\r') # Erase the downloading message
            return path
        elif response.text.strip().startswith("<!DOCTYPE html>"): # The password was wrong, so we just got a webpage
            print(" " * len(download_msg), end='\r') # Erase the downloading message
            return "wrong_password"

    # If we get to this point, the download failed.
    file_name = path.split(os.sep)[-1]
    raise DownloadFailedError(f"Download failed for {file_name}.")

def get_box_token():

    @Request.application
    def receive(request):
        q.put(request.args.get('code'))
        return Response("Authentication successful. You can close this window.", 200)

    # Don't show logs from server
    log = logging.getLogger('werkzeug')
    log.disabled = True

    # Set up authentication parameters
    base_url = "https://account.box.com/api/oauth2/authorize"
    client_id = "kztczhjoq3oes38yywuyfp4t9tu11it8"
    client_secret = "a5xNE1qj4Z4H3BSJEDVfzbxtmxID6iKY"
    login_url = f"{base_url}?client_id={client_id}&response_type=code"

    q = Queue()
    s = make_server("localhost", 8003, receive)
    t = threading.Thread(target=s.serve_forever)
    t.start()

    # Send the user to the "Grant access" page
    webbrowser.open(login_url)
    login_msg = "Please login to Box on the webpage that was just opened and grant access for cptac to download files through your account. If you accidentally closed the browser window, press Ctrl+C and call the download function again."
    print(login_msg)

    temp_code = q.get(block=True)
    s.shutdown()
    t.join()

    # Use the temporary access code to get the long term access token
    token_url = "https://api.box.com/oauth2/token";

    params = {
       'grant_type': 'authorization_code',
       'code': temp_code,
       'client_id': client_id,
       'client_secret': client_secret,
    }

    auth_resp = requests.post(token_url, data=params)
    access_token = auth_resp.json()["access_token"]

    return access_token