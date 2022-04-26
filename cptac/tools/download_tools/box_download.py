import os
import requests
import pandas as pd
import requests
import threading
import warnings
import webbrowser
import bs4
import getpass
import logging

from queue import Queue
from werkzeug import Request, Response
from werkzeug.serving import make_server

from cptac.tools.download_tools.download_tools import *
from cptac.exceptions import DownloadFailedError

# making box token global here allows for it to be used repeatedly as the box_download function gets called over and over.
# eventually authentication will be deprecated, otherwise we would clean up the authentication process
__BOX_TOKEN__ = None

# Some websites don't like requests from sources without a user agent. Let's preempt that issue.
USER_AGENT = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0)'
HEADERS = {'User-Agent': USER_AGENT}
 
def box_download(cancer, source, datatypes, version, redownload):
   
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
    if datatypes != "all" and source != "pdc":
        version_index = get_filtered_version_index(version_index=version_index, datatypes=datatypes, source=dataset, version=version_number)

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
            # check if the required file is from a source whose files are stored on Box.com
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

def gather_files(version_path, version_index, redownload):
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