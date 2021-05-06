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

from flask import Flask, cli, request
from threading import Thread
from pathlib import Path

import webbrowser
import time
import logging
import os
import requests
import getpass
import bs4

from .file_tools import *
from .exceptions import InvalidParameterError, NoInternetError, DownloadFailedError

# Some websites don't like requests from sources without a user agent. Let's preempt that issue.
USER_AGENT = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0)'
HEADERS = {'User-Agent': USER_AGENT}

# For a rudimentary data sharing between processes
path_here = os.path.abspath(os.path.dirname(__file__))
LOCK_PATH = os.path.join(path_here, "lock.tmp")
CODE_FILE_PATH = os.path.join(path_here, "code.tmp")

def download(dataset, version="latest", redownload=False, _box_auth=False, _box_token=None):
    """Download data files for the specified datasets. Defaults to downloading latest version on server.

    Parameters:
    dataset (str): The name of the dataset to download data for, or "all" to download data for all datasets
    version (str, optional): Which version of the data files to download. Defaults to latest on server.
    redownload (bool, optional): Whether to redownload the data files, even if that version of the data is already downloaded. Default False.
    _box_auth (bool, optional): Whether to download the files using Box file IDs and OAuth2 authentication. Default False.
    _box_token (str, optional): The OAuth2 token for Box, if already generated. Default None.

    Returns:
    bool: Indicates whether download was successful.
    """

    dataset = dataset.lower()

    # Check that they're using the right download function
    datasets = [
        "brca",
        "ccrcc",
        "colon",
        "endometrial",
        "gbm",
        "hnscc",
        "lscc",
        "luad",
        "ovarian",
        "pdac",
        "ucecconf",
    ]

    if dataset in datasets and _box_auth:
        raise InvalidParameterError(f"You are trying to use the cptac.pancan.download function to download the {dataset} dataset. This is the wrong function; that dataset is not associated with the cptac.pancan module. Please instead use the regular cptac.download function.")

    # Process the optional "all" parameter
    if dataset == "all":
        overall_result = True
        for dataset in datasets:
            if not download(dataset, redownload=redownload):
                overall_result = False

        return overall_result

    # Get our dataset path
    dataset_path = get_dataset_path(dataset)

    # Update the index
    update_index(dataset)

    # Load the index
    index = get_index(dataset)

    # Validate the version number, including parsing if it's "latest"
    version = validate_version(version, dataset, use_context="download")

    # Construct the path to the directory for this version
    version_path = os.path.join(dataset_path, f"{dataset}_v{version}")

    # See if they've downloaded this version before. Get list of files to download.
    version_index = index.get(version)
    if os.path.isdir(version_path):
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
                return True
    else:
        os.mkdir(version_path)
        files_to_download = list(version_index.keys())

    # Download the files
    password = None
    total_files = len(files_to_download)

    if _box_auth and _box_token is None:
        _box_token = get_box_token()

    for data_file in files_to_download:

        file_index = version_index.get(data_file)
        server_hash = file_index.get("hash")
        file_url = file_index.get("url")

        file_path = os.path.join(version_path, data_file)
        file_number = files_to_download.index(data_file) + 1

        downloaded_path = download_file(file_url, file_path, server_hash, password=password, _box_token=_box_token, file_message=f"{dataset} v{version} data files", file_number=file_number, total_files=total_files)

        while downloaded_path == "wrong_password":
            if password is None:
                password = getpass.getpass(prompt=f'Password for {dataset} dataset: ') # We manually specify the prompt parameter so it shows up in Jupyter Notebooks
            else:
                password = getpass.getpass(prompt="Wrong password. Try again: ")
            print("\033[F", end='\r') # Use an ANSI escape sequence to move cursor back up to the beginning of the last line, so in the next line we can clear the password prompt
            print("\033[K", end='\r') # Use an ANSI escape sequence to print a blank line, to clear the password prompt

            downloaded_path = download_file(file_url, file_path, server_hash, password=password, _box_token=_box_token, file_message=f"{dataset} v{version} data files", file_number=file_number, total_files=total_files)

    return True

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
    download_file(index_url, index_path, server_index_hash, file_message=f"{dataset} index")

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

def download_file(url, path, server_hash, password=None, _box_token=None, file_message=None, file_number=None, total_files=None): 
    """Download a file from a given url to the specified location.

    Parameters:
    url (str): The direct download url for the file.
    path (str): The path to the file (not just the directory) to save the file to on the local machine.
    server_hash (str): The hash for the file, to check it against. If check fails, try download one more time, then throw an exception.
    password (str, optional): If the file is password protected, the password for it. Unneeded otherwise.
    _box_token (str, optional): The OAuth2 token for Box, if we're downloading a file that needs that. Default of None ignores that option.
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
            if _box_token is not None: # We are using Box OAuth2
                download_url = f"https://api.box.com/2.0/files/{url}/content/" # url is actually file ID
                headers = dict(HEADERS)
                headers["Authorization"] = f"Bearer {_box_token}"
                response = requests.get(download_url, headers=headers)

            elif password is None: # No password or OAuth2
                response = requests.get(url, headers=HEADERS, allow_redirects=True)

            else: # The file is password protected
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
        except requests.RequestException: # Parent class for all exceptions in the requests module
            raise NoInternetError("Insufficient internet. Check your internet connection.") from None
            
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

# Set up a localhost server to receive access token
app = Flask(__name__)

@app.route('/receive')
def receive():

    # Get the temporary access code
    code = request.args.get('code')

    # Create our "lock flag" file
    Path(LOCK_PATH).touch()

    # Save the code
    with open(CODE_FILE_PATH, "w") as code_file:
        code_file.write(code)

    # Remove lock flag
    os.remove(LOCK_PATH)

    # Shutdown the server. This will allow the thread it's running on to finish as well.
    request.environ.get("werkzeug.server.shutdown")()

    return "Authentication successful. You can close this window."

def get_box_token():

    # Don't show starting message from server
    cli.show_server_banner = lambda *_: None

    # Don't show logs from server
    log = logging.getLogger('werkzeug')
    log.disabled = True

    # Set up authentication parameters
    base_url = "https://account.box.com/api/oauth2/authorize"
    client_id = "kztczhjoq3oes38yywuyfp4t9tu11it8"
    client_secret = "a5xNE1qj4Z4H3BSJEDVfzbxtmxID6iKY"
    login_url = f"{base_url}?client_id={client_id}&response_type=code"

    # Start the server
    server = Thread(target=app.run, kwargs={"port": "8003"})
    server.start()

    # Send the user to the "Grant access" page
    webbrowser.open(login_url)
    login_msg = "Please login to Box on the webpage that was just opened and grant access for cptac to download files through your account. If you accidentally closed the browser window, press Ctrl+C and call the download function again."
    print(login_msg)

    # Get the temporary access code from the server on the child process
    temp_code = None
    while temp_code is None:
        if os.path.isfile(CODE_FILE_PATH) and not os.path.isfile(LOCK_PATH):
            with open(CODE_FILE_PATH) as code_file:
                temp_code = code_file.read()
            os.remove(CODE_FILE_PATH)
        else:
            time.sleep(1)

    # Wait for the server process to finish. Note in the receive endpoint handler shuts down the server automatically after receiving the code, so we don't need to worry about trying to kill it.
    server.join()

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
