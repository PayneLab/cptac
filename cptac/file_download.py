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

import os
import requests
import getpass
import bs4
from .file_tools import *

def download(dataset, version="latest", redownload=False):
    """Download data files for the specified datasets. Defaults to downloading latest version on server.

    Parameters:
    dataset (str): The name of the dataset to download data for
    version (str, optional): Which version of the data files to download. Defaults to latest on server.
    redownload (bool, optional): Whether to redownload the data files, even if that version of the data is already downloaded. Default False.

    Returns:
    bool: Indicates whether download was successful.
    """
    # Get our dataset path
    dataset = dataset.lower()
    dataset_path = get_dataset_path(dataset)
    if dataset_path is None: # Invalid dataset. get_dataset_path already printed an error message.
        return False

    # Update the index
    updated = update_index(dataset)
    if not updated:
        print("Insufficient internet. Check your internet connection.")
        return False

    # Load the index
    index = get_index(dataset_path)

    # Validate the version number, including parsing if it's "latest"
    version = validate_version(version, dataset, use_context="download")
    if version is None: # Invalid version
        return False

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
                file_index = version_index.get(data_file)
                server_hash = file_index.get("hash")
                file_path = os.path.join(version_path, data_file)
                local_hash = hash_file(file_path) 

                if local_hash != server_hash:
                    files_to_download.append(data_file)

            if len(files_to_download) == 0:
                print("All files already downloaded and correct.")
                return True
    else:
        os.mkdir(version_path)
        files_to_download = list(version_index.keys())

    # Download the files
    password_protected_datasets = [
        ]
    password = None

    total_files = len(files_to_download)
    for data_file in files_to_download:

        if (dataset in password_protected_datasets) and (password is None):
            password = getpass.getpass()
            print("\033[F", end='\r') # Use an ANSI escape sequence to move cursor back up to the beginning of the last line, so in the next line we can clear the password prompt
            print("\033[K", end='\r') # Use an ANSI escape sequence to print a blank line, to clear the password prompt

        file_index = version_index.get(data_file)
        server_hash = file_index.get("hash")
        file_url = file_index.get("url")

        file_path = os.path.join(version_path, data_file)
        file_number = files_to_download.index(data_file) + 1

        downloaded_path = download_file(file_url, file_path, server_hash, password=password, file_message="data files", file_number=file_number, total_files=total_files)

        while downloaded_path == "wrong_password":
            password = getpass.getpass(prompt="Wrong password. Try again: ")
            print("\033[F", end='\r') # Use an ANSI escape sequence to move cursor back up to the beginning of the last line, so in the next line we can clear the password prompt
            print("\033[K", end='\r') # Use an ANSI escape sequence to print a blank line, to clear the password prompt

            downloaded_path = download_file(file_url, file_path, server_hash, password=password, file_message="data files", file_number=file_number, total_files=total_files)

        if downloaded_path is None:
            print("Insufficient internet. Check your internet connection.")
            return False

    formatted_name = "RenalCcrcc" if dataset == "renalccrcc" else dataset.title()
    print(f"{formatted_name} data download successful.")
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

    checking_msg = "Checking that index is up-to-date..."
    print(checking_msg, end='\r')
    server_index_hash = download_text(index_hash_url)
    print(" " * len(checking_msg), end='\r') # Erase the checking message

    if server_index_hash is None: # I.e., we have no internet
        return False
    else:
        index_path = os.path.join(dataset_path, index_file)
        local_index_hash = hash_file(index_path)
        if local_index_hash == server_index_hash:
            return True
        else:
            index_url = urls_dict.get(index_file)
            index_downloaded_path = download_file(index_url, index_path, server_index_hash, file_message="index")
            if index_downloaded_path is None:
                return False
            else:
                return True

def download_text(url):
    """Download text from a direct download url for a text file.

    Parameters:
    url (str): The direct download url for the text.

    Returns:
    str: The downloaded text.
    """
    try:
        response = requests.get(url, allow_redirects=True)
        response.raise_for_status() # Raises a requests HTTPError if the response code was unsuccessful
    except requests.RequestException: # Parent class for all exceptions in the requests module
        return

    text = response.text.strip()
    return text

def download_file(url, path, server_hash, password=None, file_message=None, file_number=None, total_files=None): 
    """Download a file from a given url to the specified location.

    Parameters:
    url (str): The direct download url for the file.
    path (str): The path to the file (not just the directory) to save the file to on the local machine.
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
            if password is None:
                response = requests.get(url, allow_redirects=True)
            else: # The file is password protected
                with requests.Session() as session: # Use a session object to save cookies
                    # Construct the urls for our GET and POST requests
                    get_url = url
                    post_url = get_url.replace("https://byu.box.com/shared", "https://byu.app.box.com/public")

                    # Send initial GET request and parse the request token out of the response
                    get_response = session.get(get_url) 
                    soup = bs4.BeautifulSoup(get_response.text, "html.parser")
                    token_tag = soup.find(id="request_token")
                    token = token_tag.get("value")

                    # Send a POST request, with the password and token, to get the data
                    payload = {
                        'password': password,
                        'request_token': token}
                    response = session.post(post_url, data=payload)

            response.raise_for_status() # Raises a requests.HTTPError if the response code was unsuccessful
        except requests.RequestException: # Parent class for all exceptions in the requests module
            return
            
        local_hash = hash_bytes(response.content)
        if local_hash == server_hash: # Only replace the old file if the new one downloaded successfully.
            with open(path, 'wb') as dest:
                dest.write(response.content)
            print(" " * len(download_msg), end='\r') # Erase the downloading message
            return path
        elif response.text.strip().startswith("<!DOCTYPE html>"): # The password was wrong, so we just got a webpage
            print(" " * len(download_msg), end='\r') # Erase the downloading message
            return "wrong_password"
