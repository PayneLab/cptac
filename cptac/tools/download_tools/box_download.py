import bs4
import getpass
import glob
import hashlib
import logging
import os
import packaging.version
import pandas as pd
import requests
import threading
import warnings
import webbrowser

from datetime import datetime, timedelta
from importlib.resources import path
from queue import Queue
from werkzeug import Request, Response
from werkzeug.serving import make_server

import cptac
from cptac.exceptions import *

# Some websites don't like requests from sources without a user agent. Let's preempt that issue.
USER_AGENT = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0)'
HEADERS = {'User-Agent': USER_AGENT}
 
def box_download(cancer, source, datatypes, version, redownload):
   
    if source in ["harmonized", "mssm"]:
        dataset = source
    else: 
        dataset = source + "_" + cancer

    # Get our dataset path
    dataset_path = get_data_path(dataset)

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
    files_to_download = gather_files(version_path=version_path, version_index=version_index, redownload=redownload)

    # Return true if no new files to download
    if files_to_download is None: 
        return True

    # Else Download the files
    password = cptac.box_auth.get_password(dataset)
    total_files = len(files_to_download)

    for data_file in files_to_download:

        file_index = version_index.get(data_file)
        server_hash = file_index.get("hash")
        file_url = file_index.get("url")

        file_path = os.path.join(version_path, data_file)
        file_number = files_to_download.index(data_file) + 1

        downloaded_path = download_file(url=file_url, path=file_path, server_hash=server_hash,source=source, password=password, file_message=f"{dataset} v{version} data files", file_number=file_number, total_files=total_files)

        while downloaded_path == "wrong_password":
            if password is None:
                password = getpass.getpass(prompt=f'Password for {dataset} dataset: ') # We manually specify the prompt parameter so it shows up in Jupyter Notebooks
            else:
                password = getpass.getpass(prompt="Wrong password. Try again: ")
            print("\033[F", end='\r') # Use an ANSI escape sequence to move cursor back up to the beginning of the last line, so in the next line we can clear the password prompt
            print("\033[K", end='\r') # Use an ANSI escape sequence to print a blank line, to clear the password prompt

            downloaded_path = download_file(url=file_url, path=file_path, source=source, server_hash=server_hash, password=password, file_message=f"{dataset} v{version} data files", file_number=file_number, total_files=total_files)

    return True

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

    # download the files
    for i in range(2):
        try:
            # check if the required file is from a source whose files are stored on Box.com
            if source in ["bcm", "broad", "harmonized", "mssm", "pdc", "umich", "washu"]: # We are using Box OAuth2
                cptac.box_auth.refresh_token() # global box_auth object
                download_url = f"https://api.box.com/2.0/files/{url}/content" # url is actually file ID
                headers = dict(HEADERS)
                headers["Authorization"] = f"Bearer {cptac.box_auth.get_box_token()}"
                response = requests.get(download_url, headers=headers)
            
            elif password is None: # No password or OAuth2 (awg files and index files)
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
    raise DownloadFailedError(f"Download failed for {file_name}. \nL_Hash: {local_hash}\nS_Hash: {server_hash}")

def get_data_path(dataset):
    """Get the path to the main directory for a dataset.

    Parameters:
    dataset (str): The path to get the directory for. Must be all lowercase.

    Returns:
    str: The path to the main directory of the specified dataset.
    """
    dataset_dir = f"data/data_{dataset}"
    dataset_path = os.path.join(cptac.CPTAC_BASE_DIR, dataset_dir)

    if os.path.isdir(dataset_path):
        return dataset_path
    else:
        source, cancer = dataset.split("_")
        raise InvalidParameterError(f"{source} is not a valid source for {cancer}. Path: {dataset_path}")

def validate_version(version, dataset, use_context, valid_versions=None):
    """Parse and validate a given version number. If version is "latest", check that index and installed latest match.

    Parameters:
    version (str): The version number to validate.
    dataset (str): The name of the dataset we're validating the version for.
    use_context (str): Either "download" or "init", depending on whether the function is being called as part of a data download or a dataset loading. Allows for more detailed error messages. Pass None if you don't want more detailed error messages.
    valid_versions (list of str, optional): A list of versions that are valid for this dataset. This way, we can ensure that if there's a new data version from the downloaded index that this version of the package doesn't have the code to handle, it won't try to load it. Default of None will skip this check.

    Returns:
    str: The version number.
    """
    # Get our dataset path, then our dataset index
    dataset_path = get_data_path(dataset)
    index = get_index(dataset)

    # See what the highest version in the index is
    index_latest = max(index.keys(), key=packaging.version.parse)

    # Parse and validate the version they passed
    if version in index.keys():
        if version != index_latest: # Print a warning if they're using an old version
            warnings.warn(f"Old {dataset} data version. Latest is {index_latest}. This is {version}.", OldDataVersionWarning, stacklevel=4)
        return_version = version

    elif version.lower() == "latest":
        latest_installed = get_latest_installed(dataset_path)
        if (index_latest == latest_installed) or (latest_installed is None):
            return_version = index_latest
        else: # If their latest installed version is different from the latest version recorded in the index, then we don't know which one they meant when they passed "latest".
            if use_context == "download":
                warnings.warn(f"Downloading new version of {dataset} dataset: {index_latest}. This will now be the default version when the dataset is loaded. If you wish to load an older version of the data, you must specify it with the 'version' parameter when you load the dataset.", DownloadingNewLatestWarning, stacklevel=3)
                return_version = index_latest
            elif use_context == "init":
                raise AmbiguousLatestError(f"You requested to load the {dataset} dataset. Latest version is {index_latest}, which is not installed locally. To install it, call the download function (either 'cptac.download' or 'cptac.pancan.download', depending on which module you're using) and pass '{dataset}' to the 'dataset' parameter. This will download the latest version, and you will then be able to load that version of the dataset. To skip this and instead load the older version that is already installed, re-call the function that generated this error, but pass '{latest_installed}' to the 'version' parameter.")
    else:
        raise InvalidParameterError(f"{version} is an invalid version for the {dataset} dataset. Valid versions: {', '.join(index.keys())}")

    if valid_versions is not None:
        if return_version not in valid_versions:
            raise PackageCannotHandleDataVersionError(f"You tried to load data version {return_version}, but your version of cptac can only handle these versions: {valid_versions}. Update your package to be able to load the new data. Or, if you cannot currently update, manually specify the old data version using the 'version' parameter when you load the dataset.")

    return return_version

def get_latest_installed(dataset_path):
    """Return the latest version number installed in a dataset directory.

    Parameters:
    dataset_path (str): The path to the dataset of interest.

    Returns:
    str: The latest version installed locally. Returns None if no versions are installed.
    """
    dirs = [dir.strip() for dir in os.listdir(dataset_path)
                if os.path.isdir(os.path.join(dataset_path, dir))]

    dataset_dir = dataset_path.split(os.sep)[-1]
    dataset_name = dataset_dir.split('_')[1] # For example, we get 'endometrial' from 'data_endometrial'
    version_dir_prefix = dataset_name + '_v'
    versions = [dir.replace(version_dir_prefix, '') for dir in dirs if dir.startswith(version_dir_prefix)]
    if len(versions) == 0:
        return
        
    latest_installed = max(versions, key=packaging.version.parse)
    return latest_installed

def get_index(dataset):
    """Get the index for a dataset, as a nested dictionary

    Parameters:
    dataset(str): The name of dataset you want the index of.

    Returns:
    dict: The index, as a nested dictionary.
    """
    dataset_path = get_data_path(dataset)
    index_file = "index.txt"
    index_path = os.path.join(dataset_path, index_file)

    # Check that the index is installed
    if not os.path.isfile(index_path):
        dataset_version_pattern = f"{dataset}_v*" # If not, check whether we've installed any version directories, to know what type of error to raise
        dataset_version_search = os.path.join(dataset_path, dataset_version_pattern)
        version_dirs = glob.glob(dataset_version_search)
        if len(version_dirs) > 0:  
            raise MissingFileError(f"Missing file '{index_file}'. Call the download function (either 'cptac.download' or 'cptac.pancan.download', depending on which module you're using) to download it, passing '{dataset}' to the 'dataset' parameter, and passing 'True' to the 'redownload' parameter.")
        else:
            raise DatasetNotInstalledError(f"{dataset} dataset is not installed. To install, call the download function (either 'cptac.download' or 'cptac.pancan.download', depending on which module you're using), passing '{dataset}' to the 'dataset' parameter.")

    with open(index_path, 'r') as index_file:
        index_lines = index_file.readlines()

    index = {}
    version = None
    for line in index_lines:
        line = line.strip()
        if line.startswith('#'):
            version = line[1:]
            index[version] = {}
        else:
            line_list = line.split('\t')
            if len(line_list) > 1:
                file_name = line_list[0]
                file_hash = line_list[1]
                file_url = line_list[2]
                file_datatype = line_list[3]
                index[version][file_name] = {}
                index[version][file_name]["hash"] = file_hash
                index[version][file_name]["url"] = file_url
                index[version][file_name]["datatype"] = file_datatype
    return index

def get_filtered_version_index(version_index, source, datatypes, version):
    """Filter the version index to only include files of the desired types.
    
    Also checks for invalid or unavailable datatypes.

    Parameters:
    version_index (dict): The version index dictionary
    datatypes (list of str): The datatypes desired.
    source: The source from where the datatypes will be loaded 

    Returns:
    dict: The filtered version index
    """
    found_datatypes = list()
    filtered_version_index = dict()
    annotation_types = ['clinical', 'derived_molecular', 'experimental_design', 'medical_history']

    # always add mapping and definition files
    helper_types = ['mapping', 'definitions']
    datatypes.extend(helper_types)
    datatypes = [d.lower() for d in datatypes]
    

    # add the files of the desired datatypes to a new dict (filtered_version_index)
    for file_name in version_index.keys():
        
        file_type = version_index[file_name]['datatype'].lower()
        
        # add desired annotation file types
        # Developers, note that this functionality depends very much on accurately specifying datatypes in the index files
        if file_type == 'annotation' and any(item in annotation_types for item in datatypes):
            filtered_version_index[file_name] = version_index[file_name]
            desired_annotations = set(annotation_types).intersection(set(datatypes))
            for ft in desired_annotations:
                found_datatypes.append(ft)

        # add mapping and definition files quietly
        elif file_type in helper_types:
            filtered_version_index[file_name] = version_index[file_name]

        # add file if of a desired types 
        elif file_type in datatypes:
            filtered_version_index[file_name] = version_index[file_name]
            found_datatypes.append(version_index[file_name]['datatype'])
    
    # check for invalid and unavailable datatypes
    found_datatypes = set(found_datatypes)
    datatypes = set(datatypes) - set(helper_types)
    if set([x.lower() for x in found_datatypes]) != set([d.lower() for d in datatypes]):
        w = None
        
        if source in ['harmonized', 'mssm']:
            w = f"These {source} datatypes were not found in version v_{version}: {set(datatypes) - set(found_datatypes)}\nSee cptac.list_datasets() for more info."
        
        else:
            source, cancer = source.split("_")
            w = f"These {cancer} datatypes were not found for the {source} source in version v_{version}: {set(datatypes) - set(found_datatypes)}\nSee cptac.list_datasets() for more info."
        
        warnings.warn(w, DataTypeNotInSourceWarning, stacklevel=2)

    # return version index with only files of desired types
    return filtered_version_index

def update_index(dataset):
    """Check if the index of the given dataset is up to date with server version, and update it if needed.

    Parameters:
    dataset (str): The name of the dataset to check the index of.

    Returns:
    bool: Indicates if we were able to check the index and update if needed (i.e. we had internet)
    """
    # Get the path to our dataset
    dataset_path = get_data_path(dataset)

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

def parse_tsv_dict(path):
    """Read in a dictionary from the given two column tsv file.

    Parameters:
    path (str): The path to the two column tsv file.

    Returns:
    dict: The tsv file read into a dictionary.
    """
    if not os.path.isfile(path): 
        raise MissingFileError(f"Missing file {path}. Please update the cptac package to restore.")

    with open(path, 'r') as data_file:
        lines = data_file.readlines()

    data_dict = {}
    for line in lines:                         
        if len(line.strip()) == 0:
            continue
        line_list = line.strip().split("\t")
        key = line_list[0]
        value = line_list[1]
        data_dict[key] = value

    return data_dict

def hash_file(path):
    """Return the md5 hash for the file at the given path.

    Parameters:
    path (str): The absolute path to the file to hash.

    Returns:
    str: The hash for the file.
    """
    with open(path, 'rb') as file_obj:
        hash = hash_bytes(file_obj.read())
    return hash

def hash_bytes(bytes):
    """Hash the given bytes.

    Parameters:
    bytes (bytes): The bytes to hash.

    Returns:
    str: The hash for the bytes.
    """
    hasher = hashlib.md5()
    hasher.update(bytes)
    hash = hasher.hexdigest()
    return hash

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