import hashlib
import os
import requests
import glob
from .exceptions import InvalidVersionError

def sync(dataset, version="latest"):
    """Sync the specified version of the specified dataset.

    Parameters:
    dataset (str): The name of the dataset to sync
    version (str, optional): Which version of the dataset to sync. Defaults to latest.

    Returns:
    int: Indicates whether sync was successful.
    """
    # Get the path to our data directory
    path_here = os.path.abspath(os.path.dirname(__file__))
    datatset_dir = "data_" + dataset
    dataset_path = os.path.join(path_here, dataset_dir)

    # Define our file names we'll need
    index_urls_file = "index_urls.tsv"
    index_hash_file = "index_hash.txt"
    index_file = "index.txt"

    # Get, from the server, what the md5 hash of our index file should be
    index_urls_path = os.path.join(dataset_path, index_urls_file)
    urls_dict = parse_tsv_dict(index_urls_path)
    index_hash_url = urls_dict.get(index_hash_file)
    try:
        index_hash_response = requests.get(index_hash_url, allow_redirects=True)
    #catch no internet
    server_index_hash = index_hash_response.text

    # Check our index against the server hash, and update it if needed.
    index_path = os.path.join(dataset_path, index_file)
    local_index_hash = hash_file(index_path)
    if local_index_hash != server_index_hash:
        index_url = urls_dict.get(index_file)
        try:
            downloaded = download_file(url, index_path)
        # catch no internet

    # Read in the index
    index = parse_index(index_path)

    # If they chose to sync "latest", make sure the locally installed latest is the latest version
    if version == "latest":
        server_latest = max(index.keys(), key=float) # See what the highest version in the server index is
        local_latest = get_local_latest(dataset_path)

    # Check that they chose a valid version
    if version not in index.keys():
        raise InvalidVersionError(version, index.keys())

    # If they haven't downloaded this version previously, make a new directory for it

    # Check the files in that version of the dataset. Download if don't exist, and update if have been changed.

def get_local_latest(data_path):
    """Return the latest version number installed in a data directory.

    Parameters:
    data_path (str): The path to the directory containing the data files.

    Returns:
    str: The latest version installed locally.
    """
    dirs = [dir.strip() for dir in os.listdir(data_path)
                    if os.path.isdir(os.path.join(data_path, dir))]
    
    data_dir = data_path.split(os.sep)[-1]
    dataset_name = data_dir.split('_')[1] # For example, we get 'endometrial' from 'data_endometrial'
    version_dir_prefix = dataset_name + '_v'
    versions = [dir.replace(version_dir_prefix, '') for dir in dirs
                    if dir.startswith(version_dir_prefix)]
    local_latest = max(versions, key=float)
    return local_latest

def parse_index(index_path):
    """Read the given index file into a nested dictionary.

    Parameters:
    index_path (str): The path to the file containing the index

    Returns:
    dict: The index in a nested dictionary.
    """
    with open(index_path, 'r') as index_file:
        index_lines = index_file.readlines()

    index = {}
    version = None
    for line in index_lines:
        line = line.strip()
        if line.startswith('#'):
            version = int(line[1:])
            index[version] = {}
        else:
            line_list = line.split('\t')
            file_name = list_list[0]
            file_hash = line_list[1]
            file_url = line_list[2]
            index[version][file_name] = {}
            index[version][file_name]["hash"] = file_hash
            index[version][file_name]["url"] = file_url
    return index

def parse_tsv_dict(path):
    """Read in a dictionary from the given two column tsv file.

    Parameters:
    path (str): The path to the two column tsv file.

    Returns:
    dict: The tsv file read into a dictionary.
    """
    with open(path, 'r') as data_file:
        lines = data_file.readlines()

    data_dict = {}
    for line in lines:
        line_list = line.strip().split("\t")
        key = line_list[0]
        value = line_list[1]
        data_dict[key] = value

    return data_dict


def hash_file(path):
    """Return the md5 hash for the file at the given path. Return None if file doesn't exist.

    Parameters: 
    path (str): The absolute path to the file to hash.

    Returns:
    str: The hash for the file.
    """
    if os.path.isfile(path):
        hasher = hashlib.md5()
        with open(path, 'rb') as file_obj:
            buffer = file_obj.read()
            hasher.update(buffer)
        hash = hasher.hexdigest()
        return hash
    else:
        return None

def download_file(url, path, server_hash): 
    """Download a file from a given url to the specified location.

    Parameters:
    url (str): The direct download url for the file.
    path (str): The path to the file (not just the directory) to save the file to on the local machine.
    server_hash (str): The hash for the file, to check it against. If check fails, try download one more time, then throw an exception.

    Returns:
    str: The path the file was downloaded to.
    """
    file_name = path.split(os.sep)[-1]
    print("Downloading {}...".format(file_name), end='\r')

    counter = 0
    while counter is not None:
        if counter >= 2:
            # Throw exception

        response = requests.get(url, allow_redirects=True)
        with open(path, 'wb') as dest:
            dest.write(response.content)

        local_hash = hash_file(path)
        if local_hash == server_hash:
            counter = None
        else:
            counter += 1

    print("\033[K", end='\r') # Erase the downloading message
    return path
