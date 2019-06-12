import hashlib
import os
import requests

def sync(dataset, version="latest"):
    """Sync the specified version of the specified dataset.

    Parameters:
    dataset (str): The name of the dataset to sync
    version (str, optional): Which version of the dataset to sync. Defaults to latest.

    Returns:
    bool: Indicates whether sync was successful.
    """
    # Standardize parameter casing
    dataset = dataset.lower()
    version = version.lower()

    # Get the path to our data directory
    path_here = os.path.abspath(os.path.dirname(__file__))
    dataset_dir = f"data_{dataset}"
    dataset_path = os.path.join(path_here, dataset_dir)
    if not os.path.isdir(dataset_path):
        print(f"{dataset} is not a valid dataset.")
        return False

    updated = update_index(dataset_path)
    if not updated:
        print("Insufficient internet to sync. Check your internet connection.")
        return False

    # Load the index
    index = get_index(dataset_path)

    # If they chose to sync "latest", get the latest version number
    if version == "latest":
        version = get_latest_version_number(index, dataset_path)
        if version is None: # Latest version installed did not match latest in index. get_latest_version_number already printed error message.
            return False

    # Check that they chose a valid version
    if version not in index.keys():
        print(f"{version} is an invalid version for this dataset. Valid versions: {', '.join(index.keys())}")
        return False

    # If they haven't downloaded this version previously, make a new directory for it
    version_dir = f"{dataset}_v{version}"
    version_path = os.path.join(dataset_path, version_dir)
    if not os.path.isdir(version_path):
        os.mkdir(version_path)

    # Check the files in that version of the dataset. Download if don't exist, and update if have been changed.
    version_index = index.get(version)
    for data_file in version_index.keys():

        file_index = version_index.get(data_file)
        file_path = os.path.join(version_path, data_file)
        local_hash = hash_file(file_path)
        server_hash = file_index.get("hash")

        if local_hash != server_hash:
            file_url = file_index.get("url")
            downloaded_path = download_file(file_url, file_path, server_hash)
            if downloaded_path is None:
                print("Insufficient internet to sync. Check your internet connection.")
                return False
    return True

def get_latest_version_number(index, dataset_path):
    """Check that latest locally installed data is same version as latest version in index. If so, return latest version number; else, return None.

    Parameters:
    index (dict): The parsed index for the dataset.
    dataset_path (str): The path to the dataset's main directory

    Returns: 
    str: The latest version number.
    """
    index_latest = max(index.keys(), key=float) # See what the highest version in the index is
    latest_installed = get_latest_installed(dataset_path)
    if (index_latest == latest_installed) or (latest_installed is None):
        return index_latest
    else:
        print(f"Ambiguous request for latest version. Latest version in index is {index_latest}, but latest version installed locally is {latest_installed}. To download the latest version in the index, run cptac.sync with '{index_latest}' as the version parameter. To perform your requested action with the latest version that is installed locally, run your desired function with '{latest_installed}' as the version parameter.")
        return None

def update_index(dataset_path):
    """Check if the index of the given dataset is up to date with server version, and update it if needed.

    Parameters:
    dataset_path (str): The path to the dataset to check the index of.

    Returns:
    bool: Indicates if we were able to check the index and update if needed (i.e. we had internet)
    """
    # Define our file names we'll need
    index_urls_file = "index_urls.tsv"
    index_hash_file = "index_hash.txt"
    index_file = "index.txt"

    # Get, from the server, what the md5 hash of our index file should be
    index_urls_path = os.path.join(dataset_path, index_urls_file)
    urls_dict = parse_tsv_dict(index_urls_path)
    index_hash_url = urls_dict.get(index_hash_file)

    print(f"Checking that index is up-to-date...", end='\r')
    server_index_hash = download_text(index_hash_url).strip()
    print("\033[K", end='\r') # Erase the status message

    if server_index_hash is None: # I.e., we have no internet
        return False
    else:
        index_path = os.path.join(dataset_path, index_file)
        local_index_hash = hash_file(index_path)
        if local_index_hash == server_index_hash:
            return True
        else:
            index_url = urls_dict.get(index_file)
            index_downloaded_path = download_file(index_url, index_path, server_index_hash)
            if index_downloaded_path is None:
                return False
            else:
                return True

def get_latest_installed(dataset_path):
    """Return the latest version number installed in a dataset directory.

    Parameters:
    dataset_path (str): The path to the directory containing the data files.

    Returns:
    str: The latest version installed locally.
    """
    dirs = [dir.strip() for dir in os.listdir(dataset_path)
                if os.path.isdir(os.path.join(dataset_path, dir))]
    
    dataset_dir = dataset_path.split(os.sep)[-1]
    dataset_name = dataset_dir.split('_')[1] # For example, we get 'endometrial' from 'data_endometrial'
    version_dir_prefix = dataset_name + '_v'
    versions = [dir.replace(version_dir_prefix, '') for dir in dirs
                    if dir.startswith(version_dir_prefix)]
    if len(versions) == 0:
        return None
    latest_installed = max(versions, key=float)
    return latest_installed

def get_index(dataset_path):
    """Get the index for a dataset, as a nested dictionary

    Parameters:
    dataset_path (str): The path the dataset you want the index of.

    Returns:
    dict: The index, as a nested dictionary.
    """
    index_file = "index.txt"
    index_path = os.path.join(dataset_path, index_file)
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
            file_name = line_list[0]
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
        with open(path, 'rb') as file_obj:
            hash = hash_bytes(file_obj.read())
        return hash
    else:
        return None

def hash_bytes(bytes):
    """Hash the given bytes.

    Parameters:
    bytes (bytes): The bytes to has.

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
        response = requests.get(url, allow_redirects=True)
        response.raise_for_status() # Raises an HTTPError if the response code was unsuccessful
    except requests.RequestException:
        return None

    text = response.text
    return text

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
    if os.path.isfile(path):
        action = "Updating"
    else:
        action = "Downloading"
    print(f"{action} {file_name}...", end='\r')

    for i in range(2):
        try:
            response = requests.get(url, allow_redirects=True)
            response.raise_for_status() # Raises an HTTPError if the response code was unsuccessful
        except requests.RequestException:
            return None
            
        local_hash = hash_bytes(response.content)
        if local_hash == server_hash: # Only replace the old file if the new one downloaded successfully.
            with open(path, 'wb') as dest:
                dest.write(response.content)
            print("\033[K", end='\r') # Erase the downloading message
            return path
