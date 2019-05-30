import hashlib
import os
import wget

def check_data(data_path):
    """Checks that all data files that should be in the given directory exist and are up to date. Downloads or re-downloads files as needed.

    Parameters:
    data_path (str): The path to the data directory for the given dataset.

    Returns: None
    """
    # Get our file with the urls for each of the data files, and read it into a dict
    urls_path = os.path.join(data_path, "urls", "urls.tsv") 
    urls_dict = load_dict(urls_path) # File names are keys, urls are values

    # For each file in the url dict, put the local hash in a dict. If it doesn't exist, put None in the dict.
    local_hashes = {}
    for data_file in urls_dict:
        if data_file == "checksums":
            continue # We haven't downloaded the checksums yet, so we can't hash them. And we re-download them every time, so we don't need to hash them anyways.
        data_file_path = os.path.join(data_path, data_file)
        if os.path.isfile(data_file_path):
            hasher = hashlib.md5()
            with open(data_file_path, 'rb') as data_file_obj:
                buffer = data_file_obj.read()
                hasher.update(buffer)
            file_hash = hasher.hexdigest()
            local_hashes[data_file] = file_hash
        else:
            local_hashes[data_file] = None

    # Download server checksums file into temp file, read into a dict, and delete temp file
    server_checksums_url = urls_dict["checksums"]
    server_checksums_path = os.path.join(data_path, "server_checksums.tsv.tmp")
    wget.download(server_checksums_url, server_checksums_path)
    server_checksums = load_dict(server_checksums_path) # File names are keys, checksums are values    
    os.remove(server_checksums_path)

    # Download or re-download files as needed
    for name, server_checksum in server_checksums.items():
        if name not in local_hashes.keys():
            print("This dataset has been updated to include a new data file: {0}. Please update the cptac package in order access the {0} file.".format(name))
            continue

        local_checksum = local_hashes[name]
        if local_checksum != server_checksum: # Either the local checksum is None (the data file hasn't yet been downloaded onto this computer), or the checksums don't match (the server version has been updated)
            file_url = urls_dict[name]
            local_file_path = os.path.join(data_path, name)
            wget.download(file_url, local_file_path)

            # Tell them we changed something
            if local_checksum is None:
                print("{} data file downloaded.".format(name))
            else:
                print("{} data file updated to new version.".format(name))

def load_dict(path):
    """Read a given tsv file into a dict.

    Parameters:
    path (str): The absolute path to the file to read.

    Returns:
    dict: The information in the file, read into a dict.
    """
    # Get the file
    with open(path, 'r') as data_file:
        data_lines = data_file.readlines()

    # Read it into a dict
    data_dict = {}
    for line in data_lines:
        line_list = line.strip().split("\t")
        key = line_list[0]
        value = line_list[1]
        data_dict[key] = value

    return data_dict
