import hashlib
import os
import wget
import glob

def check_data(data_path):
    """Checks that all data files that should be in the given directory exist and are up to date. Downloads or re-downloads files as needed.

    Parameters:
    data_path (str): The path to the data directory for the given dataset.

    Returns: None
    """
    print("Checking that data files are up-to-date...")

    # Get our file with the urls for each of the data files, and read it into a dict
    urls_path = os.path.join(data_path, "urls", "urls.tsv") 
    urls_dict = load_dict(urls_path) # File names are keys, urls are values
    checksums_file_name = "server_checksums.tsv" # This is one of the values in our dict.

    # For each file in the url dict, put the local hash in a dict. If it doesn't exist, put None in the dict.
    local_hashes = {}
    for data_file in urls_dict.keys():
        if data_file == checksums_file_name:
            continue # We haven't downloaded the checksums yet, so we can't hash them. And we re-download them every time, so we don't need to hash them anyways.
        data_file_path = os.path.join(data_path, data_file)
        if os.path.isfile(data_file_path):
            hasher = hashlib.md5()
            with open(data_file_path, 'rb') as data_file_obj:
                buffer = data_file_obj.read()
                hasher.update(buffer)
            file_hash = hasher.hexdigest()
            local_hashes[data_file] = file_hash
        else: # The file doesn't exist on our local machine yet--we need to download a new version.
            local_hashes[data_file] = None

    # Download server checksums file into temp file, read them into a dict, and delete temp file
    server_checksums_path = download_from_urls_dict(checksums_file_name, urls_dict, data_path, print_download_msg=False)
    server_checksums = load_dict(server_checksums_path) # File names are keys, checksums are values    
    os.remove(server_checksums_path) # Remove the file once we have the data

    # Download or re-download files as needed
    for name, server_checksum in server_checksums.items():
        if name not in local_hashes.keys():
            print("\n ******NEW DATA******\nThis dataset has been updated to include a new data file: {0}. Please update the cptac package in order access the {0} file.\n".format(name))
            continue

        local_checksum = local_hashes[name]
        if local_checksum is None: # The file didn't previously exist on the local machine
            downloaded_path = download_from_urls_dict(name, urls_dict, data_path)

        elif local_checksum != server_checksum: # The file is on the local machine, but out of date.
            update_response = input("File '{}' is out-of-date. Would you like to update it (y/n)? ".format(os.path.join(data_path, name)))
            valid_input = False
            while not valid_input:
                if update_response == 'y':
                    valid_input = True
                    downloaded_path = download_from_urls_dict(name, urls_dict, data_path)
                    print("{} updated to most current version.\n**********".format(downloaded_path))
                elif update_response == 'n':
                    valid_input = True
                    print("\n ******WARNING******\nFile '{}' not updated. We recommend updating it as soon as possible, to have the most current data.\n".format(os.path.join(data_path, name)))
                else: 
                    update_response = input("Invalid response. Please enter 'y' for yes or 'n' for no: ")

    print("Data check complete.")

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

def download_from_urls_dict(file_name, urls_dict, dir_path, print_download_msg=True):
    """Download a file from a url and save it to the given location.

    Parameters:
    file_name (str): The name of the desired file.
    urls_dict (dict): A dict containing the desired file's name as a key, and its url as the corresponding value.
    dir_path (str): The path to the directory to download the file to.
    print_download_msg (bool): Whether to print a message telling the user which file we're downloading. (Note: Even if False, wget will still print the download status bar).

    Returns:
    str: The path the file was downloaded to.
    """
    file_url = urls_dict[file_name]
    local_file_path = os.path.join(dir_path, file_name)

    if print_download_msg:
        print("Downloading {}...".format(file_name))
    try:
        downloaded_path = wget.download(file_url, local_file_path)
        if downloaded_path != local_file_path: # This happens when we replace a file--wget appends " (1)" to the file name
            if os.path.isfile(local_file_path):
                os.remove(local_file_path) # Remove the old file, if it exists
            os.rename(downloaded_path, local_file_path) # Rename the new file to the desired name
            downloaded_path = local_file_path # Update downloaded_path to the new location of the downloaded file.

    finally: # Delete any *.tmp files, in case a download was interrupted.
        tmp_files_path = os.path.join(dir_path, "*.tmp")
        tmp_files = glob.glob(tmp_files_path)
        for tmp_file in tmp_files:
            os.remove(tmp_file)
            print("\nTemporary file '{}' removed during clean-up. A file download may have been interrupted.".format(tmp_file))
    print() # Add a newline after wget's download status bar

    return downloaded_path
