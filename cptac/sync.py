import hashlib
import os
import requests
import glob

sync(dataset, version="latest"):
    """Sync the specified version of the specified dataset.

    Parameters:
    dataset (str): The name of the dataset to sync
    version (str, optional): Which version of the dataset to sync. Defaults to latest.

    Returns:
    int: Indicates whether sync was successful.
    """
    path_here = os.path.abspath(os.path.dirname(__file__))
    datatset_dir = "data_" + dataset
    dataset_path = os.path.join(path_here, dataset_dir)

    index_urls_path = os.path.join(dataset_path, "index_urls.tsv")
    with open(index_urls_path, 'r') as urls_file:
        urls_lines = urls_file.readlines()

    urls_dict = {}
    for line in urls_lines:
        line_list = line.strip().split("\t")
        key = line_list[0]
        value = line_list[1]
        urls_dict[key] = value

    index_hash_url = urls_dict.get("index_hash.txt")
    index_hash = download_text(index_hash_url)
    if index_hash is None: # Download failed
        # Print message
        return None

download_response(url):
    """Get the response from a GET request to the given url. Return None if internet has issues.

    Parameters:
    url (str): The url to send the GET request to.

    Returns:
    requests.Response: The response to the GET request.
    """
    try:
        response = requests.get(url, allow_redirects=True)
    #catch no internet

    return response

download_text(url):
    """Get the content of a response object from a GET request to a direct download url. Return None if internet has issues.

    Parameters:
    url (str): The direct download url.

    Returns:
    str: The response's content, as text.
    """
    response = download_response(url)
    if response is None:
        return None

    return response.text

download_file(url, path):
    """Download a file from a given url to the specified location.

    Parameters:
    url (str): The direct download url for the file.
    path (str): The path to the file (not just the directory) to save the file to on the local machine.

    Returns:
    str: The path the file was downloaded to.
    """
    response = download_response(url)
    if response is None:
        return None

    with open(path, 'wb') as dest:
        dest.write(response.content)

    return path
