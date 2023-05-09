import getpass
import glob
import os
import packaging.version
import requests
import warnings
import zenodopy
import wget
import requests
import urllib.parse

import cptac
from cptac.exceptions import *

# Some websites don't like requests from sources without a user agent. Let's preempt that issue.
USER_AGENT = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0)'
HEADERS = {'User-Agent': USER_AGENT}
DATA_DIR = os.path.join(cptac.CPTAC_BASE_DIR, "data")
INDEX_FILE_NAME = 'all_index.txt'
INDEX_FILE_PATH = os.path.join(DATA_DIR, INDEX_FILE_NAME)
STATIC_DOI = '10.5281/zenodo.7897498'
ZENO_TOKEN = 'GijLB8joEFbeVEBQsjtJ8rH1uXMK8p5REgkNTfgHCMSR5LDyisZiZx1BRPQT'

def zeno_download(cancer, source, datatypes):
    """
    Downloads data files for a specific cancer, source, and datatype from Zenodo

    :param cancer: The cancer type (e.g. 'brca').
    :param source: The data source (e.g. 'harmonized').
    :param datatype: The datatype of the files to download (e.g. 'clinical')
    """

    if not cancer or not source or not datatypes:
        raise ValueError("Cancer, source, and datatypes must be provided.")

    # Download the index file if it's not already present
    index_path = download_index_file_if_needed()

    # Get the file names matching the cancer, source, and datatypes
    file_names = get_file_names(cancer, source, datatypes, index_path)

    if not file_names:
        raise FileNotFoundError(
            f"No matching files found for source='{source}', cancer='{cancer}', datatype='{datatypes}'")

    # Create the output directory if it doesn't exist
    output_folder = os.path.join(DATA_DIR, f"data_{source}_{cancer}")

    zenodo = zenodopy.Client()
    record = zenodo.get_urls_from_doi(STATIC_DOI)

    # Download each file in file_names
    for file_name in file_names:
        url_file_name = urllib.parse.quote(file_name)
        for url in record:
            if url.endswith(url_file_name):
                destination_path = os.path.join(output_folder, file_name)
                try:
                    response = requests.get(url, headers=HEADERS, allow_redirects=True)
                    response.raise_for_status()
                    with open(destination_path, "wb") as f:
                        f.write(response.content)
                except requests.exceptions.HTTPError as e:
                    if e.response.status_code == 404:
                        print(f"File not found on server: {file_name}")
                    else:
                        raise
                break

    return True

def download_index_file_if_needed():
    """
    Downloads the index file if it does not already exist in the data directory
    
    :return: The path to the index file
    """
    index_path = os.path.join(DATA_DIR, INDEX_FILE_NAME)
    if not os.path.exists(index_path):
        index_url = get_index_file_url()
        if not index_url:
            raise FileNotFoundError(f"Index file '{INDEX_FILE_NAME}' not found in Zenodo record (DOI: {STATIC_DOI})")
        wget.download(index_url, index_path)

    return index_path

def get_index_file_url():
    """
    Gets the URL of the index file from the Zenodo record.

    :return: The URL of the index file.
    :raise FileNotFoundError: If the index file is not found in the Zenodo record
    """
    zenodo = zenodopy.Client()
    record = zenodo.get_urls_from_doi(STATIC_DOI)

    if not record:
        return None

    for url in record:
        if url.endswith(INDEX_FILE_NAME):
            return url

    raise FileNotFoundError(f"Index file '{INDEX_FILE_NAME}' not found in Zenodo record (DOI: {STATIC_DOI})")

    return None

def get_file_names(cancer, source, datatypes, index_path):
    """
    Gets the URLs of the files that match the given cancer, source, and datatypes

    :param cancer: The cancer type (e.g. 'brca)
    :param source: The data source (e.g. 'harmonized)
    :param datatype: The datatype of the files to download (e.g. 'clinical')
    :param index_path: The path to the index file.
    :return: A list of file URLs that match the given parameters
    """
    file_names = []

    with open(index_path, 'r') as f:
        for line in f:
            tokens = line.strip().split('\t')
            file_identifiers = tokens[0].split('_')

            if file_identifiers[0] == source and file_identifiers[1] == cancer and file_identifiers[2] in datatypes:
                file_name = tokens[1]
                file_names.append(file_name)
    return file_names

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