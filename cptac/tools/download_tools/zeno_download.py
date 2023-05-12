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
import pandas as pd
from tqdm import tqdm
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor

import cptac
from cptac.exceptions import *

# Some websites don't like requests from sources without a user agent. Let's preempt that issue.
USER_AGENT = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0)'
HEADERS = {'User-Agent': USER_AGENT}
DATA_DIR = os.path.join(cptac.CPTAC_BASE_DIR, "data")
INDEX_FILE_NAME = 'all_index.txt'
INDEX_FILE_PATH = os.path.join(DATA_DIR, INDEX_FILE_NAME)
STATIC_DOI = '10.5281/zenodo.7897498'
RECORD_ID = STATIC_DOI.split('.')[-1]
ZENO_TOKEN = 'GijLB8joEFbeVEBQsjtJ8rH1uXMK8p5REgkNTfgHCMSR5LDyisZiZx1BRPQT'
AUTH_HEADER = {'Authorization': 'Bearer ' + ZENO_TOKEN}


def zeno_download(cancer: str, source: str, datatypes: str) -> bool:
    """
    Downloads data files for a specific cancer, source, and datatype from Zenodo

    :param cancer: The cancer type (e.g. 'brca').
    :param source: The data source (e.g. 'harmonized').
    :param datatype: The datatype of the files to download (e.g. 'clinical')
    """
    if not cancer or not source or not datatypes:
        raise InvalidParameterError("Cancer, source, and datatypes must be provided.")

    # FIXME: Change datatypes from list of strings to single string
    dtype = datatypes[0]

    try:
        bucket = get_bucket()
        # Get file name from index file
        if not os.path.exists(INDEX_FILE_PATH):
            get_data(f"{bucket}/{INDEX_FILE_NAME}", INDEX_FILE_NAME)
        with open(INDEX_FILE_PATH) as index:
            files = dict([line.split('\t') for line in index])
        file_name = files[f"{source}_{cancer}_{dtype}"].strip('\n')

        # Download requested dataframe
        output_dir = DATA_DIR + f"/data_{source}_{cancer}"  # FIXME: Requires version number
        os.makedirs(output_dir, exist_ok=True)
        output_file = file_name[len(f"{source}_{cancer}_"):]
        get_data(f"{bucket}/{file_name}", f"{output_dir}/{output_file}")

        return True

    except Exception as e:
        raise HttpResponseError(f"Failed to download data file for {source} {cancer} {dtype} with error:\n{e}") from e


def get_bucket() -> str:
    'Gets the bucket in zenodo that houses all data files.'
    projects = requests.get("https://zenodo.org/api/deposit/depositions", headers=AUTH_HEADER).json()
    for project in projects:
        if project['conceptrecid'] == RECORD_ID:
            return project['links']['bucket']
    raise CptacDevError("Failed to get bucket. Perhaps check that the token is correct?")

def download_chunk(url, start, end, file_path):
    headers = {'Range': f'bytes={start}-{end}'}
    headers.update(AUTH_HEADER)
    response = requests.get(url, headers=headers, stream=True)
    response.raise_for_status()

    with open(file_path, 'rb+') as data_file:
        data_file.seek(start)
        data_file.write(response.content)

def get_data(url: str, subfolder: str = '', num_threads: int = 4) -> str:
    if not os.path.exists(os.path.split(subfolder)[0]):
        os.makedirs(os.path.split(subfolder)[0], exist_ok=True)
    response = requests.head(url, headers=AUTH_HEADER)
    response.raise_for_status()

    file_size = int(response.headers['content-length'])
    chunk_size = file_size // num_threads

    # Create an empty file with the same size as the file to be downloaded
    with open(os.path.join(DATA_DIR, subfolder), 'wb') as data_file:
        data_file.truncate(file_size)

    try:
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = []
            for i in range(num_threads):
                start = i * chunk_size
                end = start + chunk_size - 1 if i != num_threads - 1 else file_size - 1
                futures.append(executor.submit(download_chunk, url, start, end, os.path.join(DATA_DIR, subfolder)))

            for future in tqdm(concurrent.futures.as_completed(futures),
                               desc=f"Downloading {os.path.split(subfolder)[1]}",
                               total=num_threads,
                               unit='chunk'):
                future.result()  # Raise any exception encountered during download

        return subfolder

    except:  # Yes, even keyboard interrupts
        os.remove(subfolder)
        raise



# def download_index_file_if_needed():
#     """
#     Downloads the index file if it does not already exist in the data directory

#     :return: The path to the index file
#     """
#     index_path = os.path.join(DATA_DIR, INDEX_FILE_NAME)
#     if not os.path.exists(index_path):
#         index_url = get_index_file_url()
#         if not index_url:
#             raise FileNotFoundError(f"Index file '{INDEX_FILE_NAME}' not found in Zenodo record (DOI: {STATIC_DOI})")
#         wget.download(index_url, index_path)

#     return index_path

# def get_index_file_url():
#     """
#     Gets the URL of the index file from the Zenodo record.

#     :return: The URL of the index file.
#     :raise FileNotFoundError: If the index file is not found in the Zenodo record
#     """
#     zenodo = zenodopy.Client()
#     record = zenodo.get_urls_from_doi(STATIC_DOI)

#     if not record:
#         return None

#     for url in record:
#         if url.endswith(INDEX_FILE_NAME):
#             return url

#     raise FileNotFoundError(f"Index file '{INDEX_FILE_NAME}' not found in Zenodo record (DOI: {STATIC_DOI})")

#     return None

# def get_file_names(cancer, source, datatypes, index_path):
#     """
#     Gets the URLs of the files that match the given cancer, source, and datatypes

#     :param cancer: The cancer type (e.g. 'brca)
#     :param source: The data source (e.g. 'harmonized)
#     :param datatype: The datatype of the files to download (e.g. 'clinical')
#     :param index_path: The path to the index file.
#     :return: A list of file URLs that match the given parameters
#     """
#     file_names = []

#     with open(index_path, 'r') as f:
#         for line in f:
#             tokens = line.strip().split('\t')
#             file_identifiers = tokens[0].split('_')

#             if file_identifiers[0] == source and file_identifiers[1] == cancer and file_identifiers[2] in datatypes:
#                 file_name = tokens[1]
#                 file_names.append(file_name)
#     return file_names

def download_text(url):
    """Download text from a direct download url for a text file.

    Parameters:
    url (str): The direct download url for the text.

    Returns:
    str: The downloaded text.
    """

    try:
        response = requests.get(url, headers=HEADERS, allow_redirects=True)
        response.raise_for_status()  # Raises a requests HTTPError if the response code was unsuccessful
    except requests.RequestException:  # Parent class for all exceptions in the requests module
        raise NoInternetError("Insufficient internet. Check your internet connection.") from None

    text = response.text.strip()
    return text

# def parse_tsv_dict(path):
#     """Read in a dictionary from the given two column tsv file.

#     Parameters:
#     path (str): The path to the two column tsv file.

#     Returns:
#     dict: The tsv file read into a dictionary.
#     """
#     if not os.path.isfile(path):
#         raise MissingFileError(f"Missing file {path}. Please update the cptac package to restore.")

#     with open(path, 'r') as data_file:
#         lines = data_file.readlines()

#     data_dict = {}
#     for line in lines:
#         if len(line.strip()) == 0:
#             continue
#         line_list = line.strip().split("\t")
#         key = line_list[0]
#         value = line_list[1]
#         data_dict[key] = value

#     return data_dict