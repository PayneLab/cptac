import os
import requests
from tqdm import tqdm
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor

import cptac
from cptac.exceptions import *

# Some websites don't like requests from sources without a user agent. Let's preempt that issue.
# This variable sets a user agent to be included in the request headers.
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

    dtype = datatypes[0] # FIXME: Change datatypes from list of strings to single string

    try:
        bucket = get_bucket()
        # Get file name from index file
        if not os.path.exists(INDEX_FILE_PATH):
            get_data(f"{bucket}/{INDEX_FILE_NAME}", INDEX_FILE_NAME)
        with open(INDEX_FILE_PATH) as index:
            files = dict([line.split('\t') for line in index])

        # Download requested dataframe
        if source in ['harmonized', 'mssm']:
            output_dir = DATA_DIR + f"/data_{source}"
            file_name = files[f"{source}_{dtype}"].strip('\n')
            output_file = file_name[len(f"{source}_"):]
        else:
            output_dir = DATA_DIR + f"/data_{source}_{cancer}"  # FIXME: Requires version number
            file_name = files[f"{source}_{cancer}_{dtype}"].strip('\n')
            output_file = file_name[len(f"{source}_{cancer}_"):]
        os.makedirs(output_dir, exist_ok=True)
        get_data(f"{bucket}/{file_name}", f"{output_dir}/{output_file}")

        return True

    except Exception as e:
        raise HttpResponseError(f"Failed to download data file for {source} {cancer} {dtype} with error:\n{e}") from e


def get_bucket() -> str:
    """
    Gets the bucket in Zenodo that houses all data files.

    :return: The URL of the Zenodo bucket containing the data files.
    """
    projects = requests.get("https://zenodo.org/api/deposit/depositions", headers=AUTH_HEADER).json()
    for project in projects:
        if project['conceptrecid'] == RECORD_ID:
            return project['links']['bucket']
    raise CptacDevError("Failed to get bucket. Perhaps check that the token is correct?")


def download_chunk(url, start, end, file_path):
    """
    Downloads a chunk of a file and writes it to the specified file path.
    
    :param url: THE URL of the file to download.
    :param start: The start byte of the chunk to download.
    :param end: The end byte of the chunk to download.
    :param file_path: The path to the file where the downloaded chunk should be written
    """
    headers = {'Range': f'bytes={start}-{end}'}
    headers.update(AUTH_HEADER)
    response = requests.get(url, headers=headers, stream=True)
    response.raise_for_status()

    with open(file_path, 'rb+') as data_file:
        data_file.seek(start)
        data_file.write(response.content)


def get_data(url: str, subfolder: str = '', num_threads: int = 4) -> str:
    """
    Downloads a file using multithreading and saves it to the specified subfolder.

    :param url: The URL of the file to download.
    :param subfolder: The subfolder where the downloaded file should be saved.
    :param num_threads: The number of threads to use for downloading the file (default is 4).
    :return: The path of the downloaded file.
    """
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


def download_text(url):
    """
    Download text from a direct download url for a text file.

    :param url: The direct download url for the text.
    :return: The downloaded text
    """

    try:
        response = requests.get(url, headers=HEADERS, allow_redirects=True)
        response.raise_for_status()  # Raises a requests HTTPError if the response code was unsuccessful
    except requests.RequestException:  # Parent class for all exceptions in the requests module
        raise NoInternetError("Insufficient internet. Check your internet connection.") from None

    text = response.text.strip()
    return text