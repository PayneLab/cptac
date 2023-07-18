import os
import requests
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
from hashlib import md5
from tqdm import tqdm
import cptac
from cptac.exceptions import *

# Set directory constants
DATA_DIR = os.path.join(cptac.CPTAC_BASE_DIR, "data")
STATIC_DOI = '10.5281/zenodo.7897498'
RECORD_ID = STATIC_DOI.split('.')[-1]
ZENO_TOKEN = 'GijLB8joEFbeVEBQsjtJ8rH1uXMK8p5REgkNTfgHCMSR5LDyisZiZx1BRPQT'
AUTH_HEADER = {'Authorization': 'Bearer ' + ZENO_TOKEN}
BUCKET=None


def fetch_repo_data() -> dict:
    """Fetches the repo data from Zenodo, including metadata and file links."""
    repo_link = f"https://zenodo.org/api/records/{RECORD_ID}"
    response = requests.get(repo_link, headers=AUTH_HEADER)
    response.raise_for_status()
    return response.json()


def init_files() -> None:
    """
    Initializes several files that are essential for cptac to run, such as the file index.
    Creates directory if it does not exist, fetches data from the repository, and writes
    index file.
    """
    os.makedirs(DATA_DIR, exist_ok=True)
    index_path = os.path.join(DATA_DIR, 'index.tsv')
    acetyl_mapping_path = os.path.join(DATA_DIR, 'cptac_genes.csv')
    brca_mapping_path = os.path.join(DATA_DIR, 'brca_mapping.csv')

    # Error handling for different exceptions
    try:
        repo_data = fetch_repo_data()
        global BUCKET
        BUCKET = repo_data['links']['bucket']
        index_data = []
        index_data.append(f"description\tfilename\tchecksum")
        for data_file in repo_data['files']:
            if data_file['key'].startswith('.'):
                continue # ignore hidden files
            filename_list = data_file['key'].split('-')
            description = '-'.join(filename_list[:3])
            index_data.append(f"{description}\t{'-'.join(filename_list)}\t{data_file['checksum']}")
        with open(index_path, 'w') as index_file:
            index_file.write('\n'.join(index_data))
        # Download some other necessart files
        if not os.path.isfile(acetyl_mapping_path):
            get_data(f"{BUCKET}/{'cptac_genes.csv'}", acetyl_mapping_path)
        if not os.path.isfile(brca_mapping_path):
            get_data(f"{BUCKET}/{'brca_mapping.csv'}", brca_mapping_path)
    except requests.ConnectionError:
        raise NoInternetError("Cannot initialize data files: No internet connection.")
    except requests.RequestException as e:
        raise HttpResponseError(f"Requesting data failed with the following error: {e}")
    except Exception as e:
        raise CptacError(f"ERROR: {e}")
    

def download(cancer: str, source: str, dtype: str, data_file: str) -> bool:
    """
    Downloads data files for a specific cancer, source, datatype, and file name from Zenodo

    :param cancer: The cancer type (e.g. 'brca').
    :param source: The data source (e.g. 'harmonized').
    :param dtype: The datatype of the files to download (e.g. 'clinical')
    :param data_file: The file name to download (look at index.tsv for examples).

    :return: True if data has successfully downloaded; raises error otherwise.
    """
    # Input validation
    if not cancer or not source or not dtype:
        raise InvalidParameterError("Cancer, source, and datatypes must be provided.")
    
    # Prepare for data download
    # The procedure varies depending on the source and type of data
    if source in ['harmonized', 'mssm']:
        description = f"{source}-all_cancers-{dtype}"
    elif source in ['washu'] and dtype in ['tumor_purity', 'hla_typing']:
        description = f"{source}-all_cancers-{dtype}"
    else:
        description = f"{source}-{cancer}-{dtype}"
    output_dir = '-'.join(description.split('-')[:-1])
    os.makedirs(os.path.join(DATA_DIR, output_dir), exist_ok=True)
    # Download requested dataframe
    file_name = f"{description}-{data_file}"
    output_file = os.path.join(output_dir, data_file)
    # Error handling for different exceptions
    try:
        get_data(f"{get_bucket()}/{file_name}", output_file)
        # Verify checksum
        with open(os.path.join(DATA_DIR, output_file), 'rb') as data_file:
            local_hash = md5(data_file.read()).hexdigest()
        if "md5:"+local_hash != cptac.INDEX.query("filename == @file_name")['checksum'].item():
            os.remove(os.path.join(DATA_DIR, output_file))
            raise DownloadFailedError("Download failed: local an remote files do not match. Please try again.")
        return True
    except NoInternetError as e:
        raise NoInternetError("Download failed -- No internet connection.")
    except DownloadFailedError as e:
        raise e
    except requests.RequestException as e:
        raise HttpResponseError(f"Requesting data failed with the following error: {e}")
    except Exception as e:
        raise DownloadFailedError(f"Failed to download data file for {source} {cancer} {dtype} with error:\n{e}") from e


def get_bucket() -> str:
    """
    Gets the bucket in Zenodo that houses all data files.
    :return: The URL of the Zenodo bucket containing the data files.
    """
    if BUCKET is None:
        init_files()
    return BUCKET


def download_chunk(url, start, end, file_path, pbar=None):
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
    block_size = 1024 # We can adjust this value as needed

    with open(file_path, 'rb+') as data_file:
        data_file.seek(start)
        for data in response.iter_content(block_size):
            data_file.write(data)
            if pbar is not None:
                pbar.update(len(data)) # Update the progress bar based on the size of the data block


def get_data(url: str, subfolder: str = '', num_threads: int = 4) -> str:
    """
    Downloads a file using multithreading and saves it to the specified subfolder.

    :param url: The URL of the file to download.
    :param subfolder: The subfolder where the downloaded file should be saved.
    :param num_threads: The number of threads to use for downloading the file (default is 4).
    :return: The path of the downloaded file.
    """
    response = requests.head(url, headers=AUTH_HEADER)
    response.raise_for_status()

    file_size = int(response.headers['content-length'])
    chunk_size = file_size // num_threads

    # Create an empty file with the same size as the file to be downloaded
    with open(os.path.join(DATA_DIR, subfolder), 'wb') as data_file:
        data_file.truncate(file_size)

    with ThreadPoolExecutor(max_workers=num_threads) as executor, tqdm(total=file_size, unit='B', unit_scale=True, desc=f"Downloading {os.path.split(subfolder)[1]}") as pbar:
        futures = []
        for i in range(num_threads):
            start = i * chunk_size
            end = start + chunk_size - 1 if i != num_threads - 1 else file_size - 1
            futures.append(executor.submit(download_chunk, url, start, end, os.path.join(DATA_DIR, subfolder), pbar))

        for future in concurrent.futures.as_completed(futures):
            future.result()  # Raise any exception encountered during download

    return subfolder

