import getpass
import glob
import os
import packaging.version
import requests
import warnings
import zenodopy
import wget

import cptac
from cptac.exceptions import *

STATIC_DOI = '10.5281/zenodo.7897498'
INDEX_FILE_NAME = 'all_index.txt'
DATA_DIR = os.path.join(cptac.CPTAC_BASE_DIR, "data/")

def zeno_download(cancer, source, datatype):
    """
    Downloads data files for a specific cancer, source, and datatype from Zenodo

    :param cancer: The cancer type (e.g. 'brca').
    :param source: The data source (e.g. 'harmonized').
    :param datatype: The datatype of the files to download (e.g. 'clinical')
    """
    if not cancer or not source or not datatype:
        raise ValueError("Cancer, source, and datatypes must be provided.")
    
    index_path = download_index_file_if_needed()

    file_urls = get_file_urls(cancer, source, datatype, index_path)

    if not file_urls:
        raise FileNotFoundError(f"No matching files found for cancer='{cancer}', source='{source}', datatype='{datatype}'")

    output_folder = os.path.join(DATA_DIR, f"data_{source}_{cancer}")

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for url in file_urls:
        download_file(url, output_folder)

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

def get_file_urls(cancer, source, datatype, index_path):
    """
    Gets the URLs of the files that match the given cancer, source, and datatypes

    :param cancer: The cancer type (e.g. 'brca)
    :param source: The data source (e.g. 'harmonized)
    :param datatype: The datatype of the files to download (e.g. 'clinical')
    :param index_path: The path to the index file.
    :return: A list of file URLs that match the given parameters
    """
    file_urls = []

    with open(index_path, 'r') as input:
        for line in input:
            columns = line.split('\t')
            identifiers = columns[0].split('_')

            if identifiers[0] == source and identifiers[1] == cancer and identifiers[2] == datatype:
                file_urls.append(columns[1].strip())

    return file_urls

def download_file(url, output_folder):
    """
    Downloads a file from the given URL to the specified output folder if it isn't already downloaded

    :param url: The URL of the file to download
    :param output_folder: The folder where the file should be downloaded
    """
    if not url or not output_folder:
        raise ValueError("URL and output folder must be provided.")
    
    file_name = url.split('/')[-1]
    output_path = os.path.join(output_folder, file_name)

    if not os.path.exists(output_path):
        wget.download(url, output_path)