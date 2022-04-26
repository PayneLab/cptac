#   Copyright 2018 Samuel Payne sam_payne@byu.edu
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import os

import cptac
from cptac.tools.download_tools.pdc_download import pdc_download
from cptac.tools.download_tools.box_download import box_download

from cptac.exceptions import CptacDevError, DataSourceNotFoundError, InvalidParameterError, NoInternetError




def download(sources, cancers='all', version="latest", redownload=False):
    """Download data files for the specified cancers, sources, and datatypes. Defaults to downloading latest version on server.

    Parameters:
    sources (dict): Keys are source names (str), values are the datatypes (list of str)
    cancers (list of str): The cancers for which the sources/datatypes will be downloaded
    version (str, optional): Which version of the data files to download. Defaults to latest on server.
    redownload (bool, optional): Whether to redownload the data files, even if that version of the data is already downloaded. Default False.

    Returns:
    bool: Indicates whether download was successful.
    """

    # check if sources parameter is valid
    sources = _validate_sources(sources)
    
    # check if cancers parameter is valid
    cancers = _validate_cancers(cancers)

    # iterate through cancers and sources and download corresonding data files
    success = True
    for cancer in cancers:
        for source, datatypes in sources.items():
            if source == "pdc" and not pdc_download(cancer=cancer, version=version, redownload=redownload):
                success = False
            elif not box_download(cancer, source, datatypes, version=version, redownload=redownload):
                success = False

    return success

def _validate_sources(sources):
    all_sources = cptac.get_source_options()
    if type(sources) is not dict:
        raise InvalidParameterError("Sources must be a dict of form {'source':['datatypes']}. 'all' is a valid datatype.")
    else:
        for s, datatypes in sources.items():
            if s not in all_sources:
                raise DataSourceNotFoundError(f"{s} is not a valid source! Call cptac.list_datasets() for valid options.")
            if type(datatypes) is str:
                sources[s] = list([datatypes])
        return sources

def _validate_cancers(cancers):
    all_cancers = cptac.get_cancer_options()
    if type(cancers) is str and cancers == 'all':
        return all_cancers
    elif type(cancers) is list and len(cancers) == 1 and cancers[0] == 'all':
        return all_cancers
    elif type(cancers) is str and cancers in all_cancers:
        return list([cancers])
    elif type(cancers) is list:
        invalid_cancers = list()
        for c in cancers:
            if c not in all_cancers:
                invalid_cancers.append(c)
        if len(invalid_cancers) > 0:
            raise InvalidParameterError(f"{invalid_cancers} are not a valid cancers. Run cptac.list_datasets() to see valid cancer types.")
    else: # handle case where cancers is an invalid string
        raise InvalidParameterError(f"{cancers} is not a valid cancer. Run cptac.Run cptac.list_datasets() to see valid cancer types.")

def update_index(dataset):
    """Check if the index of the given dataset is up to date with server version, and update it if needed.

    Parameters:
    dataset (str): The name of the dataset to check the index of.

    Returns:
    bool: Indicates if we were able to check the index and update if needed (i.e. we had internet)
    """
    # Get the path to our dataset
    dataset_path = get_dataset_path(dataset)

    # Define our file names we'll need
    index_urls_file = "index_urls.tsv"
    index_hash_file = "index_hash.txt"
    index_file = "index.txt"

    # Get, from the server, what the md5 hash of our index file should be
    index_urls_path = os.path.join(dataset_path, index_urls_file)
    urls_dict = parse_tsv_dict(index_urls_path)
    index_hash_url = urls_dict.get(index_hash_file)

    checking_msg = f"Checking that {dataset} index is up-to-date..."
    print(checking_msg, end='\r')
    try:
        server_index_hash = download_text(index_hash_url)
    finally:
        print(" " * len(checking_msg), end='\r') # Erase the checking message, even if there was an internet error

    index_path = os.path.join(dataset_path, index_file)

    if os.path.isfile(index_path):
        local_index_hash = hash_file(index_path)
        if local_index_hash == server_index_hash:
            return True

    index_url = urls_dict.get(index_file)
    download_file(url=index_url, path=index_path, server_hash=server_index_hash, file_message=f"{dataset} index")

    if os.path.isfile(index_path):
        local_index_hash = hash_file(index_path)
        if local_index_hash == server_index_hash:
            return True

    # If we get here, something apparently went wrong with the download.
    raise NoInternetError("Insufficient internet. Check your internet connection.")

