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

from cptac.exceptions import DataSourceNotFoundError, InvalidParameterError, NoInternetError


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

    # variable for tracking download success
    success = True

    # handle special harmonized and mssm cases
    special_cases = set(['harmonized', 'mssm']).intersection(set(sources.keys()))
    for case in special_cases:
        source = case
        datatypes = sources[case]
        if not box_download(cancer='brca', source=source, datatypes=datatypes, version=version, redownload=redownload):
            success = False
        del sources[case]

    # iterate through cancers and sources and download corresonding data files
    for cancer in cancers:
        for source, datatypes in sources.items():
            if source == "pdc":
                # download the mapping files
                if not box_download(cancer=cancer, source=source, datatypes=None, version=version, redownload=redownload):
                    success = False
                if not pdc_download(cancer=cancer, datatypes=datatypes, version=version, redownload=redownload):
                    success = False
            
            elif not box_download(cancer, source, datatypes, version=version, redownload=redownload):
                success = False

    return success

def _validate_sources(sources):
    all_sources = cptac.get_source_options()
    if type(sources) is not dict:
        raise InvalidParameterError("Sources must be a dict of form {'source':['datatypes']}. ")
    else:
        for s, datatypes in sources.items():
            if s not in all_sources:
                raise DataSourceNotFoundError(f"{s} is not a valid source! Call cptac.list_datasets() for valid options.")
            if type(datatypes) is str:
                sources[s] = list([datatypes])
        return sources

def _validate_cancers(cancers):
    all_cancers = cptac.get_cancer_options()
    if cancers in ['all', ['all']]:
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
        else:
            return cancers
    else: # handle case where cancers is an invalid string
        raise InvalidParameterError(f"{cancers} is not a valid cancer. Run cptac.Run cptac.list_datasets() to see valid cancer types.")