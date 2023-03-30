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
from unittest.loader import VALID_MODULE_NAME
import pandas as pd
import requests
import warnings

import cptac
from cptac.file_download import get_box_token
from cptac.exceptions import CptacDevError, InvalidParameterError, NoInternetError, PdcDownloadError

from .brca import SOURCES as BRCA_SOURCES
from .ccrcc import SOURCES as CCRCC_SOURCES
from .coad import SOURCES as COAD_SOURCES
from .gbm import SOURCES as GBM_SOURCES
from .hnscc import SOURCES as HNSCC_SOURCES
from .lscc import SOURCES as LSCC_SOURCES
from .luad import SOURCES as LUAD_SOURCES
from .ov import SOURCES as OV_SOURCES
from .ucec import SOURCES as UCEC_SOURCES
from .pdac import SOURCES as PDAC_SOURCES

VALID_SOURCES = ['bcm', 'broad', 'mssm', 'umich', 'washu', 'harmonized']

def download(dataset, version="latest", redownload=False, box_token=None):

    if type(dataset) is str:

        dataset = dataset.lower()
        if dataset.startswith("pancan") or dataset == "all":
            download = _download_by_cancer
    
    elif type(dataset) is dict:
        download = _download_by_source

    else:
        raise InvalidParameterError(f"the dataset parameter must be a string of the cancer type (e.g. 'pancanbrca') or a dictionary {{'source' : datatype}} (e.g. {{'broad': 'CNV'}})")

    if box_token is None:
        box_token = get_box_token()

    return download(dataset, version, redownload, box_token)

# Helper functions
def _download_by_cancer(dataset, version, redownload, box_token):
    '''Download PDC data by cancer (brca, hnscc, etc.)
    '''
    if dataset.startswith("pancan") or dataset == "all":

        if dataset == "pancanbrca":
            sources = BRCA_SOURCES
        elif dataset == "pancanccrcc":
            sources = CCRCC_SOURCES
        elif dataset == "pancancoad":
            sources = COAD_SOURCES
        elif dataset == "pancangbm":
            sources = GBM_SOURCES
        elif dataset == "pancanhnscc":
            sources = HNSCC_SOURCES
        elif dataset == "pancanlscc":
            sources = LSCC_SOURCES
        elif dataset == "pancanluad":
            sources = LUAD_SOURCES
        elif dataset == "pancanov":
            sources = OV_SOURCES
        elif dataset == "pancanucec":
            sources = UCEC_SOURCES
        elif dataset == "pancanpdac":
            sources = PDAC_SOURCES
        elif dataset == "all":
            sources = sorted(set(BRCA_SOURCES + CCRCC_SOURCES + COAD_SOURCES + GBM_SOURCES + HNSCC_SOURCES + LSCC_SOURCES + LUAD_SOURCES + OV_SOURCES + UCEC_SOURCES + PDAC_SOURCES))
        else:
            raise InvalidParameterError(f"{dataset} is not a valid dataset.")

        overall_success = True
        for source in sources:

            single_success = cptac.download(source, version=version, redownload=redownload, _box_auth=True, _box_token=box_token)

            if not single_success:
                overall_success = False

        return overall_success

    else:
        return cptac.download(dataset, version=version, redownload=redownload, _box_auth=True, _box_token=box_token)

def _download_by_source(dataset, version, redownload, box_token):
    '''Download PDC data by source (broad, washu, etc.)
    dataset: a dict of form {datatype : source} or {source : datatype}
        datatype should be a string
        source can be a string or a list of strings
    '''
    # verify all sources are valid
    for source in dataset.keys():
        invalid_sources = set()
        if source.lower() not in VALID_SOURCES:
            invalid_sources.add(source)
            
        if len(invalid_sources) > 0:
            raise InvalidParameterError(f"Invalid source(s) detected: {invalid_sources}. Valid sources include {VALID_SOURCES}.")

    # all valid sources for all cancer types
    SOURCES = sorted(set(BRCA_SOURCES + CCRCC_SOURCES + COAD_SOURCES + GBM_SOURCES + HNSCC_SOURCES + LSCC_SOURCES + LUAD_SOURCES + OV_SOURCES + UCEC_SOURCES + PDAC_SOURCES))
    
    overall_success = True
    for source, datatype in dataset.items():

        if type(dataset) is not list:
            datatypes = list([datatype])
        
        # get data for all cancer types from the given source
        for s in SOURCES:
            if s.startswith(source):
                single_success = cptac.download(dataset='all', source=source, datatypes=datatypes, version=version, redownload=redownload, _box_auth=True, _box_token=box_token)
                if not single_success:
                    overall_success = False

    return overall_success
