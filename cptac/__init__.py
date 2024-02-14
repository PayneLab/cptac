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

import os.path as path
import sys
import threading
import warnings
import pandas as pd
import requests
import cptac

# cptac base path
CPTAC_BASE_DIR = path.abspath(path.dirname(__file__))
# Some websites don't like requests from sources without a user agent. Let's preempt that issue.
USER_AGENT = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.10; rv:39.0)'
HEADERS = {'User-Agent': USER_AGENT}

# Function imports
from cptac.tools.download_tools import download, init_files, DATA_DIR
from cptac.exceptions import CptacError, CptacWarning, NoInternetError, OldPackageVersionWarning
from cptac.utils.other_utils import df_to_tree

# Dataset imports
from cptac.cancers.brca import Brca
from cptac.cancers.ccrcc import Ccrcc
from cptac.cancers.coad import Coad
from cptac.cancers.gbm import Gbm
from cptac.cancers.hnscc import Hnscc
from cptac.cancers.lscc import Lscc
from cptac.cancers.luad import Luad
from cptac.cancers.ov import Ov
from cptac.cancers.pdac import Pdac
from cptac.cancers.ucec import Ucec


#### Create custom exception and warning hooks to simplify error messages for new users
def _exception_handler(exception_type, exception, traceback, default_hook=sys.excepthook): 
    """Catch cptac-generated exceptions, and make them prettier."""
    if issubclass(type(exception), CptacError):
        print(f"cptac error: {str(exception)} ({traceback.tb_frame.f_code.co_filename}, line {traceback.tb_lineno})", file=sys.stderr)
    else:
        default_hook(exception_type, exception, traceback)

def _warning_displayer(message, category, filename, lineno, file=None, line=None, default_displayer=warnings.showwarning):
    """Catch cptac-generated warnings and make them prettier."""
    if issubclass(category, CptacWarning):
        print(f"cptac warning: {str(message)} ({filename}, line {lineno})", file=sys.stderr)
    else:
        default_displayer(message, category, filename, lineno, file, line)

sys.excepthook = _exception_handler
warnings.showwarning = _warning_displayer
warnings.simplefilter("always", category=CptacWarning)

#### Create the index for the data files (file lookup table)
try:
    init_files()
except NoInternetError:
    if path.exists(path.join(CPTAC_BASE_DIR, 'data', 'index.tsv')):
        pass
    else:
        raise NoInternetError("Unable to initialize cptac without index file. Please run the package at least once with an internet connection.")

INDEX = pd.read_csv(path.join(CPTAC_BASE_DIR, 'data', 'index.tsv'), sep='\t')
#### Generates the OPTIONS dataframe which shows all possible cancer, source, datatype combinations
def _load_options():
    """Load the tsv file with all the possible cancer, source, datatype combinations"""
    options_df = pd.DataFrame(INDEX['description'].str.split('-').tolist())
    options_df.columns = ['Source', 'Cancer', 'Datatype']
    options_df = options_df[['Cancer', 'Source', 'Datatype']]
    # options_df.loc[options_df.iloc[:2].str.contains('miRNA')] = 'miRNA' # condense all forms of micro RNA
    # options_df = options_df.unique().reset_index(drop=True)
    return options_df

OPTIONS = _load_options()

def list_datasets(*, condense_on = None, column_order = None, print_tree=False):
    """
    List all available datasets.
    
    :param condense_on (list): A list of column names. Values in selected columns will be aggregated into a list.
    :param print_tree (bool): If True, returns the database split in a pretty tree.
    """
    df = OPTIONS.copy()
    df = df[df['Datatype'] != 'mapping'].reset_index(drop=True)
    if column_order is None:
        column_order = df.columns
    if type(condense_on) == list:
        group_on_cols = [col for col in column_order if col not in condense_on]
        df = df.groupby(group_on_cols).agg({col: lambda x: list(set(x)) for col in condense_on})
    else:
        df = df[column_order]

    return df if not print_tree else df_to_tree(df)

def get_cancer_options():
    return list_datasets(condense_on=['Datatype'])

def get_cancer_info():
    cancer_abbreviations = {
        "brca": "Breast invasive carcinoma",
        "ccrcc": "Clear cell renal cell carcinoma",
        "coad": "Colon adenocarcinoma",
        "gbm": "Glioblastoma multiforme",
        "hnscc": "Head and Neck squamous cell carcinoma",
        "lscc": "Lung squamous cell carcinoma",
        "luad": "Lung adenocarcinoma",
        "ov": "Ovarian serous cystadenocarcinoma",
        "pda": "Pancreatic ductal adenocarcinoma",
        "pdac": "Pancreatic ductal adenocarcinoma",
        "ucec": "Uterine Corpus Endometrial Carcinoma",
        # Add more if needed
    }
    return cancer_abbreviations
    
def get_source_options():
    return list_datasets(condense_on=['Cancer'], column_order=['Source', 'Datatype', 'Cancer'])

def get_datatype_options():
    return list_datasets(condense_on=['Cancer'], column_order=['Datatype', 'Source', 'Cancer'])

#### End __OPTIONS__ code

# This website no longer works
# def embargo():
#     """Open CPTAC embargo details in web browser."""
#     message = "Opening embargo details in web browser..."
#     print(message, end = '\r')
#     webbrowser.open("https://proteomics.cancer.gov/data-portal/about/data-use-agreement")
#     print(" " * len(message), end='\r') # Erase the message

def version():
    """Return version number of cptac package."""
    version = {}
    version_path = path.join(CPTAC_BASE_DIR, "version.py")
    with open(version_path) as fp:
        exec(fp.read(), version)
    return version['__version__']

def how_to_cite():
    """Give instructions for citing CPTAC datasets."""
    print("If you use the API to generate results, please cite our manuscript describing the API - Lindgren et al. 2021, PMID:33560848, https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00919")
    print('\n')
    print("For instructions on how to cite a specific dataset, please call its how_to_cite method, e.g. cptac.Endometrial().how_to_cite()")


def _get_cancer_class(cancer_str):
    """
    Converts a string to a corresponding cancer class.

    Args:
    cancer_str (str): A string identifying the cancer type. This should match one of the following: 
        'brca', 'ccrcc', 'coad', 'gbm', 'hnscc', 'lscc', 'luad', 'ov', 'pdac', 'ucec', 'all_cancers'.
    
    Returns:
    A cancer class from the cptac module.

    Raises:
    ValueError: If the provided string does not correspond to a known cancer class.
    """
    # This dictionary should be updated as necessary
    mapping = {
        "brca": cptac.Brca,
        "ccrcc": cptac.Ccrcc,
        "coad": cptac.Coad,
        "gbm": cptac.Gbm,
        "hnscc": cptac.Hnscc,
        "lscc": cptac.Lscc,
        "luad": cptac.Luad,
        "ov": cptac.Ov,
        "pdac": cptac.Pdac,
        "ucec": cptac.Ucec,
        "all_cancers" : cptac.Ucec
    }

    try:
        return mapping[cancer_str.lower()]
    except KeyError:
        raise ValueError(f"'{cancer_str}' is not a known cancer class. Valid options are: {list(mapping.keys())}")

def download_cancer(cancer):
    """Downloads all datasets for a given cancer"""
    datasets = list_datasets()
    for _, row in datasets.iterrows():
        if row['Cancer'] == cancer:
            cancer_class = _get_cancer_class(row['Cancer'])
            cancer_instance = cancer_class()

            data = cancer_instance.get_dataframe(row['Datatype'], row['Source'])

    return 

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


def check_version():
    """Check in background whether the package is up-to-date"""
    version_url = "https://byu.box.com/shared/static/kbwivmqnrdnn5im2gu6khoybk5a3rfl0.txt"
    try:
        remote_version = download_text(version_url)
    except NoInternetError:
        pass
    else:
        local_version = version()
        if remote_version != local_version:
            warnings.warn(f"Your version of cptac ({local_version}) is out-of-date. Latest is {remote_version}. Please run 'pip install --upgrade cptac' to update it.", OldPackageVersionWarning, stacklevel=2)

version_check_thread = threading.Thread(target=check_version)
version_check_thread.start() # We don't join because we want this to just finish in the background and not block the main thread
