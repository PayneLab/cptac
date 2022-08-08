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
import webbrowser
import pandas as pd

# cptac base path
CPTAC_BASE_DIR = path.abspath(path.dirname(__file__))

# Function imports
from cptac.tools.download_tools.download import download
from cptac.tools.download_tools.box_download import download_text as _download_text
from cptac.exceptions import CptacError, CptacWarning, InvalidParameterError, NoInternetError, OldPackageVersionWarning

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

# auth import
from cptac.tools.auth_tools.box_auth import BoxAuth
box_auth = BoxAuth()

#### This code generates the __OPTIONS__ dataframe which shows all possible cancer, source, datatype combinations
def _load_options():
    """Load the tsv file with all the possible cancer, source, datatype combinations"""
    options_file = path.join(CPTAC_BASE_DIR, "options.tsv")
    df = pd.read_csv(options_file, sep="\t")
    return df

__OPTIONS__ = _load_options()

def get_options():
    return __OPTIONS__.copy()

def get_cancer_options():
    df = __OPTIONS__.copy()
    return df["Cancers"].unique()

def get_source_options():
    df = __OPTIONS__.copy()
    return df["Sources"].unique()

def list_datasets():
    """List all available datasets."""
    df = __OPTIONS__.\
    copy().\
    drop("Loadable datatypes", axis=1)

    df = df.\
    assign(Datatypes=df["Datatypes"].str.split("\ *,\ *", expand=False, regex=True)).\
    explode("Datatypes").\
    reset_index(drop=True)

    # Print our dataframe as a pretty tree structure
    info = {}
    for row in df.set_index(["Cancers", "Sources", "Datatypes"]).index.values:
        if row[0] not in info.keys():
            info[row[0]] = {}
        if row[1] not in info[row[0]].keys():
            info[row[0]][row[1]] = []
        info[row[0]][row[1]].append(row[2])

    df_tree = _tree(info)
    print(df_tree)

    #if the dataframe is needed it can be returned. If not,
    #the python interpreter will print anything that is returned, so no return for now
    #return df

def _tree(nest, prepend=""):
    """Recursively build a formatted string to represent a dictionary"""
    tree_str = ""
    if isinstance(nest, dict):
        for i, (k, v) in enumerate(nest.items()):
            if i == len(nest.keys()) - 1:
                branch = "└"
                newprepend = prepend + "    "
            else:
                branch = "├"
                newprepend = prepend + "│   "
            tree_str += f"{prepend}{branch}── {k}\n"
            tree_str += _tree(nest=v, prepend=newprepend)
    elif isinstance(nest, list):
        for i, v in enumerate(nest):
            if i == len(nest) - 1:
                branch = "└"
            else:
                branch = "├"
            tree_str += f"{prepend}{branch}── {v}\n"
    else:
        raise ValueError(f"Unexpected type '{type(nest)}'")

    return tree_str

#### End __OPTIONS__ code

def embargo():
    """Open CPTAC embargo details in web browser."""
    message = "Opening embargo details in web browser..."
    print(message, end = '\r')
    webbrowser.open("https://proteomics.cancer.gov/data-portal/about/data-use-agreement")
    print(" " * len(message), end='\r') # Erase the message

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

#### Helper functions for handling exceptions and warnings

# Because Python binds default arguments when the function is defined,
# default_hook's default will always refer to the original sys.excepthook
def _exception_handler(exception_type, exception, traceback, default_hook=sys.excepthook): 
    """Catch cptac-generated exceptions, and make them prettier."""
    if issubclass(type(exception), CptacError):
        print(f"cptac error: {str(exception)} ({traceback.tb_frame.f_code.co_filename}, line {traceback.tb_lineno})", file=sys.stderr) # We still send to stderr
    else:
        default_hook(exception_type, exception, traceback) # This way, exceptions from other packages will still be treated the same way

def _warning_displayer(message, category, filename, lineno, file=None, line=None, default_displayer=warnings.showwarning): # Python binds default arguments when the function is defined, so default_displayer's default will always refer to the original warnings.showwarning
    """Catch cptac-generated warnings and make them prettier."""
    if issubclass(category, CptacWarning):
        print(f"cptac warning: {str(message)} ({filename}, line {lineno})", file=sys.stderr) # We still send to stderr
    else:
        default_displayer(message, category, filename, lineno, file, line) # This way, warnings from other packages will still be displayed the same way

sys.excepthook = _exception_handler # Set our custom exception hook
warnings.showwarning = _warning_displayer # And our custom warning displayer
warnings.simplefilter("always", category=CptacWarning) # Edit the warnings filter to show multiple occurences of cptac-generated warnings

def check_version():
    """Check in background whether the package is up-to-date"""
    version_url = "https://byu.box.com/shared/static/kbwivmqnrdnn5im2gu6khoybk5a3rfl0.txt"
    try:
        remote_version = _download_text(version_url)
    except NoInternetError:
        pass
    else:
        local_version = version()
        if remote_version != local_version:
            warnings.warn(f"Your version of cptac ({local_version}) is out-of-date. Latest is {remote_version}. Please run 'pip install --upgrade cptac' to update it.", OldPackageVersionWarning, stacklevel=2)

version_check_thread = threading.Thread(target=check_version)
version_check_thread.start() # We don't join because we want this to just finish in the background and not block the main thread
