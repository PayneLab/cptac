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

import pandas as pd
import webbrowser
import os.path as path
import sys
import warnings
from .file_download import download
from .exceptions import CptacError, CptacWarning
from .brca import Brca
from .endometrial import Endometrial
from .colon import Colon
from .luad import Luad
from .ovarian import Ovarian
from .renalccrcc import RenalCcrcc

def list_datasets():
    """List all available datasets."""
    col_names = ["Description", "Data reuse status", "Publication link"]
    col_index = pd.Index(data=col_names, name="Dataset name")
    datasets = {
        "Brca": ["breast cancer", "password access only", "unpublished"],
        "Colon": ["colorectal cancer", "no restrictions", "https://www.ncbi.nlm.nih.gov/pubmed/31031003"],
        "Endometrial": [ "endometrial carcinoma (uterine)", "no restrictions", "unpublished"],
        "Luad": ["lung adenocarcinoma", "password access only", "unpublished"],
        "Ovarian": ["high grade serous ovarian cancer", "no restrictions", "unpublished"],
        "RenalCcrcc": ["clear cell renal cell carcinoma (kidney)", "no restrictions", "unpublished"],
        }
    dataset_df = pd.DataFrame(data=datasets, index=col_index)
    dataset_df = dataset_df.transpose()
    dataset_df.index.name = "" # Giving the index a name, even though it's an emtpy string, causes a space to be printed between the column names and the first row, which improves readability.
    print(f"Available datasets:\n\n{dataset_df}")

def embargo():
    """Open CPTAC embargo details in web browser."""
    message = "Opening embargo details in web browser..."
    print(message, end = '\r')
    webbrowser.open("https://proteomics.cancer.gov/data-portal/about/data-use-agreement")
    print(" " * len(message), end='\r') # Erase the message

def version():
    """Return version number of cptac package."""
    version = {}
    path_here = path.abspath(path.dirname(__file__))
    version_path = path.join(path_here, "version.py")
    with open(version_path) as fp:
        exec(fp.read(), version)
    return(version['__version__'])

def how_to_cite():
    """Give instructions for citing CPTAC datasets."""
    print("For instructions on how to cite a specific dataset, please call its how_to_cite method, e.g. cptac.Endometrial().how_to_cite()")

def _exception_handler(exception_type, exception, traceback, default_hook=sys.excepthook): # Because Python binds default arguments when the function is defined, default_hook's default will always refer to the original sys.excepthook
    """We're going to catch cptac-generated exceptions, and make them prettier."""
    if issubclass(type(exception), CptacError):
        print(str(exception))
    else:
        default_hook(exception_type, exception, traceback) # This way, exceptions from other packages will still be treated the same way

def _warning_displayer(message, category, filename, lineno, file=None, line=None, default_displayer=warnings.showwarning): # Because Python binds default arguments when the function is defined, default_displayer's default will always refer to the original warnings.showwarning
    """We're also going to catch cptac-generated warnings and make them prettier."""
    if issubclass(category, CptacWarning):
        print("Warning: " + str(message))
    else:
        default_displayer(message, category, filename, lineno, file, line) # This way, warnings from other packages will still be displayed the same way

sys.excepthook = _exception_handler # Set our custom exception hook
warnings.showwarning = _warning_displayer # And our custom warning displayer