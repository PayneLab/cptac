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
import io
import sys
import warnings

# Function imports
from .file_download import download
from .file_download import download_text as _download_text
from .exceptions import CptacError, CptacWarning, InvalidParameterError, NoInternetError, OldPackageVersionWarning

# Dataset imports
from .brca import Brca
from .ccrcc import Ccrcc
from .colon import Colon
from .endometrial import Endometrial
from .gbm import Gbm
from .hnscc import Hnscc
from .lscc import Lscc
from .luad import Luad
from .ovarian import Ovarian

def list_datasets():
    """List all available datasets."""

    dataset_list_url = "https://byu.box.com/shared/static/5vwsvgu8fyzao0pb7rx8lt26huofcdax.tsv"

    try:
        dataset_list_text = _download_text(dataset_list_url)
    except NoInternetError:
        raise NoInternetError("Insufficient internet to download available dataset info. Check your internet connection.") from None

    return pd.read_csv(io.StringIO(dataset_list_text), sep="\t", index_col=0)

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

# Helper functions for handling exceptions and warnings
def _exception_handler(exception_type, exception, traceback, default_hook=sys.excepthook): # Because Python binds default arguments when the function is defined, default_hook's default will always refer to the original sys.excepthook
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

# Check whether the package is up-to-date
version_url = "https://byu.box.com/shared/static/kbwivmqnrdnn5im2gu6khoybk5a3rfl0.txt"
try:
    remote_version = _download_text(version_url)
except NoInternetError:
    pass
else:
    local_version = version()
    if remote_version != local_version:
        warnings.warn(f"Your version of cptac ({local_version}) is out-of-date. Latest is {remote_version}. Please run 'pip install --upgrade cptac' to update it.", OldPackageVersionWarning, stacklevel=2)
