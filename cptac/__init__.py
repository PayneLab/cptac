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

import webbrowser
import textwrap
import os

def list_data():
    """List all available datasets."""
    print("Available datasets:")
    print("endometrial")
    print("ovarian")
    print("colon")

def list_api():
    """Print docstrings for all accessible functions."""
    help(__name__)

def embargo():
    """Open CPTAC embargo details in web browser."""
    print("Opening embargo details in web browser...")
    webbrowser.open("https://proteomics.cancer.gov/data-portal/about/data-use-agreement")

def version():
    """Print version number of cptac package."""
    version = {}
    with open(dir_path + os.sep + "version.py") as fp:
    	exec(fp.read(), version)
    return(version['__version__'])

dir_path = os.path.dirname(os.path.realpath(__file__))
message = "Welcome to the cptac data service package. Available datasets may be viewed using cptac.list_data(). In order to access a specific data set, import a cptac subfolder using either \'import cptac.dataset\' or \'from cptac import dataset\'.\n"
wrapped_list = textwrap.wrap(message)
for line in wrapped_list:
    print(line)
print("******")
print("Version:",version())
print("******")
