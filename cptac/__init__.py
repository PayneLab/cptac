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
import os.path as path
from .endometrial import Endometrial
from .colon import Colon
from .ovarian import Ovarian

def list_data():
    """List all available datasets."""
    print("Available datasets:")
    datasets = [
        "Colon",
        "Endometrial",
        "Ovarian",]
    for dataset in sorted(datasets):
        print("\t" + dataset)

def embargo():
    """Open CPTAC embargo details in web browser."""
    print("Opening embargo details in web browser...")
    webbrowser.open("https://proteomics.cancer.gov/data-portal/about/data-use-agreement")

def version():
    """Return version number of cptac package."""
    version = {}
    path_here = path.abspath(path.dirname(__file__))
    version_path = path.join(path_here, "version.py")
    with open(version_path) as fp:
        exec(fp.read(), version)
    return(version['__version__'])

message = "Welcome to cptac, a python package for disseminating cancer proteogenomics data. To view available datasets, enter 'cptac.list_data()'. Extensive tutorials are available at https://github.com/PayneLab/cptac/tree/master/doc\n"
wrapped_list = textwrap.wrap(message)
for line in wrapped_list:
    print(line)

print("\n******\nVersion: {}\n******".format(version()))
