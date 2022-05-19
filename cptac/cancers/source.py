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

from os import path

import cptac

from cptac import CPTAC_BASE_DIR
from cptac.exceptions import DataTypeNotInSourceError, InvalidDataVersionError, MissingFileError
from cptac.tools.dataframe_tools import *

class Source:

    def __init__(self, cancer_type, source, version, valid_versions, data_files, load_functions, no_internet):
        self.no_internet = no_internet
        self.source = source
        self.cancer_type = cancer_type
        self._data = list()
        self.data_files = data_files
        self.load_functions = load_functions
        self.version = self.set_version(version)
        self.valid_versions = valid_versions
        
    def get_df(self, df_type):

        # if that df hasn't been loaded yet, load it
        if df_type not in self._data:

            # check to see if the df type requested is availabe from this source
            if df_type not in self.load_functions:
                raise DataTypeNotInSourceError(f"The {self.source} source does not have {df_type} data for {self.cancer_type} cancer.")
            else:
                # print loading message
                loading_msg = f"Loading {df_type} dataframe for {self.source} {self.cancer_type} (v{self.version()})"
                loading_msg = loading_msg + "."
                print(loading_msg, end='\r')
                
                # call the load function for that df type
                self.load_functions[df_type]()

                # Call function from dataframe_tools.py to standardize the names of the index and column axes
                if self.source in ['awg', 'awgconf']:
                    standardize_axes_names(self._data[df_type])

                # Erase the loading message
                print(' ' * len(loading_msg), end='\r') 

        return self._data[df_type]

    def set_version(self, version):

        # check if version is valid
        if version not in self.valid_versions:
            raise InvalidDataVersionError(f"{version} is not a valid version. These are the valid versions: {self.valid_versions}")
        else:
            self.version = version

    def version(self):
        return self.version

    
def get_file_path(self, df_type, data_file):
    """For dataset loading. Check that a version is installed, then return the paths to the data files for that version.

    Parameters:
    data_files (list of strings): The file names to get paths for.

    Returns:
    string: The path to the given data file for the currently set version of the source.
    """
    # Get our dataset path and index
    file_path = path.join(CPTAC_BASE_DIR, f"data/data_{self.source}_{self.cancer_type}/v_{self.version}/{data_file}")

    if path.isdir(file_path):
        return file_path

    elif not path.isdir(file_path) and not self.no_internet:
        return "not downloaded"

    # Raise error if file is not installed and they don't have an internet connection
    if not path.isdir(file_path) and self.no_internet:
        raise MissingFileError(f"The {self.source} {df_type} file for the {self.cancer_type} is not downloaded and you are running cptac in no_internet mode.")

def perform_initial_checks(self, df_type):
    # check if df_type is valid for the set version
    if self.version not in self.data_files:
        raise InvalidDataVersionError(f"{df_type} is not available in the data freeze for version v_{self.version} of {self.source} {self.cancer_type} data.")
    
    # get the file path
    f = self.data_files[self.version][df_type]
    file_path = self.get_file_path(f)
    
    # if the file hasn't been downloaded and they have internet, download it
    if file_path == "not downloaded":
        cptac.download(sources={self.source : [df_type]}, cancers=self.cancer_type, version=self.version)
        file_path = self.get_file_path(f)

    return file_path