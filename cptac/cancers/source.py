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
        self._data = {}
        self._helper_tables = {}
        self.data_files = data_files
        self.load_functions = load_functions
        self.set_version(version)
        self.valid_versions = valid_versions

    def get_df(self, df_type):
        """Get the dataframe of the specified data type
        I'm thinking about having tissue_type do its thing here. If all dfs are parsed to have normal samples clearly identified,
        then this function can return the whole data frame for "both", or a subset of it for "tumor" or "normal"

        Parameters:
        df_type (str): Name of datatype to return e.g. "proteomics".
        tissue_type (str): Acceptable values in ["tumor","normal","both"]. Specifies the desired tissue type desired in the dataframe. Defaults to "both".

        Returns:
        pandas.DataFrame: The dataframe of the desired datatype.
        """

        # if that df hasn't been loaded yet, load it
        if df_type not in self._data:

            # check to see if the df type requested is availabe from this source
            if df_type not in self.load_functions:
                raise DataTypeNotInSourceError(f"The {self.source} source does not have {df_type} data for {self.cancer_type} cancer.")
            else:
                # print loading message
                loading_msg = f"Loading {df_type} dataframe for {self.source} {self.cancer_type} (v{self.get_version()})"
                #print(loading_msg, end='\r')

                # call the load function for that df type
                self.load_functions[df_type]()

                # Erase the loading message
                print(' ' * len(loading_msg), end='\r')

        return self._data[df_type]


    def save_df(self, df_type, df):
        """Perform final formatting so all cptac data is consistent, and save the dataframe in the self._data dictionary with df_type as the key"""

        # set the index name to "Patient_ID" and the columns name to "Name"
        standardize_axes_names(df)

        # Sort the dataframe based off sample status (tumor or normal), then alphabetically
        df = df.sort_index()
        #'.N' for normal, '.C' for cored normals (in HNSCC)
        normal = df.loc[df.index.str.contains('\.[NC]$', regex = True, na = False)]
        # Tumor samples don't have any special endings except in the awg confirmatory cohorts for now
        tumor = df.loc[~ df.index.str.contains('\.[NC]$', regex = True, na = False)]
        df = pd.concat([tumor, normal])

        self._data[df_type] = df


    def set_version(self, version):
        # check if version is valid
        if version == "latest":
            version = sorted(self.valid_versions)[-1]
        if version not in self.valid_versions:
            raise InvalidDataVersionError(f"{version} is not a valid version. These are the valid versions: {self.valid_versions}")
        else:
            self.version = version


    def get_version(self):
        return self.version


    def get_file_path(self, data_file):
        """Return the path to a specific data file

        Parameters:
        data_file (str): The file name to get a filepath for

        Returns:
        string: The path to the given data file for the currently set version of the source.
        """
        # Get our dataset path and index
        if self.source in ["harmonized", "mssm"]:
            dataset = self.source
        else:
            dataset = self.source + "_" + self.cancer_type

        file_path = path.join(CPTAC_BASE_DIR, f"data/data_{dataset}/{dataset}_v{self.version}/{data_file}")

        if path.isfile(file_path):
            return file_path

        elif not path.isdir(file_path) and not self.no_internet:
            return "not downloaded"

        # Raise error if file is not installed and they don't have an internet connection
        if not path.isdir(file_path) and self.no_internet:
            raise MissingFileError(f"The {self.source} {data_file} file for the {self.cancer_type} is not downloaded and you are running cptac in no_internet mode.")


    def locate_files(self, datatype):
        """Checks if the datatype is valid for the source's set version. 
        If the datatype is valid, it finds the file path, downloading the files if necessary

        Parameters:
            datatype (str): The datatype to get all filepaths for.

        Returns:
            A single file path or a list of file paths to all files of the given datatype
        """
        # check if datatype is valid for the set version
        if self.version not in self.data_files:
            raise InvalidDataVersionError(f"{datatype} is not available in the data freeze for version v_{self.version} of {self.source} {self.cancer_type} data.")

        # pull the file name or list of file names from the self.data_files dict
        f = self.data_files[self.version][datatype]

        # get the file path(s)
        if type(f) == list:
            file_paths = list()
            for file_name in f:
                path = self.get_file_path(file_name)
                 # if the file hasn't been downloaded and they have internet, download it
                if path == "not downloaded":
                    cptac.download(sources={self.source : [datatype]}, cancers=self.cancer_type, version=self.version)
                    path = self.get_file_path(file_name)
                # append path to file paths
                file_paths.append(path)

            return file_paths

        else:
            file_path = self.get_file_path(f)
            if file_path == "not downloaded":
                cptac.download(sources={self.source : [datatype]}, cancers=self.cancer_type, version=self.version)
                file_path = self.get_file_path(f)

            return file_path
