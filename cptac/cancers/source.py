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
import pandas as pd
from hashlib import md5
from warnings import warn

from cptac import CPTAC_BASE_DIR
from cptac.exceptions import DataTypeNotInSourceError, MissingFileError, FailedChecksumWarning
from cptac.tools.dataframe_tools import standardize_axes_names

class Source:
    """
    The Source class is a base class that provides methods to manage and interact with different data sources.
    
    Attributes:
        no_internet (bool): If set to True, the data source objects are initialized
                            without downloading data. Default is False.
        source (str): The name of the data source.
        cancer_type (str): The type of cancer that the data pertains to.
        _data (dict): Dictionary to hold the loaded dataframes.
        _helper_tables (dict): Dictionary to hold helper tables that support
                                certain functions.
        data_files (dict): Dictionary that holds the data file names for each
                            data type.
        load_functions (dict): Dictionary that holds the load function for each
                                data type.
    """

    def __init__(self, cancer_type, source, data_files, load_functions, no_internet):
        """
        The constructor for the Source class.

        Parameters:
            cancer_type (str): The type of cancer that the data pertains to.
            source (str): The name of the data source.
            data_files (dict): Dictionary that holds cthe data file names for each
                                data type.
            load_functions (dict): Dictionary that holds the load function for each
                                    data type.
            no_internet (bool): If set to True, the data source objects are
                                initialized without downloading data. Default is False.
        """
        self.no_internet = no_internet
        self.source = source
        self.cancer_type = cancer_type
        self._data = {}
        self._helper_tables = {}
        self.data_files = data_files
        self.load_functions = load_functions

    def get_df(self, df_type):
        """Get the dataframe of the specified data type

        Parameters:
        df_type (str): Name of datatype to return e.g. "proteomics".

        Returns:
        pandas.DataFrame: The dataframe of the desired datatype.
        """
        # if that df hasn't been loaded yet, load it
        if df_type not in self._data:
            # check to see if the df type requested is availabe from this source
            if df_type not in self.load_functions:
                raise DataTypeNotInSourceError(f"The {self.source} source does not have {df_type} data for {self.cancer_type} cancer.")
            self.load_functions[df_type]()
        return self._data[df_type]

    def save_df(self, df_type, df):
        """Perform final formatting so all cptac data is consistent, and save the dataframe in the self._data dictionary with df_type as the key"""

        # set the index name to "Patient_ID" and the columns name to "Name"
        standardize_axes_names(df)

        # Sort the dataframe based off sample status (tumor or normal), then alphabetically
        df = df.sort_index()
        #'.N' for normal, '.C' for cored normals (in HNSCC)
        normal = df.loc[df.index.str.contains(r'\.[NC]$', regex = True, na = False)]
        # Tumor samples don't have any special endings cohorts for now
        tumor = df.loc[~ df.index.str.contains(r'\.[NC]$', regex = True, na = False)]
        df = pd.concat([tumor, normal])

        self._data[df_type] = df


    def locate_files(self, datatype):
        """Checks if the datatype is valid
        If the datatype is valid, it finds the file path, downloading the files if necessary

        Parameters:
            datatype (str): The datatype to get all filepaths for.

        Returns:
            A single file path or a list of file paths to all files of the given datatype
        """
        # pull the file name or list of file names from the self.data_files dict
        data_files = self.data_files[datatype]
        # self.data_files can either contain a single string or a list of strings. Let's work with all as a list for now.
        if type(data_files) != list:
            data_files = [data_files]
        file_paths = []
        # Locate and download each data_file
        for data_file in data_files:
            # dataset = self.source if self.source in ['harmonized', 'mssm'] else f"{self.source}_{self.cancer_type}"
            # This should eventually be handled within the respective sources, but this will do for now
            cancer_type = "all_cancers" if self.source in ['mssm', 'harmonized'] or self.source in['washu'] and datatype in ['tumor_purity', 'hla_typing'] else self.cancer_type
            dataset = f"{self.source}-{cancer_type}"
            file_path = os.path.join(CPTAC_BASE_DIR, f"data/{dataset}/{data_file}")
            prefixed_file = f"{self.source}-{cancer_type}-{datatype}-{data_file}"
            # Ensure data is not corrupted, download files if needed
            if os.path.isfile(file_path) and not self.no_internet: # It's pointless to check the checksum if we can't redownload it
                with open(file_path, 'rb') as in_file:
                    local_hash = f"md5:{md5(in_file.read()).hexdigest()}"
                if local_hash != cptac.INDEX.loc[cptac.INDEX['filename']==prefixed_file, 'checksum'].item():
                    warn(FailedChecksumWarning("Local file and online file have different checksums; redownloading data"))
                    os.remove(file_path)
                
            if not os.path.isfile(file_path) and not self.no_internet:
                cptac.download(self.cancer_type, self.source, datatype, data_file)
            elif not os.path.isfile(file_path) and self.no_internet:
                raise MissingFileError(f"The {self.source} {data_file} file for the {self.cancer_type} is not downloaded and you are running cptac in no_internet mode.")

            file_paths.append(file_path)
        
        return file_paths if len(file_paths) >= 2 else file_paths[0]
