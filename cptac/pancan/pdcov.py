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
import numpy as np
import os
import warnings
import datetime

from cptac.dataset import Dataset
from cptac.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError

class PdcOv(Dataset):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        valid_versions = ["1.0"]

        data_files = {
            "1.0": [
                "clinical.tsv.gz",
                "phosphoproteome.tsv.gz",
                "proteome.tsv.gz",
                "OV_sample_TMT_annotation_UMich_GENCODE34_0315.csv"
            ]
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="pdcov", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet, attempt_update_index=False)

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

            if file_name == "clinical.tsv.gz":
                df = pd.read_csv(file_path, sep="\t", index_col=0)
                clin_drop_rows = ['JHU QC', 'PNNL-JHU Ref']
                df = df.drop(clin_drop_rows, axis = 'index') # drop quality control and reference
                self._data["clinical"] = df

            if file_name == "phosphoproteome.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                self._data["phosphoproteomics"] = df

            if file_name == "proteome.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                self._data["proteomics"] = df
                
            # This file maps Ov aliquots to patient IDs (case ID with tissue type) and 
            # can be found on Box under CPTAC/cptac/pancan/helper_files
            elif file_name == "OV_sample_TMT_annotation_UMich_GENCODE34_0315.csv":
                ov_map = pd.read_csv(file_path, sep = ",", usecols = ['specimen', 'sample'])
                ov_map = ov_map.loc[~ ov_map['sample'].str.contains('JHU', regex = True)] # drop quality control rows
                ov_map = ov_map.set_index('specimen')
                map_dict = ov_map.to_dict()['sample'] # create dictionary with aliquots as keys and patient IDs as values
                self._helper_tables["map_ids"] = map_dict
             

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = f"Formatting {self.get_cancer_type()} dataframes..."
        print(formatting_msg, end='\r')

        
        # Get dictionary with aliquots as keys and patient IDs as values
        mapping_dict = self._helper_tables["map_ids"]
        
        # Proteomics
        prot = self._data["proteomics"]
        prot['Patient_ID'] = prot['aliquot_submitter_id'].replace(mapping_dict) # replace aliquot with patient IDs 
        prot = prot.drop(['aliquot_submitter_id', 'case_submitter_id'], axis = 'columns')
        prot = prot.set_index('Patient_ID')
        prot.index = prot.index.str.replace('-T$','', regex = True)
        prot.index = prot.index.str.replace('-N$','.N', regex = True)  
        self._data["proteomics"] = prot

        # Phosphoproteomics 
        phos_df = self._data["phosphoproteomics"]
        phos_df['Patient_ID'] = phos_df['aliquot_submitter_id'].replace(mapping_dict) # replace aliquot with patient IDs 
        phos_df = phos_df.set_index('Patient_ID')
        phos_df = phos_df.drop(['aliquot_submitter_id', 'case_submitter_id'], axis = 'columns')
        phos_df.index = phos_df.index.str.replace('-T$','', regex = True)
        phos_df.index = phos_df.index.str.replace('-N$','.N', regex = True)  
        phos_df = map_database_to_gene_pdc(phos_df, 'refseq') # map refseq IDs to gene names
        self._data["phosphoproteomics"] = phos_df
        
        
        # Sort rows (tumor first then normal) and columns by first level (protein/gene name)
        self._data = sort_all_rows_pancan(self._data) 


        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
