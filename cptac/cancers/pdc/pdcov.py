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

from cptac.cancers.source import Source
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError

class PdcOv(Source):

    def __init__(self, version="latest", no_internet=False):
        """Define which dataframes as are available in the self.load_functions dictionary variable, with names as keys.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest datafreeze. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        self.valid_versions = ["1.0"]

        self.data_files = {
            "1.0": {
                "clinical"             : "clinical.tsv.gz",
                "phosphoproteomics"    : "phosphoproteomics.tsv.gz",
                "proteomics"           : "proteomics.tsv.gz",
                "mapping"              : "OV_sample_TMT_annotation_UMich_GENCODE34_0315.csv"
            }
        }
        
        self.load_functions = {
            'clinical' : self.load_clinical,
            'phosphoproteomics' : self.load_phosphoproteomics,
            'proteomics' : self.load_proteomics,
        }
        
        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        # Call the parent class __init__ function
        super().__init__(cancer_type="ov", source='pdc', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)
        
        
    def load_clinical(self):
        df_type = 'clinical'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep='\t', index_col=0)
            clin_drop_rows = ['JHU QC', 'PNNL-JHU Ref']
            df = df.drop(clin_drop_rows, axis = 'index') # drop quality control and reference
            
            # save df in self._data
            self.save_df(df_type, df)
            
            
    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep='\t')
            
            # Get dictionary with aliquots as keys and patient IDs as values
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]
            
            df_df['Patient_ID'] = df_df['aliquot_submitter_id'].replace(mapping_dict) # replace aliquot with patient IDs 
            df_df = df_df.set_index('Patient_ID')
            df_df = df_df.drop(['aliquot_submitter_id', 'case_submitter_id'], axis = 'columns')
            df_df.index = df_df.index.str.replace('-T$','', regex = True)
            df_df.index = df_df.index.str.replace('-N$','.N', regex = True)  
            df_df = map_database_to_gene_pdc(df_df, 'refseq') # map refseq IDs to gene names
            self._data["phosphoproteomics"] = phos_df
            
            # save df in self._data
            self.save_df(df_type, df)
            
            
    def load_proteomics(self):
        df_type = 'proteomics'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep='\t')
            
            # Get dictionary with aliquots as keys and patient IDs as values
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]

            df['Patient_ID'] = df['aliquot_submitter_id'].replace(mapping_dict) # replace aliquot with patient IDs 
            df = df.drop(['aliquot_submitter_id', 'case_submitter_id'], axis = 'columns')
            df = df.set_index('Patient_ID')
            df.index = df.index.str.replace('-T$','', regex = True)
            df.index = df.index.str.replace('-N$','.N', regex = True)  
            self._data["proteomics"] = prot
            
            # save df in self._data
            self.save_df(df_type, df)
            
            
    def load_mapping(self):
        df_type = 'mapping'
        
        if not self._helper_tables:
            file_path = self.locate_files(df_type)
            
            # This file maps Ov aliquots to patient IDs (case ID with tissue type) and 
            # can be found on Box under CPTAC/cptac/pancan/helper_files
            ov_map = pd.read_csv(file_path, sep = ",", usecols = ['specimen', 'sample'])
            ov_map = ov_map.loc[~ ov_map['sample'].str.contains('JHU', regex = True)] # drop quality control rows
            ov_map = ov_map.set_index('specimen')
            map_dict = ov_map.to_dict()['sample'] # create dictionary with aliquots as keys and patient IDs as values
            self._helper_tables["map_ids"] = map_dict
