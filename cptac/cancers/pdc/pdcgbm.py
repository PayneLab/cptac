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
import mygene

from cptac.cancers.source import Source
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError

class PdcGbm(Source):

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
                "acetylproteomics"     : "acetylproteomics.tsv.gz",
                "clinical"             : "clinical.tsv.gz",
                "phosphoproteomics"    : "phosphoproteomics.tsv.gz",
                "proteomics"           : "proteomics.tsv.gz",
                "mapping"              : "GBM_normal_sample_mapping.xlsx"
            }
        }
        
        self.load_functions = {
            'acetylproteomics' : self.load_acetylproteomics,
            'clinical' : self.load_clinical,
            'phosphoproteomics' : self.load_phosphoproteomics,
            'proteomics' : self.load_proteomics,
        }
        
        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        # Call the parent class __init__ function
        super().__init__(cancer_type="gbm", source='pdc', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)
        
        
    def load_acetylproteomics(self):
        df_type = 'acetylproteomics'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep='\t')
            
            # Get dictionary to map GTEX IDs to patient IDs (for GTEX normal samples)
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]
            
            df['Patient_ID'] = df['case_submitter_id'].replace(mapping_dict) # GTEX IDs to patient IDs for normal samples
            df = df.set_index('Patient_ID')
            df = df.drop(['aliquot_submitter_id', 'case_submitter_id'], axis = 'columns')
            df = map_database_to_gene_pdc(df, 'refseq') # map refseq IDs to gene names
            
            # save df in self._data
            self.save_df(df_type, df)
            
            
    def load_clinical(self):
        df_type = 'clinical'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.drop(['ref'], axis = 'index') # Drop quality control and ref intensity
            
            # Get dictionary to map GTEX IDs to patient IDs (for GTEX normal samples)
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]
            
            df = df.reset_index()
            df['Patient_ID'] = df['case_submitter_id'].replace(mapping_dict) # GTEX IDs to patient IDs for normal samples
            df = df.set_index('Patient_ID')
            df = df.drop(['case_submitter_id'], axis = 'columns') 
            
            # save df in self._data
            self.save_df(df_type, df)
            
            
    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep='\t')
            
            # Get dictionary to map GTEX IDs to patient IDs (for GTEX normal samples)
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]
            
            df['Patient_ID'] = df['case_submitter_id'].replace(mapping_dict) # GTEX IDs to patient IDs for normal samples
            df = df.set_index('Patient_ID')
            df = df.drop(['aliquot_submitter_id', 'case_submitter_id'], axis = 'columns')
            df = map_database_to_gene_pdc(df, 'refseq') # map refseq IDs to gene names
            
            # save df in self._data
            self.save_df(df_type, df)
            
            
    def load_proteomics(self):
        df_type = 'proteomics'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep='\t')
            
            # Get dictionary to map GTEX IDs to patient IDs (for GTEX normal samples)
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]
            
            # Proteomics
            df['Patient_ID'] = df['case_submitter_id'].replace(mapping_dict) # GTEX IDs to patient IDs for normal samples
            df = df.set_index('Patient_ID')
            df = df.drop(['aliquot_submitter_id', 'case_submitter_id'], axis = 'columns')
            
            # save df in self._data
            self.save_df(df_type, df)
            
    
    def load_mapping(self):
        df_type = 'mapping'
        
        if not self._helper_tables:
            file_path = self.locate_files(df_type)

            # This file maps GTEX normal samples to our desired 'PT-' subject identifier and
            # can be found on Box under CPTAC/cptac/pancan/helper_files
            df = pd.read_excel(file_path, index_col = 'Original Id', usecols = ['Original Id', 'Subject ID'])
            df['Subject ID'] = df['Subject ID'].apply(lambda x: x+'.N' if 'PT-' in x else x) # add normal identifier 
            map_dict = df.to_dict()['Subject ID']
            self._helper_tables["map_ids"] = map_dict
            