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

class PdcLuad(Source):

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
                "mapping"              : "aliquot_to_patient_ID.tsv"
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
        super().__init__(cancer_type="luad", source='pdc', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)
        
        
    def load_acetylproteomics(self):
        df_type = 'acetylproteomics'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep='\t')
            
            # There were 2 duplicate IDs (same case ID and aliquot) in the phosphoproteomics and acetylproteomics data.
            # I used the mapping file "aliquot_to_patient_ID.tsv" to determine the tissue type for these duplicates. 
            # They were all tumor samples. Next, I ran a pearson correlation to check how well the values from each 
            # aliquot correlated to its respective flagship sample.  I also created a scatterplot for each aliquot 
            # and flagship pair. For phosphoproteomics, each aliquot correlated well with its flagship sample
            # so we averaged them. For acetylproteomics, the second occurrence of the duplicates correlated better 
            # with the flagship values, so we dropped the first ocurrence. 
            # A file containing the correlations can be downloaded at: 
            # https://byu.box.com/shared/static/jzsq69bd079oq0zbicw4w616hyicd5ev.xlsx

            # Common rows to drop
            drop_rows = ['Normal Only IR', 'Taiwanese IR', 'Tumor Only IR']

            # Get dictionary with aliquot_ID as keys and patient_ID as values
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]

            df['Patient_ID'] = df['aliquot_submitter_id'].replace(mapping_dict) # aliquots to patient IDs (normals have '.N')
            df = df.set_index('Patient_ID')
            df = df.drop(['aliquot_submitter_id', 'case_submitter_id'], axis = 'columns') 
            df = df.drop(drop_rows, axis = 'index')
            df = map_database_to_gene_pdc(df, 'refseq') # map refseq IDs to gene names
            df = rename_duplicate_labels(df, 'index') # add ".1" to the second ocurrence of the IDs with duplicates
            # drop 1st occurrence - didn't correlate as well with flagship
            df = df.drop(['C3N-02379', 'C3N-02587'], axis = 'index') 
            df = df.rename(index = {'C3N-02379.1':'C3N-02379', 'C3N-02587.1': 'C3N-02587'})# delete '.1' duplicate identifier 
            
            # save df in self._data
            self.save_df(df_type, df)
            
            
    def load_clinical(self):
        df_type = 'clinical'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep="\t", index_col=0)
            clin_drop_rows = ['Internal Reference - Pooled Sample', 'Normal Only IR', 'Taiwanese IR', 'Tumor Only IR']
            df = df.drop(clin_drop_rows, axis = 'index') # drop quality control and ref intensity
            
            # save df in self._data
            self.save_df(df_type, df)
            
            
    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep='\t')
            
            # There were 2 duplicate IDs (same case ID and aliquot) in the phosphoproteomics and acetylproteomics data.
            # I used the mapping file "aliquot_to_patient_ID.tsv" to determine the tissue type for these duplicates. 
            # They were all tumor samples. Next, I ran a pearson correlation to check how well the values from each 
            # aliquot correlated to its respective flagship sample.  I also created a scatterplot for each aliquot 
            # and flagship pair. For phosphoproteomics, each aliquot correlated well with its flagship sample
            # so we averaged them. For acetylproteomics, the second occurrence of the duplicates correlated better 
            # with the flagship values, so we dropped the first ocurrence. 
            # A file containing the correlations can be downloaded at: 
            # https://byu.box.com/shared/static/jzsq69bd079oq0zbicw4w616hyicd5ev.xlsx

            # Common rows to drop
            drop_rows = ['Normal Only IR', 'Taiwanese IR', 'Tumor Only IR']

            # Get dictionary with aliquot_ID as keys and patient_ID as values
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]

            df['Patient_ID'] = df['aliquot_submitter_id'].replace(mapping_dict) # map aliquots to patient IDs (normals have '.N')
            df = df.set_index('Patient_ID')
            df = df.drop(['aliquot_submitter_id', 'case_submitter_id'], axis = 'columns') 
            df = df.drop(drop_rows, axis = 'index')
            # these duplicates had high correlations with their tumor flagship samples, so we average them (see long comment above)
            df = average_replicates(df, ['C3N-02379', 'C3N-02587']) 
            df = map_database_to_gene_pdc(df, 'refseq') # map refseq IDs to gene names
            
            # save df in self._data
            self.save_df(df_type, df)
            
            
    def load_proteomics(self):
        df_type = 'proteomics'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep='\t')
            
            # There were 2 duplicate IDs (same case ID and aliquot) in the phosphoproteomics and acetylproteomics data.
            # I used the mapping file "aliquot_to_patient_ID.tsv" to determine the tissue type for these duplicates. 
            # They were all tumor samples. Next, I ran a pearson correlation to check how well the values from each 
            # aliquot correlated to its respective flagship sample.  I also created a scatterplot for each aliquot 
            # and flagship pair. For phosphoproteomics, each aliquot correlated well with its flagship sample
            # so we averaged them. For acetylproteomics, the second occurrence of the duplicates correlated better 
            # with the flagship values, so we dropped the first ocurrence. 
            # A file containing the correlations can be downloaded at: 
            # https://byu.box.com/shared/static/jzsq69bd079oq0zbicw4w616hyicd5ev.xlsx

            # Common rows to drop
            drop_rows = ['Normal Only IR', 'Taiwanese IR', 'Tumor Only IR']

            # Get dictionary with aliquot_ID as keys and patient_ID as values
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]

            df['Patient_ID'] = df['aliquot_submitter_id'].replace(mapping_dict) # map aliquots to patient IDs (normals have '.N')
            df = df.set_index('Patient_ID')
            df = df.drop(['aliquot_submitter_id', 'case_submitter_id'], axis = 'columns')
            
            # save df in self._data
            self.save_df(df_type, df)
            
            
    def load_mapping(self):
        df_type = 'mapping'
        
        if not self._helper_tables:
            file_path = self.locate_files(df_type)
            
            # aliquot_to_patient_ID.tsv contains only unique aliquots (no duplicates), 
            # so there is no need to slice out cancer specific aliquots
            df = pd.read_csv(file_path, sep = '\t', index_col = 'aliquot_ID', usecols = ['aliquot_ID', 'patient_ID'])
            map_dict = df.to_dict()['patient_ID'] # create dictionary with aliquot_ID as keys and patient_ID as values
            self._helper_tables["map_ids"] = map_dict
            