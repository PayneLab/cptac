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

from cptac.cancer import Cancer
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError

class PdcPdac(Cancer):

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
                "clinical.tsv.gz", # error with download function
                "phosphoproteome.tsv.gz",
                "proteome.tsv.gz",
                "aliquot_to_patient_ID.tsv"
            ]
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="pdcpdac", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet, attempt_update_index=False)

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
                self._data["clinical"] = df

            if file_name == "phosphoproteome.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                self._data["phosphoproteomics"] = df

            if file_name == "proteome.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                self._data["proteomics"] = df
                
            # aliquot_to_patient_ID.tsv contains only unique aliquots (no duplicates), 
            # so there is no need to slice out cancer specific aliquots
            elif file_name == "aliquot_to_patient_ID.tsv":
                df = pd.read_csv(file_path, sep = "\t", index_col = 'aliquot_ID', usecols = ['aliquot_ID', 'patient_ID'])
                map_dict = df.to_dict()['patient_ID'] # create dictionary with aliquots as keys and patient IDs as values
                self._helper_tables["map_ids"] = map_dict

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = f"Formatting {self.get_cancer_type()} dataframes..."
        print(formatting_msg, end='\r')
        
        
        # Quality control and reference rows 
        drop_rows = ['KoreanReference1', 'KoreanReference2', 'KoreanReference3',
             'QC1', 'QC2', 'QC3', 'QC4', 'QC5', 'QC6', 'WU-PDA1']
        
        # These 8 aliquots were not in the mapping file. Yize said they are all normal samples.
        manually_mapped = {'CPT0347760002': 'C3L-07032.N', 'CPT0347790002': 'C3L-07033.N',
            'CPT0347820002': 'C3L-07034.N', 'CPT0347850002': 'C3L-07035.N', 'CPT0347880002': 'C3L-07036.N',
            'CPT0355180003': 'C3L-03513.N', 'CPT0355190003': 'C3L-03515.N', 'CPT0355200003': 'C3L-03514.N'}
        
        # Get dictionary with aliquots as keys and patient IDs as values
        mapping_dict = self._helper_tables["map_ids"]
        
        # Clinical
        clin = self._data["clinical"]
        clin = clin.drop(drop_rows + ['WU-pooled sample', 'pooled sample'], axis = 'index') # drop quality control and references
        self._data["clinical"] = clin
        
        # Proteomics
        prot = self._data["proteomics"]        
        prot['Patient_ID'] = prot['aliquot_submitter_id'].replace(mapping_dict) # aliquots to patient IDs (normals have '.N')
        prot = prot.set_index('Patient_ID')
        prot = prot.rename(index = manually_mapped) # map 8 aliquots that were not in the mapping file
        prot = prot.drop(['aliquot_submitter_id', 'case_submitter_id'], axis = 'columns')
        prot = prot.drop(drop_rows, axis = 'index') # drop quality control and references
        self._data["proteomics"] = prot
        
        # Phosphoproteomics
        phos = self._data["phosphoproteomics"]
        phos['Patient_ID'] = phos['aliquot_submitter_id'].replace(mapping_dict) # aliquots to patient IDs (normals have '.N')
        phos = phos.set_index('Patient_ID')
        phos = phos.rename(index = manually_mapped) # map 8 aliquots that were not in the mapping file
        phos = phos.drop(['aliquot_submitter_id', 'case_submitter_id'], axis = 'columns')        
        phos = phos.drop(drop_rows + ['WU-pooled sample', 'pooled sample'], axis = 'index') # drop quality control and references
        phos = map_database_to_gene_pdc(phos, 'refseq') # map refseq IDs to gene names
        self._data["phosphoproteomics"] = phos
        
        
        # Sort rows (tumor first then normal) and columns by first level (protein/gene name)
        self._data = sort_all_rows_pancan(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
