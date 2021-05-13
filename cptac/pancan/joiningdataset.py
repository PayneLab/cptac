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
from cptac.file_tools import *
from cptac.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError


class JoiningDataset(Dataset):

    def __init__(self, input_data):
        """Load all of the dataframes as values in the self._data dict variable, with names as keys, and format them properly.
        
        Parameters:
        input_data (dict of dict): Data dictionaries to combine for joining
        """

        # Call the parent class __init__ function
        super().__init__(cancer_type='joining', version=None, valid_versions=None, data_files=None, no_internet=None, skip_init=True)

        # Put all the input data tables into self._data
        # Each table's name needs to be its source followed by the data type, e.g. the proteomics table from bcm would be "bcm_transcriptomics"
               
        for source in input_data.keys():
            source_data = input_data[source]
            for data_type in source_data.keys():   
                self._data[source + "_" + data_type] = input_data[source][data_type]
        
        # create list of valid datasets ex bcm_transcriptomics, bcm_proteomics ect and put into self._valid_omics_dfs
        valid = []
        for source in input_data.keys():
            for df in self._valid_omics_dfs:
                valid.append(f"{source}_{df}")
                
        self._valid_omics_dfs = valid
        
        # create list of valid meta data sets ex: bcm_derived_molecular and put into self._valid_metadata_dfs
        meta_valid = []
        for source in input_data.keys():
            for df in self._valid_metadata_dfs:
                meta_valid.append(f"{source}_{df}")
                
        self._valid_metadata_dfs = meta_valid
        
  # for the join function to work the clincal data needs to be named clincal not mssm_clinical  
        self._data["clinical"] = self._data["mssm_clinical"]
        self._data["somatic_mutation"] = self._data["harmonized_somatic_mutation"]
        
        #create unionized indices for clinical
        
        master_index = unionize_indices(self._data).dropna()
        new_clinical = self._data["clinical"]
        new_clinical = new_clinical.reindex(master_index)
       
        sample_status_col = generate_sample_status_col(new_clinical, normal_test=lambda sample: sample.endswith('.N') if type(sample) is str else sample[0].endswith(".N"))
        new_clinical = new_clinical.drop(columns=['Sample_Tumor_Normal'])
        
        new_clinical.insert(0, "Sample_Tumor_Normal", sample_status_col)
        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = new_clinical

        