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


class Harmonized(Dataset):

    def __init__(self, no_internet, version, filter_type): 
        """Load all of the mssmclinical dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        valid_versions = ["1.0"]

        data_files = {
            "1.0": [
                "PanCan_Union_Maf_Broad_WashU.maf"                
            ]
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type='harmonized', version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet) # changed 'mssmclinical' to cancer_type

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below
            
            # Get tumor_code
            tumor_codes = {'pancanbrca': 'BRCA', 'pancanccrcc':'CCRCC', 
                           'pancanucec':'UCEC','pancangbm':'GBM','pancanhnscc':'HNSCC',
                           'pancanlscc': 'LSCC','pancanluad':'LUAD', 'pancanpda':'PDA',
                           'pancanhcc':'HCC','pancancoad':'CO','pancanov':'OV'}

            if file_name == "PanCan_Union_Maf_Broad_WashU.maf":
                df = pd.read_csv(file_path, sep="\t", low_memory = False)
                df = df.loc[df['COHORT'] == tumor_codes[filter_type]] 
                df['Patient_ID'] = df.loc[:, 'Tumor_Sample_Barcode']
                df = df.rename(columns={
                         "Hugo_Symbol":"Gene",
                         "Variant_Classification":"Mutation",
                         "Protein_Change":"Location"})

                df = df.set_index("Patient_ID")
                df = df[ ['Gene'] + ["Mutation"] + ["Location"] + [ col for col in df.columns if col not in ["Gene","Mutation","Location"] ] ]
                df.index = df.index.str.replace(r"_T", "", regex=True) # data based on Tumor and Normal. Remove _T 
                self._data["somatic_mutation"] = df                      
        
                
        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')
       
        

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
