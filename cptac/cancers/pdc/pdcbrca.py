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

class PdcBrca(Dataset):

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
                "acetylome.tsv.gz",
                "clinical.tsv.gz",
                "phosphoproteome.tsv.gz",
                "proteome.tsv.gz",
            ]
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="pdcbrca", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet, attempt_update_index=False)

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
                df = df.drop(['Internal Reference - Pooled Sample', 'RetroIR'])
                df = df.rename(index={'604':'CPT000814'}) # use the aliquot for 604
                self._data["clinical"] = df

            if file_name == "acetylome.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.set_index(["case_submitter_id", "aliquot_submitter_id"]) 
                df = df.drop('RetroIR', level = 1)
                self._data["acetylproteomics"] = df

            if file_name == "phosphoproteome.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.set_index(["case_submitter_id", "aliquot_submitter_id"])
                df = df.drop('RetroIR', level = 1)
                self._data["phosphoproteomics"] = df

            if file_name == "proteome.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.set_index(["case_submitter_id", "aliquot_submitter_id"])
                self._data["proteomics"] = df

                
        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = f"Formatting {self.get_cancer_type()} dataframes..."
        print(formatting_msg, end='\r')
        
        # prosp-brca-all-samples.txt shows which patient IDs have normal samples and which have replicates. This file
        # can be found on Box under CPTAC/cptac/pancan/helper_files 
        # 7 IDs with replicates: '11BR031', '11BR053', '11BR036', '11BR060', '14BR005', '11BR011', '21BR010'
        replicates = ['11BR031', '11BR053', '11BR036', '11BR060', '14BR005', '11BR011', '21BR010']
        
        # Drop normal aliquots 
        # There were 2 aliquots for 17 patients in the proteomics data. I calculated pearson correlations 
        # and created scatterplots to compare the values from each aliquot to its respective patient_ID in the 
        # original cptac data (flagship values). One aliquot for each patient did not correlate well 
        # (correlations between 0.001 and 0.4), while the other did (correlations between 0.7 and 0.9). 
        # Karsten Krug from the Broad suggested that the aliquots that did not correlate well are likely normal samples 
        # (18 patients had normal samples) which were dropped in downstream analysis because of quality control issues. 
        # The Patient_IDs with 2 aliquots did have a normal sample run (see mapping file 'prosp-brca-all-samples.txt').
        # Therefore, we drop them here. We kept the other aliquots that did correlate well. Of the 18 IDs with normal 
        # samples, 21BR010 only had one aliquot and it did not correlate well so it was dropped. I checked that the 
        # aliquots we dropped were normal samples using the biospecimen manifest for Brca proteomics on the PDC website 
        # under the biospecimen tab. I checked that the normal aliquots are the same for all omics, so we also dropped
        # them for phosphoproteomics and acetylproteomics. 
        # A file containing the correlations can be downloaded at: 
        # https://byu.box.com/shared/static/jzsq69bd079oq0zbicw4w616hyicd5ev.xlsx
        
        drop_normals = ['64ee175f-f3ce-446e-bbf4-9b6fa8_D1', '7ac27de9-0932-4ff5-aab8-29c527',
            '3208e021-1dae-42fd-bd36-0f3c3d', '6c660b6b-bfda-47b0-9499-160d49','241ecd0e-89bd-4d3a-81b3-55a250',
            '428de0d4-7f84-4075-bae1-352af6', '0a80d3c4-0758-447a-958c-ea868c', '53723086-8858-4395-93d7-0baa68',
            '1740224c-32d1-4c9f-98c6-653363', '885fe794-a98e-4f81-a284-ac4bb8', '4749ba99-d3b8-4ae3-b6f6-458bc7',
            '81116212-b7e6-454b-9579-105cf3', '1664b920-5e60-4e3b-9aab-fe121c', '3367406e-d39c-4641-a3e7-44e1f3',
            'e3d45dc6-66ef-4e0b-9d96-1b5db5', '33adae13-5dbd-4530-a5d5-3763e4', 'acf022b3-7f01-43b3-ac14-86f97d',
            '39f81c85-1832-45eb-829a-3040ad']
        
        # Proteomics
        prot = self._data["proteomics"]
        prot = prot.drop(drop_normals, level = 'aliquot_submitter_id') # drop normal aliquots (QC issues)
        prot = prot.rename(index={'604':'CPT000814'}) # use the aliquot for 604 
        prot.index = prot.index.droplevel('aliquot_submitter_id')
        self._data["proteomics"] = prot
        
        # Note -> At the time this code was added, phospho had all NaN values. In an email recieved 06/29/21, 
        # Paul Rudnick said the phospho study was broken. Whenever the phospho data is fixed, this 
        # code should work to drop normal aliquots and average replicates. It would be a good idea to check when 
        # the new data is ready though. 
        # Phosphoproteomics
        phos = self._data["phosphoproteomics"]
        phos = phos.drop(drop_normals, level = 'aliquot_submitter_id') # drop normal aliquots (QC issues)
        phos = phos.rename(index={'604':'CPT000814'}) # use the aliquot for 604
        phos.index = phos.index.droplevel('aliquot_submitter_id')
        phos = average_replicates(phos, id_list = replicates) # average replicates
        phos = map_database_to_gene_pdc(phos, 'refseq') # map refseq IDs to gene names
        self._data["phosphoproteomics"] = phos
        
        # Acetylproteomics
        acetyl = self._data["acetylproteomics"]
        acetyl = acetyl.drop(drop_normals, level = 'aliquot_submitter_id') # drop normal aliquots (QC issues)
        acetyl = acetyl.rename(index={'604':'CPT000814'}) # use the aliquot for 604
        acetyl.index = acetyl.index.droplevel('aliquot_submitter_id') # drop aliquots
        acetyl = average_replicates(acetyl, id_list = replicates) # average replicates
        acetyl = map_database_to_gene_pdc(acetyl, 'refseq') # map refseq IDs to gene names
        self._data["acetylproteomics"] = acetyl

        
        # Sort rows (tumor first then normal) and columns by first level (protein/gene name)
        self._data = sort_all_rows_pancan(self._data)  


        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
