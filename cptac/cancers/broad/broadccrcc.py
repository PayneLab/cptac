# Copyright 2018 Samuel Payne sam_payne@byu.edu
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pandas as pd
import os
from pyranges import read_gtf
from cptac.cancers.source import Source

class BroadCcrcc(Source):
    def __init__(self, no_internet=False):
        """
        Initialize the BroadCcrcc object.

        Parameters:
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. 
        This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. 
        Default is False.
        """
        # Set some needed variables, and pass them to the parent Dataset class __init__ function
        self.data_files = {
            "transcriptomics" : "CCRCC.rsem_transcripts_tpm.txt.gz",
            "mapping" : ["sample_descriptions.tsv.gz", "gencode.v34.GRCh38.genes.collapsed_only.gtf.gz"]
        }
        
        self.load_functions = {
            'transcriptomics' : self.load_transcriptomics,
        }
        
        # Call the parent class __init__ function
        super().__init__(cancer_type="ccrcc", source='broad', data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

    def load_mapping(self):
        """
        Load mapping from the provided files.
        """
        df_type = 'mapping'
        
        # If _helper_tables are empty, populate them
        if not self._helper_tables:
            file_path_list = self.locate_files(df_type)
            for file_path in file_path_list:
                path_elements = file_path.split('/') # Get a list of the levels of the path
                file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below
                
                if file_name == "sample_descriptions.tsv.gz":
                    # Load broad keys specific for CCRCC
                    broad_key = pd.read_csv(file_path, sep="\t")
                    broad_key = broad_key.loc[broad_key['cohort'] == "CCRCC"]
                    broad_key = broad_key[["sample_id","GDC_id","tissue_type"]]
                    broad_key = broad_key.set_index("sample_id") # set broad id as index
                    
                    # Add tumor type identification to end and standardize patient id format
                    broad_key["Patient_ID"] = broad_key["GDC_id"] + broad_key["tissue_type"]
                    broad_key.Patient_ID = broad_key.Patient_ID.str.replace(r"Tumor", "", regex=True)
                    broad_key.Patient_ID = broad_key.Patient_ID.str.replace(r"Normal", ".N", regex=True)
                    
                    # Convert df to dictionary
                    broad_dict = broad_key.to_dict()["Patient_ID"]
                    self._helper_tables["broad_key"] = broad_dict
                    
                elif file_name == "gencode.v34.GRCh38.genes.collapsed_only.gtf.gz":
                    # Load broad gene names
                    broad_gene_names = read_gtf(file_path)
                    broad_gene_names = broad_gene_names.as_df()
                    broad_gene_names = broad_gene_names[["gene_name","gene_id"]]
                    broad_gene_names = broad_gene_names.rename(columns = {"gene_name":"Name"}) #change name to merge 
                    broad_gene_names = broad_gene_names.set_index("gene_id")
                    broad_gene_names = broad_gene_names.drop_duplicates()                
                    self._helper_tables["broad_gene_names"] = broad_gene_names
        
    def load_transcriptomics(self):
        """
        Load transcriptomics data from the provided file.
        """
        df_type = 'transcriptomics'

        if df_type not in self._data:
            # Perform initial checks and get file path
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep="\t")
            df = df.set_index(["transcript_id","gene_id"])
            
            # Add gene names to transcriptomic data
            self.load_mapping()
            broad_gene_names = self._helper_tables["broad_gene_names"]
            broad_dict = self._helper_tables["broad_key"]        
            df = broad_gene_names.join(df, how = "left") # Merge in gene names, keep transcripts that have a gene name
            df = df.reset_index()
            df = df.rename(columns= {"transcript_id": "Transcript_ID","gene_id":"Database_ID"})
            df = df.set_index(["Name","Transcript_ID","Database_ID"])
            df = df.rename(columns = broad_dict) # Rename columns with CPTAC IDs
            df = df.sort_index() 
            df = df.T
            df.index.name = "Patient_ID"
            
            # Save df in self._data
            self.save_df(df_type, df)
