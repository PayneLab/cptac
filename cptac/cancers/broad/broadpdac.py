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
from gtfparse import read_gtf

from cptac.cancers.source import Source
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError


class BroadPdac(Source):

    def __init__(self, version="latest", no_internet=False):
        """Define which broadpdac dataframes as are available in the self.load_functions dictionary variable, with names as keys.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest datafreeze. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        self.valid_versions = ["1.0"]

        self.data_files = {
            "1.0": {
                "transcriptomics" : "PDAC.rsem_transcripts_tpm.txt.gz",
                "mapping" : ["sample_descriptions.tsv", "gencode.v34.GRCh38.genes.collapsed_only.gtf"]
            }
        }
        
        self.load_functions = {
            'transcriptomics' : self.load_transcriptomics,
        }
        
        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        # Call the parent class __init__ function
        super().__init__(cancer_type="pdac", source='broad', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

        
    def load_mapping(self):
        df_type = 'mapping'
        
        # Since this is the only location where things are added to _helper_tables, just check if they are empty
        # If they are empty, populate them
        if not self._helper_tables:
            file_path_list = self.locate_files(df_type)
            for file_path in file_path_list:
                path_elements = file_path.split(os.sep) # Get a list of the levels of the path
                file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below
                
                if file_name == "sample_descriptions.tsv":
                    broad_key = pd.read_csv(file_path, sep="\t")
                    broad_key = broad_key.loc[broad_key['cohort'] == "PDAC"] #get only PDAC keys
                    broad_key = broad_key[["sample_id","GDC_id","tissue_type"]]
                    broad_key = broad_key.set_index("sample_id")#set broad id as index
                    #add tumor type identification to end
                    broad_key["Patient_ID"] = broad_key["GDC_id"] + broad_key["tissue_type"] 
                    #change so tumor samples have nothing on end of id and .N for normal samples
                    broad_key.Patient_ID = broad_key.Patient_ID.str.replace(r"Tumor", "", regex=True)
                    broad_key.Patient_ID = broad_key.Patient_ID.str.replace(r"Normal", ".N", regex=True)
                    #covert df to dictionary
                    broad_dict = broad_key.to_dict()["Patient_ID"]
                    self._helper_tables["broad_key"] = broad_dict
                    
                elif file_name == "gencode.v34.GRCh38.genes.collapsed_only.gtf":
                    broad_gene_names = read_gtf(file_path)
                    broad_gene_names = broad_gene_names[["gene_name","gene_id"]]
                    broad_gene_names = broad_gene_names.rename(columns= {"gene_name":"Name"}) #change name to merge 
                    broad_gene_names = broad_gene_names.set_index("gene_id")
                    broad_gene_names = broad_gene_names.drop_duplicates()
                    self._helper_tables["broad_gene_names"] = broad_gene_names
        
        
    def load_transcriptomics(self):
        df_type = 'transcriptomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep="\t")
            df = df.set_index(["transcript_id","gene_id"])
            
            # Add gene names to transcriptomic data
            self.load_mapping()
            broad_gene_names = self._helper_tables["broad_gene_names"]
            broad_dict = self._helper_tables["broad_key"]        
            df = broad_gene_names.join(df, how = "left") #merge in gene names keep transcripts that have a gene name
            df = df.reset_index()
            df = df.rename(columns= {"transcript_id": "Transcript_ID","gene_id":"Database_ID"})
            df = df.set_index(["Name","Transcript_ID","Database_ID"])
            df = df.rename(columns = broad_dict)# rename columns with CPTAC IDs
            df = df.sort_index() 
            df = df.T
            df.index.name = "Patient_ID"
            # save df in self._data
            self.save_df(df_type, df)