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

from cptac.cancers.source import Source
from cptac.tools.dataframe_tools import *

class BcmHnscc(Source):

    def __init__(self, version="latest", no_internet=False):
        """Define which bcmhnscc dataframes as are available in the self.load_functions dictionary variable, with names as keys.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest datafreeze. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        self.valid_versions = ["1.0"]

        self.data_files = {
            "1.0": {
                "circular_RNA" : "HNSCC-circRNA_rsem_tumor_normal_UQ_log2(x+1)_BCM.txt",
                "mapping" : "gencode.v34.basic.annotation-mapping.txt",
                "transcriptomics" : "HNSCC-gene_rsem_removed_circRNA_tumor_normal_UQ_log2(x+1)_BCM.txt"
            }
        }
        
        self.load_functions = {
            'circular_RNA' : self.load_circular_RNA,
            'transcriptomics' : self.load_transcriptomics,
        }
        
        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        # Call the parent class __init__ function
        super().__init__(cancer_type="hnscc", source='bcm', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

        
    def load_circular_RNA(self):
        df_type = 'circular_RNA'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep="\t")
            df = df.rename_axis('INDEX').reset_index()
            df[["circ","chrom","start","end","gene"]] = df.INDEX.str.split('_', expand=True)
            df["circ_chromosome"] = df["circ"] +"_" + df["chrom"]
            df = df.set_index('gene')
            
            # Add gene names to circular RNA data
            self.load_mapping()
            gene_key = self._helper_tables["gene_key"]
            df = gene_key.join(df, how = "inner")
            df = df.reset_index()
            df = df.rename(columns= {"gene_name": "Name","gene":"Database_ID"}) # change names to match cptac package
            df = df.set_index(["Name","circ_chromosome", "start","end","Database_ID"]) #create multi-index
            df.drop(['INDEX', 'circ', 'chrom'], axis=1, inplace=True) 
            df = df.sort_index()
            df = df.T
            df.index = df.index.str.replace(r"_T", "", regex=True) # remove Tumor label
            df.index = df.index.str.replace(r"_A", ".N", regex=True)# Normal samples labeled with .N
            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)
            
        
    def load_mapping(self):
        df_type = 'mapping'
        # self._helper_tables is a dictionary of helpful dataframes that the user does not need to access
        # dataframes here are used to load the other data types, but don't show up when the user lists available data
        # this way mapping only needs to be loaded once and all other types can use it when they are loaded
        if not self._helper_tables:
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep="\t")
            df = df[["gene","gene_name"]] #only need gene (database gene id) and gene_name (common gene name)
            df = df.set_index("gene")
            df = df.drop_duplicates()
            self._helper_tables["gene_key"] = df
            
        
    def load_transcriptomics(self):
        df_type = 'transcriptomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep="\t")
            df.index.name = 'gene'

            # Add gene names to transcriptomic data
            self.load_mapping()
            gene_key = self._helper_tables["gene_key"]
            transcript = gene_key.join(df, how = "inner") #keep only gene_ids with gene names
            transcript = transcript.reset_index()
            transcript = transcript.rename(columns={"gene_name":"Name","gene":"Database_ID"})
            transcript = transcript.set_index(["Name", "Database_ID"])
            transcript = transcript.sort_index() #alphabetize
            transcript = transcript.T
            transcript.index = transcript.index.str.replace(r"_T", "", regex=True)
            transcript.index = transcript.index.str.replace(r"_A", ".N", regex=True)# Normal samples labeled with .N
            transcript.index.name = "Patient_ID"

            df = transcript
            # save df in self._data
            self.save_df(df_type, df)