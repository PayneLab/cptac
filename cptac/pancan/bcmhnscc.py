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


class BcmHnscc(Dataset):

    def __init__(self, no_internet, version):
        """Load all of the bcmbrca dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        valid_versions = ["1.0"]

        data_files = {
            "1.0": [
                "HNSCC-circRNA_rsem_tumor_normal_UQ_log2(x+1)_BCM.txt",
                "gencode.v34.basic.annotation-mapping.txt",
                "HNSCC-gene_rsem_removed_circRNA_tumor_normal_UQ_log2(x+1)_BCM.txt"
            ]
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="bcmhnscc", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

            if file_name == "HNSCC-gene_rsem_removed_circRNA_tumor_normal_UQ_log2(x+1)_BCM.txt":
                df = pd.read_csv(file_path, sep="\t")
                df.index.name = 'gene'
                self._data["transcriptomics"] = df
                
            if file_name == "gencode.v34.basic.annotation-mapping.txt":
                df = pd.read_csv(file_path, sep="\t")
                df = df[["gene","gene_name"]] #only need gene (database gene id) and gene_name (common gene name)
                df = df.set_index("gene")
                df = df.drop_duplicates()
                self._helper_tables["gene_key"] = df 
            
            if file_name == "HNSCC-circRNA_rsem_tumor_normal_UQ_log2(x+1)_BCM.txt":
                df = pd.read_csv(file_path, sep="\t")
                df = df.rename_axis('INDEX').reset_index()
                df[["circ","chrom","start","end","gene"]] = df.INDEX.str.split('_', expand=True)
                df["circ_chromosome"] = df["circ"] +"_" + df["chrom"]
                df = df.set_index('gene')
                self._data["circular_RNA"] = df
    
            
           

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        
        # Add gene names to transcriptomic data 
        
        prot = self._data["transcriptomics"]
        gene_key = self._helper_tables["gene_key"]
        transcript = gene_key.join(prot,how = "inner") #keep only gene_ids with gene names
        transcript = transcript.reset_index()
        transcript = transcript.rename(columns={
         "gene_name":"Name","gene":"Database_ID"})
        transcript = transcript.set_index(["Name", "Database_ID"])
        transcript = transcript.sort_index() #alphabetize
        transcript = transcript.T
        transcript.index = transcript.index.str.replace(r"_T", "", regex=True)
        transcript.index = transcript.index.str.replace(r"_A", ".N", regex=True)# Normal samples labeled with .N
        transcript.index.name = "Patient_ID"
        
          # Sort values
        normal = transcript.loc[transcript.index.str.contains('\.N$', regex = True)]
        normal = normal.sort_values(by=["Patient_ID"])
        tumor = transcript.loc[~ transcript.index.str.contains('\.N$', regex = True)]
        tumor = tumor.sort_values(by=["Patient_ID"])
        all_transcript = tumor.append(normal)  
        
        self._data["transcriptomics"] = all_transcript
        
        # Add gene names to circular RNA data 
        circRNA = self._data["circular_RNA"]
        
        df = gene_key.join(circRNA, how = "inner")
        df = df.reset_index()
        df = df.rename(columns= {"gene_name": "Name","gene":"Database_ID"}) # change names to match cptac package
        df = df.set_index(["Name","circ_chromosome", "start","end","Database_ID"]) #create multi-index
        df.drop(['INDEX', 'circ', 
                'chrom'], axis=1, inplace=True) 
        df = df.sort_index()
        df = df.T
        df.index = df.index.str.replace(r"_T", "", regex=True) # remove Tumor label
        df.index = df.index.str.replace(r"_A", ".N", regex=True)# Normal samples labeled with .N
        df.index.name = "Patient_ID"
        
        # Sort values
        normal = df.loc[df.index.str.contains('\.N$', regex = True)]
        normal = normal.sort_values(by=["Patient_ID"])
        tumor = df.loc[~ df.index.str.contains('\.N$', regex = True)]
        tumor = tumor.sort_values(by=["Patient_ID"])
        all_prot = tumor.append(normal)  
        
        self._data["circular_RNA"] = all_prot
      


       

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
