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
from pyranges import read_gtf
from cptac.cancers.source import Source

class BroadLuad(Source):
    def __init__(self, no_internet=False):
        """Defines the available dataframes in the self.load_functions dictionary variable with names as keys.

        Parameters:
        no_internet (bool, optional): If True, the index update step which requires an internet connection is skipped. Default is False.
        """
        
        # Initialize necessary variables and pass them to the parent Dataset class __init__ function
        self.data_files = {
            "transcriptomics" : "LUAD.rsem_transcripts_tpm.txt.gz",
            "mapping" : ["sample_descriptions.tsv.gz", "gencode.v34.GRCh38.genes.collapsed_only.gtf.gz", "aliquot_to_patient_ID.tsv.gz"]
        }

        self.load_functions = {
            'transcriptomics' : self.load_transcriptomics,
        }

        super().__init__(cancer_type="luad", source='broad', data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

    def load_mapping(self):
        """Loads and preprocesses the mapping files. Generates helper tables for later use."""
        df_type = 'mapping'
        
        # If _helper_tables is empty, populate them
        if not self._helper_tables:
            file_path_list = self.locate_files(df_type)
            for file_path in file_path_list:
                file_name = file_path.split('/')[-1] # The last element will be the name of the file
                
                # Sample descriptions mapping
                if file_name == "sample_descriptions.tsv.gz":
                    broad_key = pd.read_csv(file_path, sep="\t")
                    broad_key = broad_key.loc[broad_key['cohort'] == "LUAD"] #get only LUAD keys
                    broad_key = broad_key[["sample_id","GDC_id","tissue_type"]]
                    broad_key = broad_key.set_index("sample_id")#set broad id as index
                    #key only the main Patient IDs (first 9 characters) and remove the -0N aliqout IDs
                    broad_key['GDC_id'] = broad_key['GDC_id'].str[:9]
                    #add tumor type identification to end
                    broad_key["Patient_ID"] = broad_key["GDC_id"] + broad_key["tissue_type"] 
                    #change so tumor samples have nothing on end of id and .N for normal samples
                    broad_key.Patient_ID = broad_key.Patient_ID.str.replace(r"Tumor", "", regex=True)
                    broad_key.Patient_ID = broad_key.Patient_ID.str.replace(r"Normal", ".N", regex=True)
                    #covert df to dictionary
                    broad_dict = broad_key.to_dict()["Patient_ID"]
                    self._helper_tables["broad_key"] = broad_dict
                    
                # Gene name mapping   
                elif file_name == "gencode.v34.GRCh38.genes.collapsed_only.gtf.gz":
                    broad_gene_names = read_gtf(file_path)
                    broad_gene_names = broad_gene_names.as_df()
                    broad_gene_names = broad_gene_names[["gene_name","gene_id"]]
                    broad_gene_names = broad_gene_names.rename(columns= {"gene_name":"Name"}) #change name to merge 
                    broad_gene_names = broad_gene_names.set_index("gene_id")
                    broad_gene_names = broad_gene_names.drop_duplicates()
                    self._helper_tables["broad_gene_names"] = broad_gene_names
                    
                # Aliquot to patient ID mapping    
                elif file_name == "aliquot_to_patient_ID.tsv.gz":
                    df = pd.read_csv(file_path, sep = "\t", index_col = 0)
                    self._helper_tables["map_ids"] = df

    def load_transcriptomics(self):
        """Loads and preprocesses the transcriptomics data."""
        df_type = 'transcriptomics'

        if df_type not in self._data:
            # Perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep="\t")
            df = df.set_index(["transcript_id","gene_id"])
            
            # Add gene names to transcriptomic data
            self.load_mapping()
            broad_gene_names = self._helper_tables["broad_gene_names"]
            broad_dict = self._helper_tables["broad_key"]
            mapping_df = self._helper_tables["map_ids"]
            aliquot_dict = mapping_df.to_dict()["patient_ID"]        
            df = broad_gene_names.join(df,how = "left") #merge in gene names keep transcripts that have a gene name
            df = df.reset_index()
            df = df.rename(columns= {"transcript_id": "Transcript_ID","gene_id":"Database_ID"})
            df = df.set_index(["Name","Transcript_ID","Database_ID"])
            df = df.rename(columns = broad_dict)# rename columns with CPTAC IDs
            df = df.rename(columns = aliquot_dict) 
            df = df.T
            # average duplicates: C3N-00545-03 and C3N-00545-01 were seperated out as two different aliqout ids. 
            # They are from the same sample, so we average them. 
            df = df.groupby("index", level = 0).mean() 
            df.index.name = "Patient_ID"     
            # save df in self._data
            self.save_df(df_type, df)
