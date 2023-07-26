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
import os
from pyranges import read_gtf
from cptac.cancers.source import Source

class BroadBrca(Source):
    def __init__(self, no_internet=False):
        """Initialization of BroadBrca class.

        Args:
            no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. Default is False.
        """
        self.data_files = {
            "transcriptomics" : "BRCA.rsem_transcripts_tpm.txt.gz",
            "mapping" : ["sample_descriptions.tsv.gz", "gencode.v34.GRCh38.genes.collapsed_only.gtf.gz"]
        }

        self.load_functions = {
            'transcriptomics' : self.load_transcriptomics,
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="brca", source='broad', data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

    def load_mapping(self):
        """Load mapping file."""

        df_type = 'mapping'

        # Check if helper tables are empty, if so, populate them
        if not self._helper_tables:
            file_path_list = self.locate_files(df_type)
            for file_path in file_path_list:
                path_elements = file_path.split('/') # Split path into components
                file_name = path_elements[-1] # Last element will be the file name

                # Parse files based on file name
                if file_name == "sample_descriptions.tsv.gz":
                    # Processing for the sample descriptions file
                    self._process_sample_descriptions(file_path)
                    
                elif file_name == "gencode.v34.GRCh38.genes.collapsed_only.gtf.gz":
                    # Processing for the gtf file
                    self._process_gtf(file_path)

    def _process_sample_descriptions(self, file_path):
        """Process the sample descriptions file."""
        broad_key = pd.read_csv(file_path, sep="\t")
        broad_key = broad_key.loc[broad_key['cohort'] == "BRCA"]
        broad_key = broad_key[["sample_id","GDC_id","tissue_type"]].set_index("sample_id")
        broad_key["Patient_ID"] = broad_key["GDC_id"] + broad_key["tissue_type"] 
        broad_key.Patient_ID = broad_key.Patient_ID.str.replace(r"Tumor", "", regex=True)
        broad_key.Patient_ID = broad_key.Patient_ID.str.replace(r"Normal", ".N", regex=True)
        broad_dict = broad_key.to_dict()["Patient_ID"]
        self._helper_tables["broad_key"] = broad_dict

    def _process_gtf(self, file_path):
        """Process the gtf file."""
        broad_gene_names = read_gtf(file_path)
        broad_gene_names = broad_gene_names.as_df()
        broad_gene_names = broad_gene_names[["gene_name","gene_id"]].rename(columns = {"gene_name":"Name"}).set_index("gene_id")
        broad_gene_names = broad_gene_names.drop_duplicates()                
        self._helper_tables["broad_gene_names"] = broad_gene_names

    def load_transcriptomics(self):
        """Load transcriptomics data."""
        df_type = 'transcriptomics'

        if df_type not in self._data:
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep="\t").set_index(["transcript_id","gene_id"])

            self.load_mapping()
            broad_gene_names = self._helper_tables["broad_gene_names"]
            broad_dict = self._helper_tables["broad_key"]        
            df = broad_gene_names.join(df, how = "left")
            df = df.reset_index().rename(columns= {"transcript_id": "Transcript_ID","gene_id":"Database_ID"}).set_index(["Name","Transcript_ID","Database_ID"])
            df = df.rename(columns = broad_dict).sort_index() 
            df = df.T
            df.index.name = "Patient_ID"
            self.save_df(df_type, df)
