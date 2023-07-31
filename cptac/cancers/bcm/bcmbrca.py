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

# Importing necessary libararies
import pandas as pd
from cptac.cancers.source import Source

class BcmBrca(Source):
    def __init__(self, no_internet=False):
        """
        Define which bcmbrca dataframes as are available in the self.load_functions dictionary variable, with names as keys.

        Parameters:
        no_internet (bool, optional): If set to True, skips the index update step which requires an internet connection.
        Default is False.
        """

        # Define the needed variables and pass them to the parent Dataset __init__ function
        self.data_files = {
            "transcriptomics" : "BRCA-gene_RSEM_tumor_normal_UQ_log2(x+1)_BCM.txt.gz", 
            "mapping" : "gencode.v34.basic.annotation-mapping.txt.gz",
            "proteomics" : "BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt.gz",
            "phosphoproteomics" : "BRCA_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt"
        }
        
        self.load_functions = {
            'transcriptomics' : self.load_transcriptomics,
            'proteomics' : self.load_proteomics,
            'phosphoproteomics' : self.load_phosphoproteomics
        }
        
        # Call the parent class __init__ function
        super().__init__(cancer_type="brca", source='bcm', data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)
        
    def load_mapping(self):
        """
        Load and parse all files for mapping.
        """
        df_type = 'mapping'

        # Check if helper tables have already been loaded
        if not self._helper_tables:
            file_path = self.locate_files(df_type)
            
            # Load the file, select needed columns, set index, remove duplicates
            df = pd.read_csv(file_path, sep='\t')
            df = df[["gene","gene_name"]]
            df = df.set_index("gene")
            df = df.drop_duplicates()

            # Save mapping in helper tables
            self._helper_tables["gene_key"] = df 
            
    def load_transcriptomics(self):
        """
        Load and parse all files for bcm brca transcriptomics data
        """
        df_type = 'transcriptomics'

        # Check if data is already loaded
        if df_type not in self._data:
            # Get file path to the correct data
            file_path = self.locate_files(df_type)

            # Load and process the file
            df = pd.read_csv(file_path, sep='\t')
            df.index.name = 'gene'

            # Load mapping information
            self.load_mapping()
            gene_key = self._helper_tables["gene_key"]

            # Join gene_key to df, reset index, rename columns, set new index and sort
            transcript = gene_key.join(df,how = "inner") #keep only gene_ids with gene names
            transcript = transcript.reset_index()
            transcript = transcript.rename(columns={"gene_name":"Name","gene":"Database_ID"})
            transcript = transcript.set_index(["Name", "Database_ID"])
            transcript = transcript.sort_index() #alphabetize
            transcript = transcript.T
            transcript.index.name = "Patient_ID"

            df = transcript

            # Save df in data
            self.save_df(df_type, df)
            
    
    def load_proteomics(self):
        """
        Load and parse all files for bcm brca proteomics data
        """
        df_type = 'proteomics'

        # Check if data is already loaded
        if df_type not in self._data:
            # Get file path to the correct data
            file_path = self.locate_files(df_type)

            # Load and process the file
            df = pd.read_csv(file_path, sep='\t')
            df.index.name = 'gene'

            df.set_index('idx', inplace=True)
            # Load mapping information
            self.load_mapping()
            gene_key = self._helper_tables["gene_key"]

            # Join gene_key to df, reset index, rename columns, set new index and sort
            proteomics = gene_key.join(df, how='inner')
            proteomics = gene_key.join(df, how='inner')
            proteomics = proteomics.reset_index()
            proteomics = proteomics.rename(columns={"index": "Database_ID", "gene_name": "Name"})
            proteomics = proteomics.set_index(["Name", "Database_ID"])
            proteomics = proteomics.sort_index()  # alphabetize
            proteomics = proteomics.T
            proteomics.index.name = "Patient_ID"

            df = proteomics

            # Save df in data
            self.save_df(df_type, df)

    def load_phosphoproteomics(self):
        """
        Load and parse all files for bcm brca phosphoproteomics data
        """
        df_type = 'phosphoproteomics'

        # Check if data is already loaded
        if df_type not in self._data:
            # Get file path to the correct data
            file_path = self.locate_files(df_type)

            # Load and process the file
            df = pd.read_csv(file_path, sep='\t')
            df.index.name = 'gene'

            df.set_index('idx', inplace=True)
            # Load mapping information
            self.load_mapping()
            gene_key = self._helper_tables["gene_key"]

            # Join gene_key to df, reset index, rename columns, set new index and sort
            phosphoproteomics = gene_key.join(df, how='inner')
            phosphoproteomics = gene_key.join(df, how='inner')
            phosphoproteomics = phosphoproteomics.reset_index()
            phosphoproteomics = phosphoproteomics.rename(columns={"index": "Database_ID", "gene_name": "Name"})
            phosphoproteomics = phosphoproteomics.set_index(["Name", "Database_ID"])
            phosphoproteomics = phosphoproteomics.sort_index()  # alphabetize
            phosphoproteomics = phosphoproteomics.T
            phosphoproteomics.index.name = "Patient_ID"

            df = phosphoproteomics

            # Save df in data
            self.save_df(df_type, df)