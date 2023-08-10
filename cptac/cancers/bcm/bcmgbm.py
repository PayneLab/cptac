#   Copyright 2023 Samuel Payne sam_payne@byu.edu
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# Importing necessary libraries
import pandas as pd
from cptac.cancers.source import Source

class BcmGbm(Source):
    def __init__(self, no_internet=False):
        """
        Define bcmgbm dataframes available in self.load_functions with names as keys.
        
        Parameters:
        no_internet (bool, optional): If set to True, skips the index update step which requires an internet connection.
        Default is False.
        """
        
        # Define the needed variables and pass them to the parent Dataset class __init__ function
        self.data_files = {
            "transcriptomics" : "GBM-gene_rsem_removed_circRNA_tumor_normal_UQ_log2(x+1)_BCM.txt.gz",
            "mapping" : "gencode.v34.basic.annotation-mapping.txt.gz",
            "circular_RNA" : "GBM-circRNA_rsem_tumor_normal_UQ_log2(x+1)_BCM.txt.gz",
            "proteomics": "GBM_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt.gz",
            "phosphoproteomics" : "GBM_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt",
            "CNV" : "GBM_WES_CNV_gene_ratio_log2.txt.gz",
            "miRNA" : "GBM_miRNAseq_mature_miRNA_RPM_log2_Tumor.txt.gz"
        }

        self.load_functions = {
            'circular_RNA' : self.load_circular_RNA,
            'transcriptomics' : self.load_transcriptomics,
            'proteomics' : self.load_proteomics,
            'phosphoproteomics' : self.load_phosphoproteomics,
            'CNV' : self.load_CNV,
            'miRNA' : self.load_miRNA
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="gbm", source='bcm', data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)


    def load_circular_RNA(self):
        """
        Load and parse all files for bcm gbm circular RNA data.
        """
        df_type = 'circular_RNA'
        
        # Check if data is already loaded
        if df_type not in self._data:
            # Get file path to the correct data
            file_path = self.locate_files(df_type)

            # Load and parse the file 
            df = pd.read_csv(file_path, sep="\t")
            df = df.rename_axis('INDEX').reset_index()
            df[["circ","chrom","start","end","gene"]] = df.INDEX.str.split('_', expand=True)
            df["circ_chromosome"] = df["circ"] +"_" + df["chrom"]
            df = df.set_index('gene')

            # Load mapping information and add it to circular RNA data
            self.load_mapping()
            gene_key = self._helper_tables["gene_key"]
            df = gene_key.join(df, how = "inner")
            df = df.reset_index()
            df = df.rename(columns= {"gene_name": "Name","gene":"Database_ID"}) # change names to match cptac package
            df = df.set_index(["Name","circ_chromosome", "start","end","Database_ID"]) #create multi-index
            df.drop(['INDEX', 'circ', 'chrom'], axis=1, inplace=True) 
            df = df.sort_index()
            df = df.T
            df.index = df.index.str.replace(r"_T", "", regex=True) # remove Tumor label. All samples are tumor samples
            df.index.name = "Patient_ID"

            # Save df in data
            self.save_df(df_type, df)


    def load_mapping(self):
        """
        Load and parse all files for mapping. These will be used for circular RNA and transcriptomics loading.
        """
        df_type = 'mapping'
        
        # Check if helper tables have already been loaded
        if not self._helper_tables:
            file_path = self.locate_files(df_type)

            # Load the file, select needed columns, set index, remove duplicates
            df = pd.read_csv(file_path, sep="\t")
            df = df[["gene","gene_name"]] 
            df = df.set_index("gene")
            df = df.drop_duplicates()

            # Save mapping in helper tables
            self._helper_tables["gene_key"] = df 


    def load_transcriptomics(self):
        """
        Load and parse all files for bcm gbm transcriptomics data.
        """
        df_type = 'transcriptomics'
        
        # Check if data is already loaded
        if df_type not in self._data:
            # Get file path to the correct data
            file_path = self.locate_files(df_type)

            # Load the file 
            df = pd.read_csv(file_path, sep="\t")
            df.index.name = 'gene'

            # Load mapping information and add it to transcriptomics data
            self.load_mapping()
            gene_key = self._helper_tables["gene_key"]
            transcript = gene_key.join(df, how="inner") 
            transcript = transcript.reset_index()
            transcript = transcript.rename(columns={"gene_name":"Name","gene":"Database_ID"})
            transcript = transcript.set_index(["Name", "Database_ID"])
            transcript = transcript.sort_index() 
            transcript = transcript.T
            transcript.index = transcript.index.str.replace(r"_T", "", regex=True)
            transcript.index.name = "Patient_ID" 

            # Save transcriptomics in data
            self.save_df(df_type, transcript)

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

            # Extract Database_ID, gene name, site, and peptide from 'idx' column
            df[['Database_ID', 'Gene_Key', 'Site', 'Peptide']] = df['idx'].str.split('|', expand=True).iloc[:,
                                                                 [0, 1, 2, 3]]

            # Load mapping information
            self.load_mapping()
            gene_key_df = self._helper_tables["gene_key"]
            mapping = gene_key_df['gene_name'].to_dict()

            # Map gene_key to get gene name
            df['Name'] = df['Gene_Key'].map(mapping)

            # Drop the 'idx' and 'Gene_Key' columns
            df.drop(columns=['idx', 'Gene_Key'], inplace=True)

            # Set the 'Name', 'Site', 'Peptide', and 'Database_ID' columns as index in this order
            df.set_index(['Name', 'Site', 'Peptide', 'Database_ID'], inplace=True)

            # Transpose the dataframe so that the patient IDs are the index
            df = df.transpose()

            # Rename the index to 'Patient_ID'
            df.index.name = 'Patient_ID'

            # Save df in data
            self.save_df(df_type, df)

    def load_CNV(self):
        """Load and process the CNV data file, and save it in the _data dictionary."""
        df_type = 'CNV'
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
            df = gene_key.join(df, how='inner')
            df = df.reset_index()
            df = df.rename(columns={"index": "Database_ID", "gene_name": "Name"})
            df = df.set_index(["Name", "Database_ID"])
            df = df.sort_index()  # alphabetize
            df = df.T
            df.index.name = "Patient_ID"

            # Save df in data
            self.save_df(df_type, df)

    def load_miRNA(self):
        """
        Load and parse all files for miRNA data
        """
        df_type = 'miRNA'

        # Check if data is already loaded
        if df_type not in self._data:
            # Get file path to the correct data
            file_path = self.locate_files(df_type)

            # Load and process the file
            df = pd.read_csv(file_path, sep='\t')

            # Here the idx will be the miRNA names and columns will be the patient IDs.
            df.set_index('idx', inplace=True)
            df.index.name = 'Name'  # Rename idx to Name
            df = df.transpose()  # Transpose the data frame to have miRNA as columns and patients as rows
            df.index.name = "Patient_ID"

            # Save df in data
            self.save_df(df_type, df)