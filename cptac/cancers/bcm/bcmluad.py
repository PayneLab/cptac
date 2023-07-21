# Copyright 2018 Samuel Payne sam_payne@byu.edu
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pandas as pd
from cptac.cancers.source import Source

class BcmLuad(Source):
    """Defines the BcmLuad class, which handles the loading of BCM LUAD data sets.

    Attributes:
        data_files (dict): File names for different types of data.
        load_functions (dict): Functions to load each type of data.
    """
    
    def __init__(self, no_internet=False):
        """Initialize BcmLuad object.

        Parameters:
            no_internet (bool, optional): If True, skip the index update step that requires an internet connection.
                                         Default is False.
        """
        self.data_files = {
            "circular_RNA" : "LUAD-circRNA_rsem_tumor_normal_UQ_log2(x+1)_BCM.txt.gz",
            "mapping" : "gencode.v34.basic.annotation-mapping.txt.gz",
            "transcriptomics" : "LUAD-gene_rsem_removed_circRNA_tumor_normal_UQ_log2(x+1)_BCM.txt.gz"
        }
        
        self.load_functions = {
            'circular_RNA' : self.load_circular_RNA,
            'transcriptomics' : self.load_transcriptomics,
        }

        super().__init__(cancer_type="luad", source='bcm', data_files=self.data_files, 
                         load_functions=self.load_functions, no_internet=no_internet)

    def load_circular_RNA(self):
        """Load the circular RNA data."""
        df_type = 'circular_RNA'

        if df_type not in self._data:
            file_path = self.locate_files(df_type)
            df = self.prepare_circular_RNA_df(file_path)
            self.save_df(df_type, df)

    def load_mapping(self):
        """Load the mapping data."""
        df_type = 'mapping'

        if not self._helper_tables:
            file_path = self.locate_files(df_type)
            df = self.prepare_mapping_df(file_path)
            self._helper_tables["gene_key"] = df

    def load_transcriptomics(self):
        """Load the transcriptomics data."""
        df_type = 'transcriptomics'

        if df_type not in self._data:
            file_path = self.locate_files(df_type)
            df = self.prepare_transcriptomics_df(file_path)
            self.save_df(df_type, df)

    def prepare_circular_RNA_df(self, file_path):
        """Prepare the circular RNA data frame.

        Parameters:
            file_path (str): Path to the circular RNA data file.

        Returns:
            df (DataFrame): Prepared circular RNA data frame.
        """
        df = pd.read_csv(file_path, sep='\t')
        df = df.rename_axis('INDEX').reset_index()
        df[["circ","chrom","start","end","gene"]] = df.INDEX.str.split('_', expand=True)
        df["circ_chromosome"] = df["circ"] +"_" + df["chrom"]
        df = df.set_index('gene')

        # Add gene names to circular RNA data
        self.load_mapping()
        gene_key = self._helper_tables["gene_key"]
        df = gene_key.join(df, how = "inner")
        df = df.reset_index()
        df = df.rename(columns= {"gene_name": "Name", "gene": "Database_ID"}) 
        df = df.set_index(["Name","circ_chromosome", "start","end","Database_ID"])
        df.drop(['INDEX', 'circ', 'chrom'], axis=1, inplace=True) 
        df = df.sort_index()
        df = df.T
        df.index = df.index.str.replace(r"_T", "", regex=True)
        df.index = df.index.str.replace(r"_A", ".N", regex=True)
        df.index.name = "Patient_ID"

        return df

    def prepare_mapping_df(self, file_path):
        """Prepare the mapping data frame.

        Parameters:
            file_path (str): Path to the mapping data file.

        Returns:
            df (DataFrame): Prepared mapping data frame.
        """
        df = pd.read_csv(file_path, sep='\t')
        df = df[["gene","gene_name"]] 
        df = df.set_index("gene")
        df = df.drop_duplicates()

        return df

    def prepare_transcriptomics_df(self, file_path):
        """Prepare the transcriptomics data frame.

        Parameters:
            file_path (str): Path to the transcriptomics data file.

        Returns:
            df (DataFrame): Prepared transcriptomics data frame.
        """
        df = pd.read_csv(file_path, sep='\t')
        df.index.name = 'gene'
        
        # Add gene names to transcriptomic data
        self.load_mapping()
        gene_key = self._helper_tables["gene_key"]
        transcript = gene_key.join(df, how="inner") 
        transcript = transcript.reset_index()
        transcript = transcript.rename(columns={"gene_name":"Name","gene":"Database_ID"})
        transcript = transcript.set_index(["Name", "Database_ID"])
        transcript = transcript.sort_index() 
        transcript = transcript.T
        transcript.index = transcript.index.str.replace(r"_T", "", regex=True)
        transcript.index = transcript.index.str.replace(r"_A", ".N", regex=True)
        transcript.index.name = "Patient_ID"

        return transcript
