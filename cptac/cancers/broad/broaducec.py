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
from pyranges import read_gtf
from cptac.cancers.source import Source

class BroadUcec(Source):
    """
    This class handles the loading and providing of data for the Uterine Corpus Endometrial Carcinoma (UCEC) from the Broad Institute.
    It inherits from the Source parent class.
    """

    def __init__(self, no_internet=False):
        """
        Initializes the BroadUcec object. It specifies the available dataframes and their respective load functions.
        
        Parameters:
        no_internet (bool): If True, the index update step (which requires internet connection) will be skipped. Defaults to False.
        """

        self.data_files = {
            "transcriptomics" : "UCEC.rsem_transcripts_tpm.txt.gz",
            "mapping" : ["sample_descriptions.tsv.gz", "gencode.v34.GRCh38.genes.collapsed_only.gtf.gz", "aliquot_to_patient_ID.tsv.gz"]
        }
        
        self.load_functions = {
            'transcriptomics' : self.load_transcriptomics,
        }
        
        super().__init__(cancer_type="ucec", source='broad', data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

    def load_mapping(self):
        """
        This method populates the _helper_tables dictionary with dataframes necessary for mapping and identifying data.
        It loads three types of mapping data - broad_keys, broad_gene_names, and map_ids.
        """
        if not self._helper_tables:
            self._load_broad_keys()
            self._load_broad_gene_names()
            self._load_map_ids()

    def load_transcriptomics(self):
        """
        This method loads the transcriptomics data and combines it with the mapping data to create a dataframe with appropriate gene and patient identifiers.
        """
        df_type = 'transcriptomics'

        if df_type not in self._data:
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep="\t").set_index(["transcript_id","gene_id"])
            
            # Add gene names to transcriptomic data
            self.load_mapping()
            df = self._add_gene_names(df)
            df = self._rename_with_identifiers(df)
            self.save_df(df_type, df)

    # Additional methods to split the loading process
    def _load_broad_keys(self):
        """
        Helper method to load the broad keys dataframe.
        """
        file_path = self.locate_files("sample_descriptions.tsv.gz")
        broad_key = pd.read_csv(file_path, sep="\t")
        broad_key = broad_key.loc[broad_key['cohort'] == "UCEC"][["sample_id","GDC_id","tissue_type"]].set_index("sample_id")
        broad_key['GDC_id'] = broad_key['GDC_id'].str[:9]
        broad_key["Patient_ID"] = broad_key["GDC_id"] + broad_key["tissue_type"]
        broad_key.Patient_ID = broad_key.Patient_ID.str.replace(r"Tumor", "", regex=True)
        broad_key.Patient_ID = broad_key.Patient_ID.str.replace(r"Normal", ".N", regex=True)
        self._helper_tables["broad_key"] = broad_key.to_dict()["Patient_ID"]

    def _load_broad_gene_names(self):
        """
        Helper method to load the broad gene names dataframe.
        """
        file_path = self.locate_files("gencode.v34.GRCh38.genes.collapsed_only.gtf.gz")
        broad_gene_names = read_gtf(file_path).as_df()
        broad_gene_names = broad_gene_names[["gene_name","gene_id"]].rename(columns= {"gene_name":"Name"}).set_index("gene_id").drop_duplicates()
        self._helper_tables["broad_gene_names"] = broad_gene_names

    def _load_map_ids(self):
        """
        Helper method to load the map ids dataframe.
        """
        file_path = self.locate_files("aliquot_to_patient_ID.tsv.gz")
        self._helper_tables["map_ids"] = pd.read_csv(file_path, sep = "\t", index_col = 0)

    def _add_gene_names(self, df):
        """
        Helper method to join the dataframe with the broad gene names.
        """
        df = self._helper_tables["broad_gene_names"].join(df, how = "left").reset_index()
        df = df.rename(columns= {"transcript_id": "Transcript_ID","gene_id":"Database_ID"}).set_index(["Name","Transcript_ID","Database_ID"])
        return df

    def _rename_with_identifiers(self, df):
        """
        Helper method to rename the dataframe with the identifiers.
        """
        broad_dict = self._helper_tables["broad_key"]
        aliquot_dict = self._helper_tables["map_ids"].to_dict()["patient_ID"]
        df = df.rename(columns = broad_dict).rename(columns = aliquot_dict).sort_index()
        df = df.T.groupby("index", level = 0).mean()  # Average duplicates 
        df.index.name = "Patient_ID"
        return df
