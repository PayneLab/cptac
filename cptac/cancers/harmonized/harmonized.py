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
from cptac.cancers.source import Source

class Harmonized(Source):
    TUMOR_CODES = {
        'brca': 'BRCA', 'ccrcc': 'CCRCC', 'ucec': 'UCEC',
        'gbm': 'GBM', 'hnscc': 'HNSCC', 'lscc': 'LSCC',
        'luad': 'LUAD', 'pdac': 'PDA', 'hcc': 'HCC',
        'coad': 'CO', 'ov': 'OV'
    }

    def __init__(self, filter_type, no_internet=False):
        """
        Constructor method that initializes the Harmonized object.

        Parameters:
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. 
        Default is False.
        filter_type (str): The cancer type for which you want information. Harmonized keeps all data in a single table, 
        so to get data on a single cancer type all other types are filtered out.
        """

        self.data_files = {
            "somatic_mutation" : "PanCan_Union_Maf_Broad_WashU_v1.1.maf.gz",
            "ancestry_prediction" : "washu_mssm_consensus_ancestries.tsv.gz"
        }

        self.load_functions = {
            'somatic_mutation' : self.load_somatic_mutation,
            'ancestry_prediction': self.load_ancestry_prediction
        }

        super().__init__(cancer_type=filter_type, source='harmonized', data_files=self.data_files, 
                         load_functions=self.load_functions, no_internet=no_internet)

    def load_somatic_mutation(self):
        """
        Load the somatic mutation data for the specific cancer type.

        Returns:
        None. The data is stored internally in the _data dictionary.
        """
        df_type = 'somatic_mutation'

        if df_type not in self._data:
            # Locate the appropriate file
            file_path = self.locate_files(df_type)

            tumor_codes = {'brca':'BRCA', 'ccrcc':'CCRCC',
                           'ucec':'UCEC', 'gbm':'GBM', 'hnscc':'HNSCC',
                           'lscc':'LSCC', 'luad':'LUAD', 'pdac':'PDA',
                           'hcc':'HCC', 'coad':'CO', 'ov':'OV'}

            # Load the data, filtering for the specific tumor type
            df = pd.read_csv(file_path, sep='\t', low_memory = False)
            df = df.loc[df['COHORT'] == self.TUMOR_CODES[self.cancer_type]]

            # Add patient ID from the sample barcode
            df['Patient_ID'] = df.loc[:, 'Tumor_Sample_Barcode']

            # Rename certain columns for clarity
            df = df.rename(columns={
                     "Hugo_Symbol":"Gene",
                     "Variant_Classification":"Mutation",
                     "Protein_Change":"Location"})

            # Set patient ID as index and order columns
            df = df.set_index("Patient_ID")
            cols_order = ['Gene', 'Mutation', 'Location'] + [col for col in df.columns if col not in {'Gene', 'Mutation', 'Location'}]
            df = df[cols_order]
            
            # Remove trailing "_T" from patient IDs
            df.index = df.index.str.replace(r"_T", "", regex=True)

            # Save the dataframe internally
            self.save_df(df_type, df)

    def load_ancestry_prediction(self):
        """
        Load the ancestry prediction data for the specific cancer type.
        """
        df_type = 'ancestry_prediction'

        if df_type not in self._data:
            file_path = self.locate_files(df_type)

            tumor_codes = {'brca':'BR', 'ccrcc':'CCRCC',
                           'ucec':'UCEC', 'gbm':'GBM', 'hnscc':'HNSCC',
                           'lscc':'LSCC', 'luad':'LUAD', 'pdac':'PDA',
                           'hcc':'HCC', 'coad':'CO', 'ov':'OV'}

            df = pd.read_csv(file_path, sep='\t')
            df = df.loc[df['cancer_type'] == tumor_codes[self.cancer_type]]
            df = df.set_index('case_id')
            df.index.name = 'Patient_ID'
            df = df.sort_values(by=["Patient_ID"])
            df = df.drop(columns='cptac_cohort')  # drop unnecessary column
            self.save_df(df_type, df)

        return self._data[df_type]