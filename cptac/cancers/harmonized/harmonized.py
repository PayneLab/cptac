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
from cptac.cancers.source import Source

class Harmonized(Source):
    def __init__(self, filter_type, no_internet=False):
        """Define which dataframes as are available in the self.load_functions dictionary variable, with names as keys.

        Parameters:
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        filter_type (str): The cancer type for which you want information. Harmonized keeps all data in a single table, so to get data on a single cancer type all other types are filtered out.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        self.data_files = {
            "somatic_mutation" : "Union_Maf_Broad_WashU.maf.gz",
            "ancestry_prediction" : "washu_mssm_consensus_ancestries.tsv.gz"
        }

        self.load_functions = {
            'somatic_mutation' : self.load_somatic_mutation,
            'ancestry_prediction': self.load_ancestry_prediction
        }

        # Call the parent class __init__ function, cancer_type is dynamic and based on whatever cancer is being filtered for
        super().__init__(cancer_type=filter_type, source='harmonized', data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

    def load_somatic_mutation(self):
        df_type = 'somatic_mutation'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # which cancer_type goes with which cancer in the harmonized table
            tumor_codes = {'brca':'BRCA', 'ccrcc':'CCRCC',
                           'ucec':'UCEC', 'gbm':'GBM', 'hnscc':'HNSCC',
                           'lscc':'LSCC', 'luad':'LUAD', 'pdac':'PDA',
                           'hcc':'HCC', 'coad':'CO', 'ov':'OV'}

            df = pd.read_csv(file_path, sep='\t', low_memory = False)
            df = df.loc[df['COHORT'] == tumor_codes[self.cancer_type]]
            df['Patient_ID'] = df.loc[:, 'Tumor_Sample_Barcode']
            df = df.rename(columns={
                     "Hugo_Symbol":"Gene",
                     "Variant_Classification":"Mutation",
                     "Protein_Change":"Location"})

            df = df.set_index("Patient_ID")
            df = df[ ['Gene'] + ["Mutation"] + ["Location"] + [ col for col in df.columns if col not in ["Gene","Mutation","Location"] ] ]
            df.index = df.index.str.replace(r"_T", "", regex=True) # data based on Tumor and Normal. Remove _T

            # save df in self._data
            self.save_df(df_type, df)


        # I'm not sure what this comment is, but perhaps the harmonized parsing is not finished and it will come in handy
        # Since the cancer type is the same as filter type this would now be self.cancer_type and would not include pancan
        '''
        if filter_type == 'pancanucec':  
            print("True")
            mut_df = self._data["somatic_mutation"]
            mut_df = mut_df.loc[mut_df.index[~ mut_df.index.str.contains('NX', regex = True)]] # Drop quality control 
            self._data["somatic_mutation"] = mut_df
        '''

    def load_ancestry_prediction(self):
        df_type = 'ancestry_prediction'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # which cancer_type goes with which cancer in the mssm table
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