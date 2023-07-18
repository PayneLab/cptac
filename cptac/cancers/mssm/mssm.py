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

class Mssm(Source):
    def __init__(self, filter_type, no_internet=False):
        """Define which dataframes as are available in the self.load_functions dictionary variable, with names as keys.

        Parameters:
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        filter_type (str): The cancer type for which you want information. Mssm keeps all data in a single table, so to get data on a single cancer type all other types are filtered out.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        self.data_files = {
            "clinical" : "clinical_Pan-cancer.May2022.tsv.gz",
        }

        self.load_functions = {
            'clinical' : self.load_clinical,
            'medical_history': self.load_medical_history,
            'follow-up': self.load_followup
        }

        # Mssm is special in that it keeps information for all cancer types in a single table
        # We can still make this look like the other sources, but it works a little different under the hood
        # Essentially Mssm always loads the entire clinical table, then filters it based on filter_type

        # Call the parent class __init__ function, cancer_type is dynamic and based on whatever cancer is being filtered for
        super().__init__(cancer_type=filter_type, source='mssm', data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

    # This is all clinical data, but mssm will have other load functions to get only a subset of this data
    def load_clinical(self):
        df_type = 'clinical'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # which cancer_type goes with which cancer in the mssm table
            tumor_codes = {'brca':'BR', 'ccrcc':'CCRCC',
                           'ucec':'UCEC', 'gbm':'GBM', 'hnscc':'HNSCC',
                           'lscc':'LSCC', 'luad':'LUAD', 'pdac':'PDA',
                           'hcc':'HCC', 'coad':'CO', 'ov':'OV'}

            df = pd.read_csv(file_path, sep='\t')
            df = df.loc[df['tumor_code'] == tumor_codes[self.cancer_type]]
            df = df.loc[df['discovery_study'] != 'No'] # Only keep discovery study = 'Yes' or 'na'
            df = df.set_index("case_id")
            df.index.name = 'Patient_ID'
            df = df.sort_values(by=["Patient_ID"])

            self.save_df(df_type, df)

        return self._data[df_type]

    def load_medical_history(self):
        df_type = 'medical_history'

        if df_type not in self._data:
            # First load the clinical dataframe as you would in load_clinical
            clinical_df = self.load_clinical()

            # Then filter only the medical history columns
            medical_history_df = clinical_df[[col for col in clinical_df.columns if 'medical_history' in col]]

            # Finally, save the medical history dataframe in the _data dictionary
            self.save_df(df_type, medical_history_df)

        return self._data[df_type]

    def load_followup(self):
        df_type = 'follow-up'

        if df_type not in self._data:
            # First load the clinical dataframe as you would in load_clinical
            clinical_df = self.load_clinical()

            # Then filter only the follow-up columns
            followup_df = clinical_df[[col for col in clinical_df.columns if 'follow-up' in col]]

            # Finally, save the follow-up dataframe in the _data dictionary
            self.save_df(df_type, followup_df)

        return self._data[df_type]


