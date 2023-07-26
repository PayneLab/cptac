#   Copyright 2023 User Name
#   Licensed under the Apache License, Version 2.0 (the "License");
#   You may not use this file except in compliance with the License.
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
        """Initialize Mssm object and load the available dataframes.

        Args:
        filter_type (str): The cancer type for which you want information. Mssm keeps all data in a single table, so to get data on a single cancer type all other types are filtered out.
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. Default is False.
        """
        self.data_files = {
            "clinical" : "clinical_Pan-cancer.May2022.tsv.gz",
        }

        self.load_functions = {
            'clinical' : self.load_clinical,
            'medical_history': self.load_medical_history,
            'follow-up': self.load_followup
        }

        super().__init__(cancer_type=filter_type, source='mssm', data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

    def load_clinical(self):
        """Load the clinical data and filters it based on filter_type.
        
        Returns:
        DataFrame: The filtered clinical data.
        """
        df_type = 'clinical'
        if df_type not in self._data:
            file_path = self.locate_files(df_type)
            tumor_codes = {'brca':'BR', 'ccrcc':'CCRCC',
                           'ucec':'UCEC', 'gbm':'GBM', 'hnscc':'HNSCC',
                           'lscc':'LSCC', 'luad':'LUAD', 'pdac':'PDA',
                           'hcc':'HCC', 'coad':'CO', 'ov':'OV'}

            df = pd.read_csv(file_path, sep='\t')
            df = df[df['tumor_code'] == tumor_codes[self.cancer_type]]
            df = df[df['discovery_study'] != 'No'] # Only keep discovery study = 'Yes' or 'na'
            df.set_index("case_id", inplace=True)
            df.index.name = 'Patient_ID'
            df.sort_values(by=["Patient_ID"], inplace=True)

            self.save_df(df_type, df)

        return self._data[df_type]

    def load_medical_history(self):
        """Load the medical history data by filtering the clinical data.
        
        Returns:
        DataFrame: The filtered medical history data.
        """
        df_type = 'medical_history'
        if df_type not in self._data:
            clinical_df = self.load_clinical()
            medical_history_df = clinical_df[[col for col in clinical_df.columns if 'medical_history' in col]]
            self.save_df(df_type, medical_history_df)

        return self._data[df_type]

    def load_followup(self):
        """Load the follow-up data by filtering the clinical data.
        
        Returns:
        DataFrame: The filtered follow-up data.
        """
        df_type = 'follow-up'
        if df_type not in self._data:
            clinical_df = self.load_clinical()
            followup_df = clinical_df[[col for col in clinical_df.columns if 'follow-up' in col]]
            self.save_df(df_type, followup_df)

        return self._data[df_type]