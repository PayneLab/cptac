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
import numpy as np
import os
import warnings

from cptac.cancers.source import Source
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, ReindexMapError

class AwgCoad(Source):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the colon dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        self.valid_versions = ["0.0", "0.0.1"] # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.

        self.data_files = {
            "0.0": {
                "annotation"              : "clinical.tsi.gz",
                "miRNA"                   : "miRNA.cct.gz",
                "somatic_mutation_binary" : "mutation_binary.cbt.gz",
                "somatic_mutation"        : "mutation.txt.gz",
                "phosphoproteomics"       : ["phosphoproteomics_normal.gz", "phosphoproteomics_tumor.gz"],
                "proteomics"              : ["proteomics_normal.cct.gz", "proteomics_tumor.cct.gz"],
                "transcriptomics"         : "transcriptomics.gz"},
            "0.0.1": {
                "annotation"              : "clinical.tsi.gz",
                "followup"                : "Colon_One_Year_Clinical_Data_20160927.xls",
                "CNV"                     : "Human__CPTAC_COAD__VU__SCNA__ExomeSeq__01_28_2016__BCM__Gene__BCM_CopyWriteR_GISTIC2.cct.gz",
                "miRNA"                   : "miRNA.cct.gz",
                "somatic_mutation_binary" : "mutation_binary.cbt.gz",
                "somatic_mutation"        : "mutation.txt.gz",
                "phosphoproteomics"       : ["phosphoproteomics_normal.gz", "phosphoproteomics_tumor.gz"],
                "proteomics"              : ["proteomics_normal.cct.gz", "proteomics_tumor.cct.gz"],
                "transcriptomics"         : "transcriptomics.gz"},
        }

        self.load_functions = {
            "clinical"                  : self.load_annotation,
            "CNV"                       : self.load_CNV,
            "derived_molecular"         : self.load_annotation,
            "followup"                  : self.load_followup,
            "miRNA"                     : self.load_miRNA,
            "phosphoproteomics"         : self.load_phosphoproteomics,
            "proteomics"                : self.load_proteomics,
            "somatic_mutation"          : self.load_somatic_mutation,
            "somatic_mutation_binary"   : self.load_somatic_mutation_binary,
            "transcriptomics"           : self.load_transcriptomics,
        }

        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        super().__init__(cancer_type="coad", source='awg', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

    # Overload the default how_to_cite function, to provide the specific publication information for the Colon dataset
    def how_to_cite(self):
        """Print instructions for citing the data."""
        super().how_to_cite(cancer_type='colorectal cancer', pmid=31031003)

    def load_annotation(self):
        if 'clinical' not in self._data or 'derived_molecular' not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files('annotation')

            df = pd.read_csv(file_path, sep='\t',index_col=0)
            df = df.sort_index()
            df = df.transpose()

            # Separate clinical and derived molecular dataframes
            clinical = df.drop(columns=['StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'TumorPurity','immuneSubtype', 'CIN', 'Integrated.Phenotype', 'Transcriptomic_subtype', 'Proteomic_subtype', 'mutation_rate', 'Mutation_Phenotype'])
            derived_molecular = df[['StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'TumorPurity', 'immuneSubtype', 'CIN', 'Integrated.Phenotype', 'Transcriptomic_subtype', 'Proteomic_subtype', 'mutation_rate', 'Mutation_Phenotype']]

            # Format the dataframes
            clinical = clinical.apply(pd.to_numeric, errors="ignore")
            derived_molecular = derived_molecular.apply(pd.to_numeric, errors="ignore")
            derived_molecular = derived_molecular.sort_index(axis="columns")

            # save dataframes into self._data
            self.save_df('clinical', clinical)
            self.save_df('derived_molecular', derived_molecular)

    def load_CNV(self):
        df_type = 'CNV'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t',index_col=0)
            df = df.sort_index()
            df = df.transpose()
            # save df in self._data
            self.save_df(df_type, df)

    def load_followup(self):
        df_type = 'followup'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # parse followup file
            df = pd.read_excel(file_path)
            # Replace redundant values for "not reported" with NaN
            nan_equivalents = ['Not Reported/ Unknown', 'Reported/ Unknown', 'Not Applicable',
                'na', 'unknown', 'Not Performed', 'Unknown tumor status']

            df = df.replace(nan_equivalents, np.nan)

            # Rename and set index
            df = df.rename(columns={'PPID': 'Patient_ID'})
            df = df.set_index("Patient_ID")
            df = df.sort_index()

            # save df in self._data
            self.save_df(df_type, df)

    def load_miRNA(self):
        df_type = 'miRNA'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t',index_col=0)
            df = df.sort_index()
            df = df.transpose()
            # save df in self._data
            self.save_df(df_type, df)

    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_paths = self.locate_files(df_type)

            phosphoproteomics_dfs = {}
            for file_path in file_paths:
                file_name = file_path.split(os.sep)[-1]
                df = pd.read_csv(file_path, sep='\t',index_col=0)
                df = df.sort_index()
                df = df.transpose()

                if file_name == "phosphoproteomics_normal.gz":
                    phosphoproteomics_dfs["normal"] = df
                else:
                    phosphoproteomics_dfs["tumor"] = df

            # Get phosphoproteomics dataframes, so we can process and combine them
            phos_tumor = phosphoproteomics_dfs["tumor"]
            phos_normal = phosphoproteomics_dfs["normal"]

            # Mark entries in phosphoproteomics_normal dataframe with .N at the end of the ID
            phos_normal = phos_normal.set_index(phos_normal.index + '.N')

            # Combine the two phosphoproteomics dataframes into one dataframe
            phos_combined = pd.concat([phos_tumor, phos_normal])

            # Create our phosphoproteomics columns multiindex
            multiindex = phos_combined.columns.str.split('[_:]', expand=True) # Split the column names into their constituent parts
            multiindex = multiindex.droplevel([2, 4]) # The third level is just empty strings, and the fifth is a duplicate of the second
            multiindex = multiindex.set_names(["Name", "Site", "Database_ID"])
            phos_combined.columns = multiindex
            phos_combined = phos_combined.sort_index(axis=1) # Put all the columns in alphabetical order

            # save dataframe in self._data
            self.save_df(df_type, phos_combined)

    def load_proteomics(self):
        df_type = 'proteomics'
        if df_type not in self._data:

            file_path_list = self.locate_files(df_type)

            for file_path in file_path_list:
                path_elements = file_path.split(os.sep) # Get a list of the levels of the path
                file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

                if file_name == "proteomics_normal.cct.gz":
                    df_normal = pd.read_csv(file_path, sep='\t', index_col=0)
                    df_normal = df_normal.sort_index()
                    df_normal = df_normal.transpose()
                    # append .N to patient ids from the normal table
                    df_normal = df_normal.rename(index=lambda s: s + ".N")

                if file_name == "proteomics_tumor.cct.gz":
                    df_tumor = pd.read_csv(file_path, sep='\t', index_col=0)
                    df_tumor = df_tumor.sort_index()
                    df_tumor = df_tumor.transpose()

            # merge tumor and normal data
            df_combined = pd.concat([df_normal, df_tumor])
            df_combined.index.name = "Patient_ID"
            df_combined.columns.name = "Name"

            # save dataframe in self._data
            self.save_df(df_type, df_combined)

    def load_somatic_mutation(self):
        df_type = 'somatic_mutation'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t',index_col=0)
            df = df.sort_index()
            df = df.sort_values(by="SampleID")
            df = df.reset_index()
            df = df[["SampleID","Gene","Variant_Type","Protein_Change"]]
            df = df.drop_duplicates(keep="first") # Get rid of rows that are now duplicates since we didn't keep the mRNA column. We do this before setting the index, because drop_duplicates doesn't consider the index.
            df = df.rename(columns={"SampleID":"Patient_ID", "Variant_Type":"Mutation", "Protein_Change":"Location"})
            df = df.sort_values(by=["Patient_ID", "Gene"])
            df = df.set_index("Patient_ID") # We only do this after the drop_duplicates call above because drop_duplicates doesn't consider the index, but we of course want the Patient_ID to be considered when identifying duplicate rows to drop.

            # save df in self._data
            self.save_df(df_type, df)

    def load_somatic_mutation_binary(self):
        df_type = 'somatic_mutation_binary'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t',index_col=0)
            df = df.sort_index()
            df = df.transpose()

            # save df in self._data
            self.save_df(df_type, df)

    def load_transcriptomics(self):
        df_type = 'transcriptomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t',index_col=0)
            df = df.sort_index()
            df = df.transpose()
            # save df in self._data
            self.save_df(df_type, df)