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

class AwgOv(Source):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the ovarian dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        self.valid_versions = ["0.0", "0.0.1"] # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.

        self.data_files = {
            "0.0": {
                "clinical"          : "clinical.csv.gz",
                "CNV"               : "cnv.tsv.gz",
                "definitions"       : "definitions.txt",
                "phosphoproteomics" : "phosphoproteomics.txt.gz",
                "proteomics"        : "proteomics.txt.gz",
                "somatic_mutation"  : "somatic_38.maf.gz",
                "transcriptomics"   : "transcriptomics.tsv.gz",
                "treatment"         : "treatment.csv.gz"},
            "0.0.1": {
                "clinical"          : "clinical.csv.gz",
                "CNV"               : "cnv.tsv.gz",
                "definitions"       : "definitions.txt",
                "followup"          : "Ovary_One_Year_Clinical_Data_20160927.xls",
                "phosphoproteomics" : "phosphoproteomics.txt.gz",
                "proteomics"        : "proteomics.txt.gz",
                "somatic_mutation"  : "somatic_38.maf.gz",
                "transcriptomics"   : "transcriptomics.tsv.gz",
                "treatment"         : "treatment.csv.gz"},
        }

        self.load_functions = {
            'clinical'                : self.load_clinical,
            'CNV'                     : self.load_CNV,
            'followup'                : self.load_followup,
            'phosphoproteomics'       : self.load_phosphoproteomics,
            'proteomics'              : self.load_proteomics,
            'somatic_mutation'        : self.load_somatic_mutation,
            'transcriptomics'         : self.load_transcriptomics,
            'treatment'               : self.load_treatment,
        }

        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        super().__init__(cancer_type="ov", source='awg', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)


    def load_clinical(self):
        df_type = 'clinical'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep=',', index_col=0)
            df = df.rename(columns={"Participant_ID":"Patient_ID"})
            df = df.set_index("Patient_ID")

            # save df in self._data
            self.save_df(df_type, df)


    def load_CNV(self):
        df_type = 'CNV'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.sort_index()
            df = df.transpose()
            df = df.sort_index()

            # save df in self._data
            self.save_df(df_type, df)


    def load_somatic_mutation(self):
        df_type = 'somatic_mutation'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.reset_index()
            # The first part of the barcode is the patient id, which we need to make a Patient_ID column
            split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n=1, expand=True)
            df["Tumor_Sample_Barcode"] = split_barcode[0]
            df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]] # We only want these columns
            df = df.rename(columns={"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"})
            df = df.sort_values(by=["Patient_ID", "Gene"])
            df = df.set_index("Patient_ID")

            # save df in self._data
            self.save_df(df_type, df)


    def load_transcriptomics(self):
        df_type = 'transcriptomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.sort_index()
            df = df.transpose()
            df = df.sort_index()
            date_cols = ['1-Dec', '1-Sep', '10-Mar', '10-Sep', '11-Sep', '12-Sep',
                         '14-Sep', '15-Sep', '2-Mar', '2-Sep', '3-Mar', '3-Sep', '4-Mar',
                         '4-Sep', '5-Mar', '6-Mar', '6-Sep', '7-Mar', '7-Sep', '8-Mar',
                         '8-Sep', '9-Mar', '9-Sep']
            df = df.drop(columns=date_cols) # Drop all date values until new data is uploaded

            # save df in self._data
            self.save_df(df_type, df)


    def load_followup(self):
        df_type = 'followup'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_excel(file_path)
            # Replace redundant values for "not reported" with NaN
            nan_equivalents = ['Not Reported/ Unknown', 'Reported/ Unknown', 'Not Applicable',
                'na', 'unknown', 'Not Performed', 'Unknown tumor status', 'Unknown',
                'Unknown Tumor Status', 'Not specified']
            df = df.replace(nan_equivalents, np.nan)

            # Rename PPID to Patient_ID and set as index
            df = df.rename(columns={'PPID': 'Patient_ID'})
            df = df.set_index("Patient_ID")
            df = df.sort_index()

            # save df in self._data
            self.save_df(df_type, df)


    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df = df[df["site"].notnull()] # Drops all rows with nan values in site column

            # Create our column multiindex
            # Split the genes from the sites, splitting from the right since some genes have hyphens in their names, but the genes and sites are also separated by hyphens
            split_genes = df["site"].str.rsplit("-", n=1, expand=True)
            # hgnc_symbol is a duplicate of split_genes[0], and site is now in split_genes and will be re-inserted differently
            df = df.drop(columns=["hgnc_symbol", "site"])
            df = df.assign(Name=split_genes[0], Site=split_genes[1])
            # Get rid of all lowercase s, t, and y delimeters in the sites
            df["Site"] = df["Site"].str.replace(r"[sty]", r"", regex=True)
            df = df.rename(columns={"refseq_peptide": "Database_ID"})
            df = df.set_index(["Name", "Site", "Peptide", "Database_ID"]) # Turn these columns into a multiindex

            df = df.sort_index()
            df = df.transpose()

            # save df in self._data
            self.save_df(df_type, df)


    def load_proteomics(self):
        df_type = 'proteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df = df[df["hgnc_symbol"].notnull()] # Drops all nan values in hgnc_symbol column

            # Create our column multiindex
            df = df.rename(columns={"hgnc_symbol": "Name", "refseq_peptide": "Database_ID"})
            df = df.set_index(["Name", "Database_ID"])
            df = df.sort_index()
            df = df.transpose()

            # save df in self._data
            self.save_df(df_type, df)


    def load_treatment(self):
        df_type = 'treatment'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep=',', index_col=0)
            df = df.rename(columns={"Participant_ID":"Patient_ID"})
            df = df.set_index("Patient_ID")

            # save df in self._data
            self.save_df(df_type, df)


    def load_definitions(self):
        df_type = 'definitions'
        # TODO: low priority, make sure this actually works later
        if not self._definitions:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            with open(file_path, "r", errors="ignore") as definitions_file:
                for line in definitions_file.readlines():
                    line = line.strip()
                    line = line.split("\t")
                    term = line[0]
                    definition = line[1]
                    self._definitions[term] = definition


    # Override the save_df function from source.py so we can mark normal ov samples
    def save_df(self, datatype, df):

        # Take C prefix off of indices for those samples that have them (tumor samples have C, normal have N)
        df.index = df.index.where(~df.index.str.startswith('C'), df.index.str[1:])

        # Move the prepended N to a .N at the end to match other normal sample labeling in cptac
        df.index = df.index.where(~df.index.str.startswith('N'), df.index.str[1:] + ".N")

        # Drop all OV_QC samples--they're quality control samples not relevant for data analysis
        df = df.drop(index=df.index[df.index.str.startswith("OV_QC")])

        # Inherit the parent event
        super().save_df(datatype, df)


    def how_to_cite(self):
        return super().how_to_cite(cancer_type='high grade serous ovarian cancer', pmid=27372738)
