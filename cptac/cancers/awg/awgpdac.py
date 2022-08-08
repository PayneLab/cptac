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
import datetime

from cptac.cancers.source import Source
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError

class AwgPdac(Source):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        self.valid_versions = ["1.0"]

        self.data_files = {
            "1.0": {
            "clinical"          : "clinical_table_140.tsv.gz",
            "miRNA"             : ["microRNA_TPM_log2_Normal.cct.gz", "microRNA_TPM_log2_Tumor.cct.gz"],
            "derived_molecular" : "meta_table_140.tsv.gz",
            "transcriptomics"   : ["mRNA_RSEM_UQ_log2_Normal.cct.gz", "mRNA_RSEM_UQ_log2_Tumor.cct.gz"],
            "somatic_mutation"  : "PDAC_mutation.maf.gz",
            "phosphoproteomics" : ["phosphoproteomics_site_level_MD_abundance_normal.cct.gz", "phosphoproteomics_site_level_MD_abundance_tumor.cct.gz"],
            "proteomics"        : ["proteomics_gene_level_MD_abundance_normal.cct.gz", "proteomics_gene_level_MD_abundance_tumor.cct.gz"],
            "gene_fusion"       : ["RNA_fusion_unfiltered_normal.tsv.gz","RNA_fusion_unfiltered_tumor.tsv.gz"],
            "CNV"               : "SCNA_log2_gene_level.cct.gz"},
        }

        self.load_functions = {
            'clinical'                : self.load_clinical,
            'CNV'                     : self.load_CNV,
            'derived_molecular'       : self.load_derived_molecular,
            'gene_fusion'             : self.load_gene_fusion,
            'miRNA'                   : self.load_miRNA,
            'phosphoproteomics'       : self.load_phosphoproteomics,
            'proteomics'              : self.load_proteomics,
            'somatic_mutation'        : self.load_somatic_mutation,
            'transcriptomics'         : self.load_transcriptomics,
        }

        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        # Call the parent class __init__ function
        super().__init__(cancer_type="pdac", source='awg', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)


    def load_clinical(self):
        df_type = 'clinical'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.rename_axis("Patient_ID", axis="index")
            df = df.sort_index()
            df.columns.name = "Name"
            df["Sample_Tumor_Normal"] = "Tumor"

            # save df in self._data
            self.save_df(df_type, df)


    def load_derived_molecular(self):
        df_type = 'derived_molecular'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.sort_index()
            df.index.name = "Patient_ID"
            df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)


    def load_somatic_mutation(self):
        df_type = 'somatic_mutation'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df = df[["Hugo_Symbol", "Variant_Classification", "HGVSp_Short", "Tumor_Sample_Barcode"]]
            df = df.rename({"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
            df = df.sort_values(by=["Patient_ID", "Gene"])
            df = df.set_index("Patient_ID")
            # remove _T and similar tags from end of string
            df = df.rename(index=lambda s: s[:-2])
            df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)


    def load_CNV(self):
        df_type = 'CNV'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df.columns.name = "Name"
            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_miRNA(self):
        df_type = 'miRNA'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path_list = self.locate_files(df_type)

            for file_path in file_path_list:
                path_elements = file_path.split(os.sep) # Get a list of the levels of the path
                file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

                if file_name == "microRNA_TPM_log2_Normal.cct.gz":
                    df_normal = pd.read_csv(file_path, sep='\t', index_col=0)
                    df_normal = df_normal.sort_index()
                    df_normal = df_normal.transpose()
                    # append .N to patient ids from the normal table
                    df_normal = df_normal.rename(index=lambda s: s + ".N")

                if file_name == "microRNA_TPM_log2_Tumor.cct.gz":
                    df_tumor = pd.read_csv(file_path, sep='\t', index_col=0)
                    df_tumor = df_tumor.sort_index()
                    df_tumor = df_tumor.transpose()

            # merge tumor and normal data
            df_combined = pd.concat([df_normal, df_tumor])
            df_combined.index.name = "Patient_ID"
            df_combined.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df_combined)


    def load_transcriptomics(self):
        df_type = 'transcriptomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path_list = self.locate_files(df_type)

            for file_path in file_path_list:
                path_elements = file_path.split(os.sep) # Get a list of the levels of the path
                file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

                if file_name == "mRNA_RSEM_UQ_log2_Normal.cct.gz":
                    df_normal = pd.read_csv(file_path, sep='\t', index_col=0)
                    df_normal = df_normal.sort_index()
                    df_normal = df_normal.transpose()
                    # append .N to patient ids from the normal table
                    df_normal = df_normal.rename(index=lambda s: s + ".N")

                if file_name == "mRNA_RSEM_UQ_log2_Tumor.cct.gz":
                    df_tumor = pd.read_csv(file_path, sep='\t', index_col=0)
                    df_tumor = df_tumor.sort_index()
                    df_tumor = df_tumor.transpose()

            # merge tumor and normal data
            df_combined = pd.concat([df_normal, df_tumor])
            df_combined.index.name = "Patient_ID"
            df_combined.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df_combined)


    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path_list = self.locate_files(df_type)

            for file_path in file_path_list:
                path_elements = file_path.split(os.sep) # Get a list of the levels of the path
                file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

                if file_name == "phosphoproteomics_site_level_MD_abundance_normal.cct.gz":
                    df_normal = pd.read_csv(file_path, sep='\t')
                    column_split = df_normal["Index"].str.rsplit("_", n=1, expand=True)
                    df_normal = df_normal.assign(
                        Site = column_split[1],
                        Database_ID = column_split[0]
                    )
                    df_normal = df_normal.drop(columns="Index")
                    df_normal = df_normal.rename(columns={"Gene":"Name"})
                    df_normal = df_normal.set_index(["Name", "Site", "Peptide", "Database_ID"])
                    df_normal = df_normal.sort_index()
                    df_normal = df_normal.transpose()
                    # append .N to patient ids from the normal table
                    df_normal = df_normal.rename(index=lambda s: s + ".N")

                if file_name == "phosphoproteomics_site_level_MD_abundance_tumor.cct.gz":
                    df_tumor = pd.read_csv(file_path, sep='\t')
                    column_split = df_tumor["Index"].str.rsplit("_", n=1, expand=True)
                    df_tumor = df_tumor.assign(
                        Site = column_split[1],
                        Database_ID = column_split[0]
                    )
                    df_tumor = df_tumor.drop(columns="Index")
                    df_tumor = df_tumor.rename(columns={"Gene":"Name"})
                    df_tumor = df_tumor.set_index(["Name", "Site", "Peptide", "Database_ID"])
                    df_tumor = df_tumor.sort_index()
                    df_tumor = df_tumor.transpose()

            # merge tumor and normal data
            df_combined = pd.concat([df_normal, df_tumor])
            df_combined.index.name = "Patient_ID"
            #df_combined.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df_combined)


    def load_proteomics(self):
        df_type = 'proteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path_list = self.locate_files(df_type)

            for file_path in file_path_list:
                path_elements = file_path.split(os.sep) # Get a list of the levels of the path
                file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

                if file_name == "proteomics_gene_level_MD_abundance_normal.cct.gz":
                    df_normal = pd.read_csv(file_path, sep='\t', index_col=0)
                    df_normal = df_normal.sort_index()
                    df_normal = df_normal.transpose()
                    # append .N to patient ids from the normal table
                    df_normal = df_normal.rename(index=lambda s: s + ".N")

                if file_name == "proteomics_gene_level_MD_abundance_tumor.cct.gz":
                    df_tumor = pd.read_csv(file_path, sep='\t', index_col=0)
                    df_tumor = df_tumor.sort_index()
                    df_tumor = df_tumor.transpose()

            # merge tumor and normal data
            df_combined = pd.concat([df_normal, df_tumor])
            df_combined.index.name = "Patient_ID"
            df_combined.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df_combined)


    def load_gene_fusion(self):
        df_type = 'gene_fusion'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path_list = self.locate_files(df_type)

            for file_path in file_path_list:
                path_elements = file_path.split(os.sep) # Get a list of the levels of the path
                file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

                if file_name == "RNA_fusion_unfiltered_normal.tsv.gz":
                    df_normal = pd.read_csv(file_path, sep='\t', index_col=0)
                    df_normal = df_normal.rename(columns={"Sample": "Patient_ID"})
                    df_normal = df_normal.set_index("Patient_ID")
                    # append .N to patient ids from the normal table
                    df_normal = df_normal.rename(index=lambda s: s + ".N")

                if file_name == "RNA_fusion_unfiltered_tumor.tsv.gz":
                    df_tumor = pd.read_csv(file_path, sep='\t', index_col=0)
                    df_tumor = df_tumor.rename(columns={"Sample": "Patient_ID"})
                    df_tumor = df_tumor.set_index("Patient_ID")

            # merge tumor and normal data
            df_combined = pd.concat([df_normal, df_tumor])
            df_combined.index.name = "Patient_ID"
            df_combined.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df_combined)


    def how_to_cite(self):
        return super().how_to_cite(cancer_type='pancreatic ductal adenocarcinoma', pmid='34534465')
