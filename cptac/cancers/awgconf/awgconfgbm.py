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
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError, InvalidParameterError

class AwgConfGbm(Source):

    def __init__(self, version="latest", no_internet=False):
        """Define which dataframes as are available in the self.load_functions dictionary variable, with names as keys.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        self.valid_versions = ["0.1", "2.0"]

        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        self.data_files = {
            "0.1": {
                "acetylproteomics"    : "acetylome_pnnl_d6.v0.1.20220202.tsv.gz",
                "clinical"            : "clinical_data_core.v0.1.20220202.tsv.gz",
                "phosphoproteomics"   : "phosphoproteome_pnnl_d6.v0.1.20220202.tsv.gz",
                "proteomics"          : "proteome_pnnl_per_gene_d4.v0.1.20220202.tsv.gz",
                "transcriptomics"     : "rnaseq_washu_fpkm_uq.v0.1.20220202.tsv.gz",
                "somatic_mutation"    : "somatic_wes_mutation.v0.1.20220202.maf.gz",
                "CNV"                 : "wgs_somatic_cnv_per_gene.v0.1.20220202.tsv.gz",
            },
            "2.0": {
                "acetylproteomics"            : "acetylome_pnnl_d6.v2.0.20220408.tsv.gz",
                "clinical"                    : "clinical_data_core.v2.0.20220408.tsv.gz",
                "targeted_proteomics"         : ["Direct_SRM_corrected.v2.0.20220408.tsv.gz", "PRISM_SRM_raw.v2.0.20220408.tsv.gz"],
                "targeted_phosphoproteomics"  : "IMAC_SRM_corrected.v2.0.20220408.tsv.gz",
                "metabolomics"                : "metabolome_pnnl.v2.0.20220408.tsv.gz",
                "miRNA"                       : "mirnaseq_mirna_mature_tpm.v2.0.20220408.tsv.gz",
                "phosphoproteomics"           : "phosphoproteome_pnnl_d6.v2.0.20220408.tsv.gz",
                "lipidomics"                  : ["positive_lipidome_pnnl.v2.0.20220408.tsv.gz", "negative_lipidome_pnnl.v2.0.20220408.tsv.gz"],
                "proteomics"                  : "proteome_pnnl_per_gene_d4.v2.0.20220408.tsv.gz",
                "gene_fusion"                 : "rnaseq_gene_fusion.v2.0.20220408.tsv.gz",
                "transcriptomics"             : "rnaseq_washu_fpkm_uq.v2.0.20220408.tsv.gz",
                "somatic_mutation"            : "somatic_wes_mutation.v2.0.20220408.maf.gz",
                "CNV"                         : "wgs_somatic_cnv_gatk4_per_gene.v2.0.20220408.tsv.gz",
            },
        }

        self.load_functions = {
            'acetylproteomics'  : self.load_acetylproteomics,
            'clinical'          : self.load_clinical,
            'CNV'               : self.load_CNV,
            'phosphoproteomics' : self.load_phosphoproteomics,
            'proteomics'        : self.load_proteomics,
            'somatic_mutation'  : self.load_somatic_mutation,
            'transcriptomics'   : self.load_transcriptomics,
        }

        if version != "0.1":
            self.load_functions["targeted_proteomics"] = self.load_targeted_proteomics
            self.load_functions["targeted_phosphoproteomics"] = self.load_targeted_phosphoproteomics
            self.load_functions["metabolomics"] = self.load_metabolomics
            self.load_functions["miRNA"] = self.load_miRNA
            self.load_functions["lipidomics"] = self.load_lipidomics
            self.load_functions["gene_fusion"] = self.load_gene_fusion


        # Call the parent class __init__ function
        super().__init__(cancer_type="gbm", source='awgconf', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)


    def load_acetylproteomics(self):
        df_type = 'acetylproteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')

            # Rename columns for consistency across cptac
            df = df.rename(columns={
                "symbol": "Name",
                "site": "Site",
                "peptide": "Peptide",
                "refseq_id": "Database_ID"
            })
            # Multiindex based on renamed columns
            df = df.set_index(["Name", "Site", "Peptide", "Database_ID"])

            df = df.transpose()
            df.index.name="Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_clinical(self):
        df_type = 'clinical'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep='\t')

            # Save all sample ids that belong to the discovery cohort so they can be removed from the other dataframes
            dfd = df.loc[df['cohort'] == 'Discovery']
            dfd = dfd["preferred_sample_name"].to_frame()
            self._helper_tables["discovery_cohort_ids"] = dfd

            df = df.loc[df['cohort'] != 'Discovery']
            # Add sample tumor normal column
            df.insert(len(df.columns),'Sample_Tumor_Normal',"Tumor",)
            df.loc[df['sample_type'].str.contains("Normal"), 'Sample_Tumor_Normal'] = "Normal"

            # Fill missing tumor_occurrence_sequence info for cptac samples
            df.loc[df['sample_type'].str.contains("CPTAC") & df['preferred_sample_name'].str.contains("TP"), 'tumor_occurrence_sequence'] = "1_primary"
            df.loc[df['sample_type'].str.contains("CPTAC") & df['preferred_sample_name'].str.contains("NAT"), 'tumor_occurrence_sequence'] = "0_normal"

            df = df.set_index("preferred_sample_name")
            df = df.sort_index()
            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_targeted_proteomics(self):
        df_type = 'targeted_proteomics'
        if df_type not in self._data:
             # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path_list = self.locate_files(df_type)

            for file_path in file_path_list:
                path_elements = file_path.split(os.sep) # Get a list of the levels of the path
                file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below
                df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name, so we don't include the version

                if df_name == "Direct_SRM_corrected":
                    df_direct = pd.read_csv(file_path, sep='\t')
                    df_direct = df_direct.rename(columns={
                        "protein_name": "Name",
                        "peptide_seq": "Peptide",
                    })
                    df_direct = df_direct.set_index(["Name", "Peptide"])
                    df_direct = df_direct.transpose()
                    df_direct = df_direct.sort_index()
                    df_direct.index.name = "Patient_ID"

                if df_name == "PRISM_SRM_raw":
                    df_prism = pd.read_csv(file_path, sep='\t')
                    df_prism = df_prism.rename(columns={
                        "protein_name": "Name",
                        "peptide_seq": "Peptide",
                    })
                    df_prism = df_prism.set_index(["Name", "Peptide"])
                    df_prism = df_prism.transpose()
                    df_prism = df_prism.sort_index()
                    df_prism.index.name = "Patient_ID"

            # Combine Direct and PRISM SRM tables
            # These do not contain refseq ids
            targeted_proteomics = df_prism.join(df_direct, how="outer")
            # Add refseq_ids (saved previously in prot_refseq by load_proteomics)
            self.load_proteomics()
            targeted_proteomics = targeted_proteomics.transpose()
            targeted_proteomics = targeted_proteomics.join(self._helper_tables["prot_refseq"], how="left")
            del self._helper_tables["prot_refseq"]
            targeted_proteomics = targeted_proteomics.reset_index()
            levels = ["Name", "Peptide", "Database_ID"]
            targeted_proteomics = targeted_proteomics.set_index(levels)
            # This next line is basically transpose but it doesn't mess up the multiindex
            targeted_proteomics = targeted_proteomics.stack().unstack(levels)

            # save df in self._data
            self.save_df(df_type, targeted_proteomics)


    def load_targeted_phosphoproteomics(self):
        df_type = 'targeted_phosphoproteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df = df.drop(columns=["site"])
            df = df.rename(columns={
                "protein_name": "Name",
                "peptide_modified_seq": "Peptide",
                "refseq_id": "Database_ID"
            })
            df = df.set_index(["Name", "Peptide", "Database_ID"])
            df = df.transpose()
            df = df.sort_index()
            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_metabolomics(self):
        df_type = 'metabolomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df = df.set_index("Metabolite")
            df = df.transpose()
            # Name columns for consistency
            df.columns.name = "Name"
            df.index.name="Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_miRNA(self):
        df_type = 'miRNA'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df = df.drop(columns=["Alias", "Derives_from"])
            df = df.rename(columns={'ID': "Database_ID"})
            df = df.set_index(["Name", "Database_ID"])
            df = df.transpose()
            df.index.name="Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            # Rename columns for consistency across cptac
            df = df.rename(columns={
                "symbol": "Name",
                "site": "Site",
                "peptide": "Peptide",
                "refseq_id": "Database_ID"
            })
            # Multiindex based on renamed columns in cptac order
            df = df.set_index(["Name", "Site", "Peptide", "Database_ID"])

            df = df.transpose()
            df.index.name="Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_proteomics(self):
        df_type = 'proteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            # Save refseq information for use in targeted_proteomics
            df_ref = df[["gene", "refseq_id"]]
            df_ref = df_ref.rename(columns={
                "gene": "Name",
                "refseq_id": "Database_ID",
            })
            df_ref = df_ref.set_index("Name")
            self._helper_tables["prot_refseq"] = df_ref

            # Rename columns and multiindex
            df = df.rename(columns={"gene": "Name", 'refseq_id': "Database_ID"})
            df = df.set_index(["Name", "Database_ID"])
            df = df.transpose()
            df.index.name="Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_gene_fusion(self):
        df_type = 'gene_fusion'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df = df.set_index("preferred_sample_name")
            df.columns.name = "Name"
            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_transcriptomics(self):
        df_type = 'transcriptomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            # May need to remove discovery cohort information from this table
            columns_to_drop = ['gene_type', 'gene_status', 'havana_gene', 'full_length', 'exon_length', 'exon_num']
            df = df.drop(columns=columns_to_drop)
            df = df.rename(columns={
                "gene_name": "Name",
                "gene_id": "Database_ID",
            })
            df = df.set_index(["Name", "Database_ID"])
            df = df.transpose()
            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_somatic_mutation(self):
        df_type = 'somatic_mutation'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep='\t')

            if self.version == "0.1":
                # We don't need any of the other columns
                df = df[["preferred_sample_name", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short"]]
                df = df.rename(columns={
                    "preferred_sample_name": "Patient_ID",
                    'Hugo_Symbol': "Gene",
                    "Variant_Classification": "Mutation",
                    "HGVSp_Short": "Location"
                })
            elif self.version == "2.0":
                df = df[["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short"]]
                df = df.rename(columns={"Tumor_Sample_Barcode": "Patient_ID", 'Hugo_Symbol': "Gene", "Variant_Classification": "Mutation", "HGVSp_Short": "Location"})
                df["Patient_ID"] = df["Patient_ID"].str.rstrip("_T")

            df = df.set_index("Patient_ID")
            df = df.sort_values(by=["Patient_ID","Gene"])

            # save df in self._data
            self.save_df(df_type, df)


    def load_CNV(self):
        df_type = 'CNV'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep='\t')

            if self.version == "0.1":
                df = df.drop(columns=["Start","End", "Chr"])

            df = df.set_index("Gene")
            df = df.sort_index()
            df = df.transpose()
#             df.index.name = "Patient_ID"
#             df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)


    def load_lipidomics(self):
        df_type = 'lipidomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path_list = self.locate_files(df_type)

            for file_path in file_path_list:
                path_elements = file_path.split(os.sep) # Get a list of the levels of the path
                file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below
                df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name, so we don't include the version

                if df_name == "negative_lipidome_pnnl":
                    df_negative = pd.read_csv(file_path, sep='\t', index_col=0)
                    df_negative = df_negative.transpose()
                    df_negative.columns.name = "Name"
                    df_negative.index.name="Patient_ID"
                    df_negative = df_negative.add_suffix("_negative")

                if df_name == "positive_lipidome_pnnl":
                    df_positive = pd.read_csv(file_path, sep='\t', index_col=0)
                    df_positive = df_positive.transpose()
                    df_positive.columns.name = "Name"
                    df_positive.index.name="Patient_ID"
                    df_positive = df_positive.add_suffix("_positive")

            # Combine positive and negative lipidomics tables
            df_combined = df_positive.join(df_negative, how="outer")
            df_combined.index.name = "Patient_ID"
            df_combined.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df_combined)


    # Override the save_df function from source.py so we can do additional formatting and give an access warning the first time the data is used
    def save_df(self, datatype, df):

        # Delete discovery cohort samples from all data types
        if datatype != 'clinical':
            self.load_clinical() # loads self._helper_tables["discovery_cohort_ids"]
        dfd = self._helper_tables["discovery_cohort_ids"]
        df = df[~df.index.isin(dfd["preferred_sample_name"])]

        # Function to remove suffix "-NAT" from Patient_IDs and add .N
        # This is necessary to sort tumor and normal samples in functions written in (ie _tumor_only or _normal_only)
        df.index = df.index.where(~df.index.str.endswith('-NAT'), df.index.str[:-4] + ".N")

        if self._data == {}:
            # Print password access only warning
            warnings.warn("The GBM confirmatory data is currently strictly reserved for CPTAC investigators. Otherwise, you are not authorized to access these data. Additionally, even after these data become publicly available, they will be subject to a publication embargo (see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter cptac.embargo() to open the webpage for more details).", PublicationEmbargoWarning, stacklevel=2)

        # Inherit the parent event
        super().save_df(datatype, df)
