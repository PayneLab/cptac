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

from asyncio import CancelledError
import pandas as pd
import numpy as np
import os
import warnings
import datetime

from cptac.cancers.source import Source
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError, InvalidParameterError

class AwgConfUcec(Source):

    def __init__(self, version="latest", no_internet=False):
        """Define which dataframes as are available in the self.load_functions dictionary variable, with names as keys.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest datafreeze. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        self.valid_versions = ["1.0", "1.1", "1.2", "2.0", "2.0.1"]

        self.data_files = {
            "1.0": {
            "acetylproteomics_gene"    : "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v1.0.cct.gz",
            "acetylproteomics"         : "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
            #"not_used" : "UCEC_confirmatory_Direct_SRM_tumor_v1.0.cct.gz", #SRM not to be included in 1.0
            #"not_used" : "UCEC_confirmatory_IMAC_SRM_tumor_v1.0.cct.gz",
            "clinical"                 : "UCEC_confirmatory_meta_table_v1.0.xlsx",
            "methylation"              : "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v1.0.cct.gz",
            "miRNA"                    : "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v1.0.cct.gz",
            #"not_used" : "UCEC_confirmatory_nglycoform-site_ratio_median_polishing_log2_tumor_normal_v1.0.cct.gz",
            "phosphoproteomics_gene"   : "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
            "phosphoproteomics"        : "UCEC_confirmatory_phospho_site_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
            "proteomics"               : "UCEC_confirmatory_proteomics_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
            "circular_RNA"             : "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v1.0.cct.gz",
            "gene_fusion"              : "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v1.0.txt.gz",
            "transcriptomics"          : "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v1.0.cct.gz",
            #"not_used" : "UCEC_confirmatory_RNAseq_isoform_FPKM_removed_circRNA_log2(x+1)_tumor_normal_v1.0.cct.gz",
            "CNV_gistic"               : "UCEC_confirmatory_WES_cnv_gistic_thresholded_tumor_v1.0.cct.gz",
            "CNV_log2ratio"            : "UCEC_confirmatory_WES_cnv_log2_ratio_tumor_v1.0.cct.gz",
            "somatic_mutation_binary"  : "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.0.cbt.gz",
            "somatic_mutation"         : "UCEC_confirmatory_WES_somatic_mutation_v1.0.maf.gz",
            #"not_used" : "UCEC_confirmatory_WGS_SV_tumor_v1.0.txt.gz" #structural_variation - not to be included in 1.0
            },
        "1.1": {
            "acetylproteomics_gene"      : "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
            "acetylproteomics"           : "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
            "clinical"                   : "UCEC_confirmatory_meta_table_v1.1.xlsx",
            "methylation"                : "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v1.1.cct.gz",
            "miRNA"                      : "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v1.1.cct.gz",
            "phosphoproteomics_gene"     : "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
            "phosphoproteomics"          : "UCEC_confirmatory_phospho_site_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
            "proteomics"                 : "UCEC_confirmatory_proteomics_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
            "circular_RNA"               : "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v1.1.cct.gz",
            "gene_fusion"                : "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v1.1.txt.gz",
            "transcriptomics"            : "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v1.1.cct.gz",
            "targeted_proteomics"        : ["UCEC_confirmatory_SRM_Direct_tumor_v1.1.cct.gz", "UCEC_confirmatory_SRM_PRISM_tumor_v1.1.cct.gz"],
            "targeted_phosphoproteomics" : "UCEC_confirmatory_SRM_IMAC_tumor_v1.1.cct.gz",
            "CNV_gistic"                 : "UCEC_confirmatory_WES_cnv_gistic_thresholded_tumor_v1.1.cct.gz",
            "CNV_log2ratio"              : "UCEC_confirmatory_WES_cnv_log2_ratio_tumor_v1.1.cct.gz",
            "somatic_mutation_binary"    : "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.1.cbt.gz",
            "somatic_mutation"           : "UCEC_confirmatory_WES_somatic_mutation_v1.1.maf.gz",
            #"not_used" : "UCEC_confirmatory_nglycoform-site_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
            #"not_used" : "UCEC_confirmatory_RNAseq_isoform_FPKM_removed_circRNA_log2(x+1)_tumor_normal_v1.1.cct.gz",
            #"not_used" : "UCEC_confirmatory_WGS_SV_tumor_v1.1.txt.gz",
            },
        "1.2": {
            "clinical"                   : "UCEC_confirmatory_meta_table_v1.2.xlsx",
            "targeted_proteomics"        : ["UCEC_confirmatory_SRM_Direct_tumor_v1.2.cct.gz", "UCEC_confirmatory_SRM_PRISM_tumor_v1.2.cct.gz"],
            "targeted_phosphoproteomics" : "UCEC_confirmatory_SRM_IMAC_tumor_v1.2.cct.gz",
            "circular_RNA"               : "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v1.2.cct.gz",
            "transcriptomics"            : "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v1.2.cct.gz",
            "gene_fusion"                : "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v1.2.txt.gz",
#             "not_used" : "UCEC_confirmatory_RNAseq_isoform_FPKM_removed_circRNA_log2(x+1)_tumor_normal_v1.2.cct.gz",
            "CNV_gistic"                 : "UCEC_confirmatory_WGS_cnv_gistic_thresholded_tumor_v1.2.cct.gz",
            "CNV_log2ratio"              : "UCEC_confirmatory_WGS_cnv_log2_ratio_tumor_v1.2.cct.gz",
            "somatic_mutation_binary"    : "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.2.cbt.gz",
            "somatic_mutation"           : "UCEC_confirmatory_WES_somatic_mutation_v1.2.maf.gz",
#             "not_used" : "UCEC_confirmatory_WGS_SV_tumor_v1.2.txt.gz",
            "acetylproteomics_gene"      : "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
            "acetylproteomics"           : "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
            "methylation"                : "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v1.2.cct.gz",
            "miRNA"                      : "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v1.2.cct.gz",
#             "not_used" : "UCEC_confirmatory_nglycoform-site_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
            "phosphoproteomics_gene"     : "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
            "phosphoproteomics"          : "UCEC_confirmatory_phospho_site_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
            "proteomics"                 : "UCEC_confirmatory_proteomics_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
            },
        "2.0": {
            "acetylproteomics_gene"      : "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "acetylproteomics"           : "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "clinical"                   : "UCEC_confirmatory_meta_table_v2.0.xlsx",
            "methylation"                : "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v2.0.cct.gz",
            "miRNA"                      : "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v2.0.cct.gz",
            "phosphoproteomics_gene"     : "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "phosphoproteomics"          : "UCEC_confirmatory_phospho_site_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "proteomics"                 : "UCEC_confirmatory_proteomics_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "circular_RNA"               : "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v2.0.cct.gz",
            "gene_fusion"                : "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v2.0.txt.gz",
            "transcriptomics"            : "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v2.0.cct.gz",
            "targeted_proteomics"        : ["UCEC_confirmatory_SRM_Direct_tumor_v2.0.cct.gz", "UCEC_confirmatory_SRM_PRISM_tumor_v2.0.cct.gz"],
            "targeted_phosphoproteomics" : "UCEC_confirmatory_SRM_IMAC_tumor_v2.0.cct.gz",
#             "not_used" : "UCEC_confirmatory_WES_somatic_mutation_category_level_V1.2.txt.gz",
            "somatic_mutation_binary"    : "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.2.cbt.gz",
            "somatic_mutation"           : "UCEC_confirmatory_WES_somatic_mutation_v2.0.maf.gz",
            "CNV_gistic"                 : "UCEC_confirmatory_WGS_cnv_gistic_thresholded_tumor_v2.0.cct.gz",
            "CNV_log2ratio"              : "UCEC_confirmatory_WGS_cnv_log2_ratio_tumor_v2.0.cct.gz",
#             "not_used" : "UCEC_confirmatory_WGS_SV_tumor_v2.0.txt.gz",
            },
        "2.0.1": {
            "acetylproteomics_gene"      : "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "acetylproteomics"           : "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "clinical"                   : "UCEC_confirmatory_meta_table_v2.0.1.xlsx",
            "methylation"                : "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v2.0.cct.gz",
            "miRNA"                      : "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v2.0.cct.gz",
            "phosphoproteomics_gene"     : "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "phosphoproteomics"          : "UCEC_confirmatory_phospho_site_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "proteomics"                 : "UCEC_confirmatory_proteomics_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "circular_RNA"               : "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v2.0.cct.gz",
            "gene_fusion"                : "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v2.0.txt.gz",
            "transcriptomics"            : "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v2.0.cct.gz",
            "targeted_proteomics"        : ["UCEC_confirmatory_SRM_Direct_tumor_v2.0.cct.gz", "UCEC_confirmatory_SRM_PRISM_tumor_v2.0.cct.gz"],
            "targeted_phosphoproteomics" : "UCEC_confirmatory_SRM_IMAC_tumor_v2.0.cct.gz",
            "somatic_mutation_binary"    : "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.2.cbt.gz",
            "somatic_mutation"           : "UCEC_confirmatory_WES_somatic_mutation_v2.0.maf.gz",
            "CNV_gistic"                 : "UCEC_confirmatory_WGS_cnv_gistic_thresholded_tumor_v2.0.cct.gz",
            "CNV_log2ratio"              : "UCEC_confirmatory_WGS_cnv_log2_ratio_tumor_v2.0.cct.gz",
            },
        }

        self.load_functions = {
            'acetylproteomics_gene' : self.load_acetylproteomics_gene,
            'acetylproteomics' : self.load_acetylproteomics,
            'clinical' : self.load_clinical,
            'methylation' : self.load_methylation,
            'miRNA' : self.load_miRNA,
            'phosphoproteomics_gene' : self.load_phosphoproteomics_gene,
            'phosphoproteomics' : self.load_phosphoproteomics,
            'proteomics' : self.load_proteomics,
            'circular_RNA' : self.load_circular_RNA,
            'gene_fusion' : self.load_gene_fusion,
            'transcriptomics' : self.load_transcriptomics,
            'targeted_proteomics' : self.load_targeted_proteomics,
            'targeted_phosphoproteomics' : self.load_targeted_phosphoproteomics,
            'somatic_mutation_binary' : self.load_somatic_mutation_binary,
            'somatic_mutation' : self.load_somatic_mutation,
            'CNV_gistic' : self.load_CNV_gistic,
            'CNV_log2ratio' : self.load_CNV_log2ratio,
        }

        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        # Call the parent class __init__ function
        super().__init__(cancer_type="ucec", source='awgconf', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)


    def load_acetylproteomics_gene(self):
        df_type = 'acetylproteomics_gene'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df = df.sort_index()
            df.index.name = "Patient_ID"
            df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)


    def load_acetylproteomics(self):
        df_type = 'acetylproteomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.reset_index()
            df[['Name','Database_ID','Site']] = df.idx.str.split("@", expand=True)
            df['Site'] = df['Site'].str.rsplit('-',1,expand=True)[1]
            df = df.set_index(["Name", "Site", "Database_ID"])
            df = df.drop(columns=["idx"])
            df = df.transpose()
            df = df.sort_index()
            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_clinical(self):
        df_type = 'clinical'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_excel(file_path)
            df = df.set_index("Idx")
            df = df.rename({'Group': 'Proteomics_Tumor_Normal'}, axis=1)
            df.index.name = "Patient_ID"
            df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)


    def load_methylation(self):
        df_type = 'methylation'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0, na_values='   NA')
            df = df.transpose()
            df.index.name = "Patient_ID"
            df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)

    def load_miRNA(self):
        df_type = 'miRNA'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df = df.sort_index()
            df.index.name = "Patient_ID"
            df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)


    def load_phosphoproteomics_gene(self):
        df_type = 'phosphoproteomics_gene'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df = df.sort_index()
            df.index.name = "Patient_ID"
            df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)


    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.reset_index()
            df[['Name','Database_ID','Site']] = df.idx.str.split("@", expand=True)
            df['Site'] = df['Site'].str.rsplit('-',1,expand=True)[1]
            df = df.set_index(["Name", "Site", "Database_ID"])
            df = df.drop(columns=["idx"])
            df = df.transpose()
            df = df.sort_index()
            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)

    def load_proteomics(self):
        df_type = 'proteomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df = df.sort_index()
            df.index.name = "Patient_ID"
            df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)


    def load_circular_RNA(self):
        df_type = 'circular_RNA'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df.index.name = "Patient_ID"
            df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)


    def load_gene_fusion(self):
        df_type = 'gene_fusion'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.reset_index()
            df = df.set_index("Sample")
            df.index.name = "Patient_ID"
            df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)


    def load_transcriptomics(self):
        df_type = 'transcriptomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df = df.sort_index()
            df.index.name = "Patient_ID"
            df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)


    def load_targeted_proteomics(self):
        df_type = 'targeted_proteomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path_list = self.locate_files(df_type)
            for file_path in file_path_list:
                path_elements = file_path.split(os.sep) # Get a list of the levels of the path
                file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

                # Load and parse information from PRISM and Direct files
                if "UCEC_confirmatory_SRM_Direct_tumor" in file_name:
                    df_direct = pd.read_csv(file_path, sep='\t')
                    df_direct[['Name','Peptide']] = df_direct['idx'].str.rsplit("-", 1, expand=True)
                    df_direct = df_direct.set_index(["Name", "Peptide"])
                    df_direct = df_direct.drop(columns=["idx"])
                    df_direct = df_direct.transpose()
                    df_direct = df_direct.sort_index()
                    df_direct.index.name = "Patient_ID"

                elif "UCEC_confirmatory_SRM_PRISM_tumor" in file_name:
                    df_prism = pd.read_csv(file_path, sep='\t')
                    df_prism[['Name','Peptide']] = df_prism['idx'].str.rsplit("-", 1, expand=True)
                    df_prism = df_prism.set_index(["Name", "Peptide"])
                    df_prism = df_prism.drop(columns=["idx"])
                    df_prism = df_prism.transpose()
                    df_prism = df_prism.sort_index()
                    df_prism.index.name = "Patient_ID"

            # Merge once we have both
            df_combined = pd.concat([df_direct, df_prism])
            df_combined.index.name = "Patient_ID"
            df_combined.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df_combined)


    def load_targeted_phosphoproteomics(self):
        df_type = 'targeted_phosphoproteomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df.at[0,'idx'] = "FPSS[+80]PLRIPGGNIY[+80]ISPLK"
            df['Name'] = "RB1"
            df = df.rename(columns={"idx":"Peptide"})
            df = df.set_index(["Name", "Peptide"])
            df = df.transpose()
            df = df.sort_index()
            df.columns.name = "Name"
            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_somatic_mutation_binary(self):
        df_type = 'somatic_mutation_binary'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df.index.name = "Patient_ID"
            df.columns.name = "Name"

            # save df in self._data
            self.save_df(df_type, df)


    def load_somatic_mutation(self):
        df_type = 'somatic_mutation'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', dtype={88:object})
            df = df[['Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification','HGVSp_Short']]
            df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].apply(lambda s: s[:-2])
            df = df.rename(columns={
                "Tumor_Sample_Barcode": "Patient_ID",
                "Hugo_Symbol":"Gene",
                "Variant_Classification":"Mutation",
                "HGVSp_Short":"Location"})
            df = df.sort_values(by=["Patient_ID", "Gene"])
            df.columns.name = "Name"
            df = df.set_index("Patient_ID")

            # save df in self._data
            self.save_df(df_type, df)


    def load_CNV_gistic(self):
        df_type = 'CNV_gistic'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df[['Name','Chromosome']] = df.idx.str.split("|", expand=True)
            df = df.set_index(["Name"])
            df = df.drop(columns=["idx", "Chromosome"])
            df = df.transpose()
            df = df.sort_index()
            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_CNV_log2ratio(self):
        df_type = 'CNV_log2ratio'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df[['Name','Chromosome']] = df.idx.str.split("|", expand=True)
            df = df.set_index(["Name"])
            df = df.drop(columns=["idx", "Chromosome"])
            df = df.transpose()
            df = df.sort_index()
            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    # Override the save_df function from source.py so we can give an access warning the first time the data is used
    def save_df(self, datatype, df):

        # Mark normal samples in each datatype
        df.index = df.index.where(~df.index.str.endswith('-A'), df.index.str[:-2] + ".N")
        df.index = df.index.where(~df.index.str.startswith('N'), df.index.str[:] + ".N")

        if self._data == {}:
            # Print password access only warning
            warnings.warn("The endometrial confirmatory cohort data is currently strictly reserved for CPTAC investigators. Otherwise, you are not authorized to access these data. Additionally, even after these data become publicly available, they will be subject to a publication embargo (see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter cptac.embargo() to open the webpage for more details).", PublicationEmbargoWarning, stacklevel=2)

        # Inherit the parent event
        super().save_df(datatype, df)
