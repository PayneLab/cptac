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

from .dataset import Dataset
from .dataframe_tools import *
from .exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError, InvalidParameterError

class UcecConf(Dataset):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the dataframes as values in the self._data dict variable, with names as keys, and format them properly.
        
        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        valid_versions = ["1.0", "1.1", "1.2", "2.0", "2.0.1"]

        data_files = {
            "1.0": [
            "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
            #"UCEC_confirmatory_Direct_SRM_tumor_v1.0.cct.gz", #SRM not to be included in 1.0
            #"UCEC_confirmatory_IMAC_SRM_tumor_v1.0.cct.gz",
            "UCEC_confirmatory_meta_table_v1.0.xlsx",
            "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v1.0.cct.gz",
            "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v1.0.cct.gz",
            #"UCEC_confirmatory_nglycoform-site_ratio_median_polishing_log2_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_phospho_site_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_proteomics_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v1.0.txt.gz",
            "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v1.0.cct.gz",
            #"UCEC_confirmatory_RNAseq_isoform_FPKM_removed_circRNA_log2(x+1)_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_WES_cnv_gistic_thresholded_tumor_v1.0.cct.gz",
            "UCEC_confirmatory_WES_cnv_log2_ratio_tumor_v1.0.cct.gz",
            "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.0.cbt.gz",
            "UCEC_confirmatory_WES_somatic_mutation_v1.0.maf.gz",
            #"UCEC_confirmatory_WGS_SV_tumor_v1.0.txt.gz" #structural_variation - not to be included in 1.0
            ],
        "1.1": [
            "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
            "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
            "UCEC_confirmatory_meta_table_v1.1.xlsx",
            "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v1.1.cct.gz",
            "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v1.1.cct.gz",
            "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
            "UCEC_confirmatory_phospho_site_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
            "UCEC_confirmatory_proteomics_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
            "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v1.1.cct.gz",
            "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v1.1.txt.gz",
            "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v1.1.cct.gz",
            "UCEC_confirmatory_SRM_Direct_tumor_v1.1.cct.gz",
            "UCEC_confirmatory_SRM_IMAC_tumor_v1.1.cct.gz",
            "UCEC_confirmatory_SRM_PRISM_tumor_v1.1.cct.gz",
            "UCEC_confirmatory_WES_cnv_gistic_thresholded_tumor_v1.1.cct.gz",
            "UCEC_confirmatory_WES_cnv_log2_ratio_tumor_v1.1.cct.gz",
            "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.1.cbt.gz",
            "UCEC_confirmatory_WES_somatic_mutation_v1.1.maf.gz",
            ],
        "1.2": [
            "UCEC_confirmatory_meta_table_v1.2.xlsx",
            "UCEC_confirmatory_SRM_Direct_tumor_v1.2.cct.gz",
            "UCEC_confirmatory_SRM_IMAC_tumor_v1.2.cct.gz",
            "UCEC_confirmatory_SRM_PRISM_tumor_v1.2.cct.gz",
            "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v1.2.cct.gz",
            "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v1.2.cct.gz",
            "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v1.2.txt.gz",
#             "UCEC_confirmatory_RNAseq_isoform_FPKM_removed_circRNA_log2(x+1)_tumor_normal_v1.2.cct.gz",
            "UCEC_confirmatory_WGS_cnv_gistic_thresholded_tumor_v1.2.cct.gz",
            "UCEC_confirmatory_WGS_cnv_log2_ratio_tumor_v1.2.cct.gz",
            "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.2.cbt.gz",
            "UCEC_confirmatory_WES_somatic_mutation_v1.2.maf.gz",
#             "UCEC_confirmatory_WGS_SV_tumor_v1.2.txt.gz",
            "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
            "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
            "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v1.2.cct.gz",
            "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v1.2.cct.gz",
#             "UCEC_confirmatory_nglycoform-site_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
            "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
            "UCEC_confirmatory_phospho_site_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
            "UCEC_confirmatory_proteomics_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
            ],
        "2.0": [
            "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_meta_table_v2.0.xlsx",
            "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v2.0.cct.gz",
            "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_phospho_site_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_proteomics_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v2.0.txt.gz",
            "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_SRM_Direct_tumor_v2.0.cct.gz",
            "UCEC_confirmatory_SRM_IMAC_tumor_v2.0.cct.gz",
            "UCEC_confirmatory_SRM_PRISM_tumor_v2.0.cct.gz",
#             "UCEC_confirmatory_WES_somatic_mutation_category_level_V1.2.txt.gz",
            "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.2.cbt.gz",
            "UCEC_confirmatory_WES_somatic_mutation_v2.0.maf.gz",
            "UCEC_confirmatory_WGS_cnv_gistic_thresholded_tumor_v2.0.cct.gz",
            "UCEC_confirmatory_WGS_cnv_log2_ratio_tumor_v2.0.cct.gz",
#             "UCEC_confirmatory_WGS_SV_tumor_v2.0.txt.gz",
            ],
        "2.0.1": [
            "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_meta_table_v2.0.1.xlsx",
            "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v2.0.cct.gz",
            "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_phospho_site_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_proteomics_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v2.0.txt.gz",
            "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v2.0.cct.gz",
            "UCEC_confirmatory_SRM_Direct_tumor_v2.0.cct.gz",
            "UCEC_confirmatory_SRM_IMAC_tumor_v2.0.cct.gz",
            "UCEC_confirmatory_SRM_PRISM_tumor_v2.0.cct.gz",
            "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.2.cbt.gz",
            "UCEC_confirmatory_WES_somatic_mutation_v2.0.maf.gz",
            "UCEC_confirmatory_WGS_cnv_gistic_thresholded_tumor_v2.0.cct.gz",
            "UCEC_confirmatory_WGS_cnv_log2_ratio_tumor_v2.0.cct.gz",
            ],
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="ucecconf", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below
            
            if file_name in ["UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v1.0.cct.gz", 
                             "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
                             "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
                             "UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",]:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["acetylproteomics_gene"] = df
            
            elif file_name in ["UCEC_confirmatory_acetyl_site_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz", 
                                "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
                                "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
                                "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",]:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.reset_index()
                df[['Name','Database_ID','Site']] = df.idx.str.split("@", expand=True)
                df['Site'] = df['Site'].str.rsplit('-',1,expand=True)[1]
                df = df.set_index(["Name", "Site", "Database_ID"])
                df = df.drop(columns=["idx"])
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                self._data["acetylproteomics"] = df
            
            elif file_name in ["UCEC_confirmatory_meta_table_v1.0.xlsx", 
                               "UCEC_confirmatory_meta_table_v1.1.xlsx",
                               "UCEC_confirmatory_meta_table_v1.2.xlsx",
                               "UCEC_confirmatory_meta_table_v2.0.xlsx",
                               "UCEC_confirmatory_meta_table_v2.0.1.xlsx"]:
                df = pd.read_excel(file_path)
                df.insert(6, "Proteomics_Tumor_Normal", df["Group"])
                df.loc[df['Group'] == 'Enriched_Normal', 'Idx'] = df['Idx'] + '.N'
                df.loc[df['Group'] == 'Adjacent_normal', 'Idx'] = df['Idx'].str[:-2] + '.N'
                df = df.set_index("Idx")
                df.loc[df['Group'] != 'Tumor', 'Group'] = 'Normal'
                df = df.rename({'Group': 'Sample_Tumor_Normal'}, axis=1)
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["clinical"] = df
                
            elif file_name in ["UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v1.0.cct.gz", 
                               "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v1.1.cct.gz",
                               "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v1.2.cct.gz",
                               "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v2.0.cct.gz",]:
                df = pd.read_csv(file_path, sep='\t', index_col=0, na_values='   NA')
                df = df.transpose()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["methylation"] = df
            
            elif file_name in ["UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v1.0.cct.gz",
                               "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v1.1.cct.gz",
                               "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v1.2.cct.gz",
                               "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v2.0.cct.gz",]:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["miRNA"] = df
                
            elif file_name in ["UCEC_confirmatory_phospho_gene_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
                               "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
                               "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
                               "UCEC_confirmatory_phospho_gene_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",]:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["phosphoproteomics_gene"] = df
                
            elif file_name in ["UCEC_confirmatory_phospho_site_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
                               "UCEC_confirmatory_phospho_site_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
                               "UCEC_confirmatory_phospho_site_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
                               "UCEC_confirmatory_phospho_site_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",]:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.reset_index()
                df[['Name','Database_ID','Site']] = df.idx.str.split("@", expand=True)
                df['Site'] = df['Site'].str.rsplit('-',1,expand=True)[1]
                df = df.set_index(["Name", "Site", "Database_ID"])
                df = df.drop(columns=["idx"])
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                self._data["phosphoproteomics"] = df
                
            elif file_name in ["UCEC_confirmatory_proteomics_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
                               "UCEC_confirmatory_proteomics_ratio_median_polishing_log2_tumor_normal_v1.1.cct.gz",
                               "UCEC_confirmatory_proteomics_ratio_median_polishing_log2_tumor_normal_v1.2.cct.gz",
                               "UCEC_confirmatory_proteomics_ratio_median_polishing_log2_tumor_normal_v2.0.cct.gz",]:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["proteomics"] = df
                
            elif file_name in ["UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v1.0.cct.gz",
                               "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v1.1.cct.gz",
                               "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v1.2.cct.gz",
                               "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v2.0.cct.gz",]:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["circular_RNA"] = df
                
            elif file_name in ["UCEC_confirmatory_RNAseq_gene_fusion_tumor_v1.0.txt.gz", 
                               "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v1.1.txt.gz",
                               "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v1.2.txt.gz",
                               "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v2.0.txt.gz",]:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.reset_index()
                df = df.set_index("Sample")
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["gene_fusion"] = df
                
            elif file_name in ["UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v1.0.cct.gz", 
                               "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v1.1.cct.gz",
                               "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v1.2.cct.gz",
                               "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v2.0.cct.gz",]:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["transcriptomics"] = df
                
            # Targeted proteomics is the direct and PRISM SRM data
            elif file_name in ["UCEC_confirmatory_SRM_Direct_tumor_v1.1.cct.gz",
                               "UCEC_confirmatory_SRM_Direct_tumor_v1.2.cct.gz",
                               "UCEC_confirmatory_SRM_Direct_tumor_v2.0.cct.gz",]:
                df_direct = pd.read_csv(file_path, sep='\t')
                df_direct[['Name','Peptide']] = df_direct['idx'].str.rsplit("-", 1, expand=True)
                df_direct = df_direct.set_index(["Name", "Peptide"])
                df_direct = df_direct.drop(columns=["idx"])
                df_direct = df_direct.transpose()
                df_direct = df_direct.sort_index()
                df_direct.index.name = "Patient_ID"
                
                # Merge if we have both
                if "targeted_proteomics" in self._data:
                    df_prism = self._data["targeted_proteomics"]
                    df_combined = pd.concat([df_direct, df_prism])
                    df_combined.index.name = "Patient_ID"
                    df_combined.columns.name = "Name"
                    self._data["targeted_proteomics"] = df_combined
                else:
                    self._data["targeted_proteomics"] = df_direct
                
            elif file_name in ["UCEC_confirmatory_SRM_PRISM_tumor_v1.1.cct.gz",
                               "UCEC_confirmatory_SRM_PRISM_tumor_v1.2.cct.gz",
                               "UCEC_confirmatory_SRM_PRISM_tumor_v2.0.cct.gz",]:
                df_prism = pd.read_csv(file_path, sep='\t')
                df_prism[['Name','Peptide']] = df_prism['idx'].str.rsplit("-", 1, expand=True)
                df_prism = df_prism.set_index(["Name", "Peptide"])
                df_prism = df_prism.drop(columns=["idx"])
                df_prism = df_prism.transpose()
                df_prism = df_prism.sort_index()
                df_prism.index.name = "Patient_ID"
                
                # Merge if we have both
                if "targeted_proteomics" in self._data:
                    df_direct = self._data["targeted_proteomics"]
                    df_combined = pd.concat([df_direct, df_prism])
                    df_combined.index.name = "Patient_ID"
                    df_combined.columns.name = "Name"
                    self._data["targeted_proteomics"] = df_combined
                else:
                    self._data["targeted_proteomics"] = df_prism
                    
            elif file_name in ["UCEC_confirmatory_SRM_IMAC_tumor_v1.1.cct.gz",
                               "UCEC_confirmatory_SRM_IMAC_tumor_v1.2.cct.gz",
                               "UCEC_confirmatory_SRM_IMAC_tumor_v2.0.cct.gz",]:
                df = pd.read_csv(file_path, sep='\t')
                df.at[0,'idx'] = "FPSS[+80]PLRIPGGNIY[+80]ISPLK"
                df['Name'] = "RB1"
                df = df.rename(columns={"idx":"Peptide"})
                df = df.set_index(["Name", "Peptide"])
                df = df.transpose()
                df = df.sort_index()
                df.columns.name = "Name"
                df.index.name = "Patient_ID"
                self._data["targeted_phosphoproteomics"] = df
                
            elif file_name in ["UCEC_confirmatory_WES_cnv_gistic_thresholded_tumor_v1.0.cct.gz", 
                               "UCEC_confirmatory_WES_cnv_gistic_thresholded_tumor_v1.1.cct.gz",
                               "UCEC_confirmatory_WGS_cnv_gistic_thresholded_tumor_v1.2.cct.gz",
                               "UCEC_confirmatory_WGS_cnv_gistic_thresholded_tumor_v2.0.cct.gz",]:
                df = pd.read_csv(file_path, sep='\t')
                df[['Name','Chromosome']] = df.idx.str.split("|", expand=True)
                df = df.set_index(["Name"])
                df = df.drop(columns=["idx", "Chromosome"])
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                self._data["CNV_gistic"] = df
                
            elif file_name in ["UCEC_confirmatory_WES_cnv_log2_ratio_tumor_v1.0.cct.gz", 
                               "UCEC_confirmatory_WES_cnv_log2_ratio_tumor_v1.1.cct.gz",
                               "UCEC_confirmatory_WGS_cnv_log2_ratio_tumor_v1.2.cct.gz",
                               "UCEC_confirmatory_WGS_cnv_log2_ratio_tumor_v2.0.cct.gz",]:
                df = pd.read_csv(file_path, sep='\t')
                df[['Name','Chromosome']] = df.idx.str.split("|", expand=True)
                df = df.set_index(["Name"])
                df = df.drop(columns=["idx", "Chromosome"])
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                self._data["CNV_log2ratio"] = df
                
            elif file_name in ["UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.0.cbt.gz", 
                               "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.1.cbt.gz",
                               "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.2.cbt.gz"]:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["somatic_mutation_binary"] = df

            elif file_name in ["UCEC_confirmatory_WES_somatic_mutation_v1.0.maf.gz", 
                               "UCEC_confirmatory_WES_somatic_mutation_v1.1.maf.gz",
                               "UCEC_confirmatory_WES_somatic_mutation_v1.2.maf.gz",
                               "UCEC_confirmatory_WES_somatic_mutation_v2.0.maf.gz",]:
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
                self._data["somatic_mutation"] = df

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')


        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data)

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
        new_clinical = self._data["clinical"]
        new_clinical = new_clinical.reindex(master_index)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = new_clinical

        # Call function from dataframe_tools.py to sort all tables first by sample status, and then by the index
        self._data = sort_all_rows(self._data)

        # Call function from dataframe_tools.py to standardize the names of the index and column axes
        self._data = standardize_axes_names(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message

        # Print password access only warning
        warnings.warn("The UcecConf data is currently strictly reserved for CPTAC investigators. Otherwise, you are not authorized to access these data. Additionally, even after these data become publicly available, they will be subject to a publication embargo (see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter cptac.embargo() to open the webpage for more details).", PublicationEmbargoWarning, stacklevel=2)
        
    def get_CNV(self, algorithm=None):
        if (not algorithm):
            message = ("Please specify which type of UcecConf CNV data you want: "
            "'log2ratio' or 'gistic'. i.e. get_CNV('gistic')")
            raise InvalidParameterError(message)
        elif (algorithm == "log2ratio"):
            return super()._get_dataframe("CNV_log2ratio")
        elif (algorithm == "gistic"):
            return super()._get_dataframe("CNV_gistic")
        else: 
            message = ("Please specify a valid algorithm type for UcecConf CNV data: "
            "'log2ratio' or 'gistic'. i.e. get_CNV('gistic')")
            raise InvalidParameterError(message) 
