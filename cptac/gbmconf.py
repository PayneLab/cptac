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
from .exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError

class GbmConf(Dataset):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        valid_versions = ["0.1", "2.0"]

        data_files = {
            "0.1": [
                "acetylome_pnnl_d6.v0.1.20220202.tsv.gz",
                "clinical_data_core.v0.1.20220202.tsv.gz",
                "phosphoproteome_pnnl_d6.v0.1.20220202.tsv.gz",
                "proteome_pnnl_per_gene_d4.v0.1.20220202.tsv.gz",
                "rnaseq_washu_fpkm_uq.v0.1.20220202.tsv.gz",
                "somatic_wes_mutation.v0.1.20220202.maf.gz",
                "wgs_somatic_cnv_per_gene.v0.1.20220202.tsv.gz",
                ],
            "2.0": [
                "acetylome_pnnl_d6.v2.0.20220408.tsv.gz",
                "clinical_data_core.v2.0.20220408.tsv.gz",
                "Direct_SRM_corrected.v2.0.20220408.tsv.gz",
                "IMAC_SRM_corrected.v2.0.20220408.tsv.gz",
                "metabolome_pnnl.v2.0.20220408.tsv.gz",
                "mirnaseq_mirna_mature_tpm.v2.0.20220408.tsv.gz",
                "negative_lipidome_pnnl.v2.0.20220408.tsv.gz",
                "phosphoproteome_pnnl_d6.v2.0.20220408.tsv.gz",
                "positive_lipidome_pnnl.v2.0.20220408.tsv.gz",
                "PRISM_SRM_raw.v2.0.20220408.tsv.gz",
                "proteome_pnnl_per_gene_d4.v2.0.20220408.tsv.gz",
                "rnaseq_gene_fusion.v2.0.20220408.tsv.gz",
                "rnaseq_washu_fpkm_uq.v2.0.20220408.tsv.gz",
                "somatic_wes_mutation.v2.0.20220408.maf.gz",
                "wgs_somatic_cnv_gatk4_per_gene.v2.0.20220408.tsv.gz",
                ],
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="gbmconf", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file.
            df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name, so we don't include the version

            # Instructions for parsing each cancer type
            if df_name == "acetylome_pnnl_d6":
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

                self._data["acetylproteomics"] = df
                
            elif df_name == "clinical_data_core":
                df = pd.read_csv(file_path, sep='\t')
                
                # Save all sample ids that belong to the discovery cohort so they can be removed from the other dataframes
                dfd = df.loc[df['cohort'] == 'Discovery']
                dfd = dfd["preferred_sample_name"].to_frame()
                self._data["discovery_cohort"] = dfd
                
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
                
                self._data["clinical"] = df
                
            elif df_name == "Direct_SRM_corrected":
                df = pd.read_csv(file_path, sep='\t')
                df = df.rename(columns={
                    "protein_name": "Name",
                    "peptide_seq": "Peptide",
                })
                df = df.set_index(["Name", "Peptide"])
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                
                self._data["direct_SRM"] = df
                
            # This is out of alphabetical order because it belongs with the other SRM things
            elif df_name == "PRISM_SRM_raw":
                df = pd.read_csv(file_path, sep='\t')
                df = df.rename(columns={
                    "protein_name": "Name",
                    "peptide_seq": "Peptide",
                })
                df = df.set_index(["Name", "Peptide"])
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"

                self._data["prism_SRM"] = df
                
            elif df_name == "IMAC_SRM_corrected":
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
                df
                
                self._data["targeted_phosphoproteomics"] = df
                
            elif df_name == "metabolome_pnnl":
                df = pd.read_csv(file_path, sep='\t')
                df = df.set_index("Metabolite")
                df = df.transpose()
                # Name columns for consistency
                df.columns.name = "Name"
                df.index.name="Patient_ID"
                
                self._data["metabolomics"] = df
                
            elif df_name == "mirnaseq_mirna_mature_tpm":
                df = pd.read_csv(file_path, sep='\t')
                df = df.drop(columns=["Alias", "Derives_from"])
                df = df.rename(columns={'ID': "Database_ID"})
                df = df.set_index(["Name", "Database_ID"])
                df = df.transpose()
                df.index.name="Patient_ID"
                
                self._data["miRNA"] = df
                
            elif df_name == "negative_lipidome_pnnl":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df.columns.name = "Name"
                df.index.name="Patient_ID"
                df = df.add_suffix("_negative")
                
                self._data["lipidomics_negative"] = df
                
            elif df_name == "phosphoproteome_pnnl_d6":
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
                
                self._data["phosphoproteomics"] = df
                
            elif df_name == "positive_lipidome_pnnl":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df.columns.name = "Name"
                df.index.name="Patient_ID"
                df = df.add_suffix("_positive")
                
                self._data["lipidomics_positive"] = df
                
            elif df_name == "proteome_pnnl_per_gene_d4":
                df = pd.read_csv(file_path, sep='\t')
                # Save refseq information for use in targeted_proteomics
                df_ref = df[["gene", "refseq_id"]]
                df_ref = df_ref.rename(columns={
                    "gene": "Name",
                    "refseq_id": "Database_ID",
                })
                df_ref = df_ref.set_index("Name")
                self._data["prot_refseq"] = df_ref
                
                # Rename columns and multiindex
                df = df.rename(columns={"gene": "Name", 'refseq_id': "Database_ID"})
                df = df.set_index(["Name", "Database_ID"])
                df = df.transpose()
                df.index.name="Patient_ID"
                
                self._data["proteomics"] = df
                
            elif df_name == "rnaseq_gene_fusion":
                df = pd.read_csv(file_path, sep='\t')
                df = df.set_index("preferred_sample_name")
                df.columns.name = "Name"
                df.index.name = "Patient_ID"
                
                self._data["gene_fusion"] = df
                
            elif df_name == "rnaseq_washu_fpkm_uq":
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
                
                self._data["transcriptomics"] = df
            
            elif df_name == "somatic_wes_mutation":
                df = pd.read_csv(file_path, sep='\t')

                if self._version == "0.1":
                    # We don't need any of the other columns
                    df = df[["preferred_sample_name", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short"]]
                    df = df.rename(columns={
                        "preferred_sample_name": "Patient_ID",
                        'Hugo_Symbol': "Gene",
                        "Variant_Classification": "Mutation",
                        "HGVSp_Short": "Location"
                    })
                elif self._version == "2.0":
                    df = df[["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short"]]
                    df = df.rename(columns={"Tumor_Sample_Barcode": "Patient_ID", 'Hugo_Symbol': "Gene", "Variant_Classification": "Mutation", "HGVSp_Short": "Location"})
                    df["Patient_ID"] = df["Patient_ID"].str.rstrip("_T")

                df = df.set_index("Patient_ID")
                df = df.sort_values(by=["Patient_ID","Gene"])

                self._data["somatic_mutation"] = df

            # This is the CNV file used for 0.1
            elif df_name == "wgs_somatic_cnv_per_gene":
                df = pd.read_csv(file_path, sep='\t')
                df = df.drop(columns=["Start","End", "Chr"])
                df = df.set_index("Gene")
                df = df.sort_index()
                df = df.transpose()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"

                self._data["CNV"] = df

            # This is the CNV file used for 2.0
            elif df_name == "wgs_somatic_cnv_gatk4_per_gene":
                df = pd.read_csv(file_path, sep='\t')

                df = df.set_index("Gene")
                df = df.sort_index()
                df = df.transpose()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"

                self._data["CNV"] = df

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')
        
        if self._version in ("2.0"):
            # Combine positive and negative lipidomics tables
            lipidomics_positive = self._data["lipidomics_positive"]
            lipidomics_negative = self._data["lipidomics_negative"]

            lipidomics = lipidomics_positive.join(lipidomics_negative, how="outer")
            self._data["lipidomics"] = lipidomics

            del self._data["lipidomics_positive"]
            del self._data["lipidomics_negative"]
            
            # Combine Direct and PRISM SRM tables
            # These do not contain refseq ids
            df_direct = self._data["direct_SRM"]
            df_prism = self._data["prism_SRM"]
            
            targeted_proteomics = df_prism.join(df_direct, how="outer")
            # Add refseq_ids (saved previously in prot_refseq)
            targeted_proteomics = targeted_proteomics.transpose()
            targeted_proteomics = targeted_proteomics.join(self._data["prot_refseq"], how="left")
            targeted_proteomics = targeted_proteomics.reset_index()
            levels = ["Name", "Peptide", "Database_ID"]
            targeted_proteomics = targeted_proteomics.set_index(levels)
            # This next line is basically transpose but it doesn't mess up the multiindex
            targeted_proteomics = targeted_proteomics.stack().unstack(levels)
            self._data["targeted_proteomics"] = targeted_proteomics
            
            del self._data["direct_SRM"]
            del self._data["prism_SRM"]
            
            # Delete discovery cohort samples from all data types
            # This will probably change significantly to implement lazy loading
            self._data["discovery_cohort"] = dfd

            # acetylomics apparently has no discovery cohort ids
            # The rest have not been checked, gene_fusion for sure does, but it won't hurt a datatype to go through this process
            # gene fusion
            temp = self._data["gene_fusion"]
            temp = temp[~temp.index.isin(dfd["preferred_sample_name"])]
            self._data["gene_fusion"] = temp
            # targeted_proteomics
            temp = self._data["targeted_proteomics"]
            temp = temp[~temp.index.isin(dfd["preferred_sample_name"])]
            self._data["targeted_proteomics"] = temp
            # targeted_phosphoproteomics
            temp = self._data["targeted_phosphoproteomics"]
            temp = temp[~temp.index.isin(dfd["preferred_sample_name"])]
            self._data["targeted_phosphoproteomics"] = temp
            # metabolomics
            temp = self._data["metabolomics"]
            temp = temp[~temp.index.isin(dfd["preferred_sample_name"])]
            self._data["metabolomics"] = temp
            # miRNA
            temp = self._data["miRNA"]
            temp = temp[~temp.index.isin(dfd["preferred_sample_name"])]
            self._data["miRNA"] = temp
            # lipidomics
            temp = self._data["lipidomics"]
            temp = temp[~temp.index.isin(dfd["preferred_sample_name"])]
            self._data["lipidomics"] = temp
            # phosphoproteomics
            temp = self._data["phosphoproteomics"]
            temp = temp[~temp.index.isin(dfd["preferred_sample_name"])]
            self._data["phosphoproteomics"] = temp
            # proteomics
            temp = self._data["proteomics"]
            temp = temp[~temp.index.isin(dfd["preferred_sample_name"])]
            self._data["proteomics"] = temp
            # transcriptomics
            temp = self._data["transcriptomics"]
            temp = temp[~temp.index.isin(dfd["preferred_sample_name"])]
            self._data["transcriptomics"] = temp
            # somatic_mutation
            temp = self._data["somatic_mutation"]
            temp = temp[~temp.index.isin(dfd["preferred_sample_name"])]
            self._data["somatic_mutation"] = temp
            # CNV
            temp = self._data["CNV"]
            temp = temp[~temp.index.isin(dfd["preferred_sample_name"])]
            self._data["CNV"] = temp
            
        del self._data["discovery_cohort"]    
        del self._data["prot_refseq"]

        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data) 

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. 
        # Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
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
        warnings.warn("The GBM confirmatory data is currently strictly reserved for CPTAC investigators. Otherwise, you are not authorized to access these data. Additionally, even after these data become publicly available, they will be subject to a publication embargo (see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter cptac.embargo() to open the webpage for more details).", PublicationEmbargoWarning, stacklevel=2)
