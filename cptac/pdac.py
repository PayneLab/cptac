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

class Pdac(Dataset):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        valid_versions = ["1.0"]

        data_files = {
            "1.0": [
            "clinical_table_140.tsv.gz",
            "microRNA_TPM_log2_Normal.cct.gz",
            "microRNA_TPM_log2_Tumor.cct.gz",
            "meta_table_140.tsv.gz",
            "mRNA_RSEM_UQ_log2_Normal.cct.gz",
            "mRNA_RSEM_UQ_log2_Tumor.cct.gz",
            "PDAC_mutation.maf.gz",
            "phosphoproteomics_site_level_MD_abundance_normal.cct.gz",
            "phosphoproteomics_site_level_MD_abundance_tumor.cct.gz",
            "proteomics_gene_level_MD_abundance_normal.cct.gz",
            "proteomics_gene_level_MD_abundance_tumor.cct.gz",
            "RNA_fusion_unfiltered_normal.tsv.gz",
            "RNA_fusion_unfiltered_tumor.tsv.gz",
            "SCNA_log2_gene_level.cct.gz"],
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="pdac", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below
            mark_normal = lambda s: s + ".N"
            remove_type_tag = lambda s: s[:-2] # remove _T and similar tags from end of string

            if file_name == "clinical_table_140.tsv.gz": # Note that we use the "file_name" variable to identify files. That way we don't have to use the whole path.
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.rename_axis("Patient_ID", axis="index")
                df = df.sort_index()
                df.columns.name = "Name"
                df["Sample_Tumor_Normal"] = "Tumor"
                self._data["clinical"] = df

            elif file_name == "meta_table_140.tsv.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.sort_index()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["derived_molecular"] = df

            elif file_name == "microRNA_TPM_log2_Normal.cct.gz":
                df_normal = pd.read_csv(file_path, sep='\t', index_col=0)
                df_normal = df_normal.sort_index()
                df_normal = df_normal.transpose()
                df_normal = df_normal.rename(index=mark_normal)

                # merge tumor and normal if tumor data has already been read
                if "miRNA" in self._data:
                    df_tumor = self._data["miRNA"]
                    df_combined = pd.concat([df_normal, df_tumor])
                    df_combined.index.name = "Patient_ID"
                    df_combined.columns.name = "Name"
                    self._data["miRNA"] = df_combined
                else:
                    self._data["miRNA"] = df_normal

            elif file_name == "microRNA_TPM_log2_Tumor.cct.gz":
                df_tumor = pd.read_csv(file_path, sep='\t', index_col=0)
                df_tumor = df_tumor.sort_index()
                df_tumor = df_tumor.transpose()

                # merge tumor and normal if normal data has already been read
                if "miRNA" in self._data:
                    df_normal = self._data["miRNA"]
                    df_combined = pd.concat([df_normal, df_tumor])
                    df_combined.index.name = "Patient_ID"
                    df_combined.columns.name = "Name"
                    self._data["miRNA"] = df_combined
                else:
                    self._data["miRNA"] = df_tumor

            elif file_name == "mRNA_RSEM_UQ_log2_Normal.cct.gz":
                # create df for normal data
                df_normal = pd.read_csv(file_path, sep='\t', index_col=0)
                df_normal = df_normal.sort_index()
                df_normal = df_normal.transpose()
                df_normal = df_normal.rename(index=mark_normal)
                
                # merge tumor and normal if tumor data has already been read
                if "transcriptomics" in self._data:
                    df_tumor = self._data["transcriptomics"]
                    df_combined = pd.concat([df_normal, df_tumor])
                    df_combined.index.name = "Patient_ID"
                    df_combined.columns.name = "Name"
                    self._data["transcriptomics"] = df_combined
                else:
                    self._data["transcriptomics"] = df_normal

            elif file_name == "mRNA_RSEM_UQ_log2_Tumor.cct.gz":
                # create df for tumor data
                df_tumor = pd.read_csv(file_path, sep='\t', index_col=0)
                df_tumor = df_tumor.sort_index()
                df_tumor = df_tumor.transpose()

                # merge tumor and normal if normal data has already been read
                if "transcriptomics" in self._data:
                    df_normal = self._data["transcriptomics"]
                    df_combined = pd.concat([df_normal, df_tumor])
                    df_combined.index.name = "Patient_ID"
                    df_combined.columns.name = "Name"
                    self._data["transcriptomics"] = df_combined
                else:
                    self._data["transcriptomics"] = df_tumor

            elif file_name == "PDAC_mutation.maf.gz":
                df = pd.read_csv(file_path, sep='\t')
                df = df[["Hugo_Symbol", "Variant_Classification", "HGVSp_Short", "Tumor_Sample_Barcode"]]
                df = df.rename({"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
                df = df.sort_values(by=["Patient_ID", "Gene"])
                df = df.set_index("Patient_ID")
                df = df.rename(index=remove_type_tag)
                df.columns.name = "Name"
                self._data["somatic_mutation"] = df

            elif file_name == "phosphoproteomics_site_level_MD_abundance_normal.cct.gz":
                # create df form normal data
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
                df_normal = df_normal.rename(index=mark_normal)

                # merge tumor and normal if tumor data has already been read
                if "phosphoproteomics" in self._data:
                    df_tumor = self._data["phosphoproteomics"]
                    df_combined = pd.concat([df_normal, df_tumor])
                    df_combined.index.name = "Patient_ID"
                    #df_combined.columns.name = "Name"
                    self._data["phosphoproteomics"] = df_combined
                else:
                    self._data["phosphoproteomics"] = df_normal

            elif file_name == "phosphoproteomics_site_level_MD_abundance_tumor.cct.gz":
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
                
                # merge tumor and normal if normal data has already been read
                if "phosphoproteomics" in self._data:
                    df_normal = self._data["phosphoproteomics"]
                    df_combined = pd.concat([df_normal, df_tumor])
                    df_combined.index.name = "Patient_ID"
                    #df_combined.columns.name = "Name"
                    self._data["phosphoproteomics"] = df_combined
                else:
                    self._data["phosphoproteomics"] = df_tumor
        
            elif file_name == "proteomics_gene_level_MD_abundance_normal.cct.gz":
                df_normal = pd.read_csv(file_path, sep='\t', index_col=0)
                df_normal = df_normal.sort_index()
                df_normal = df_normal.transpose()
                df_normal = df_normal.rename(index=mark_normal)

                # merge tumor and normal if tumor data has already been read
                if "proteomics" in self._data:
                    df_tumor = self._data["proteomics"]
                    df_combined = pd.concat([df_normal, df_tumor])
                    df_combined.index.name = "Patient_ID"
                    df_combined.columns.name = "Name"
                    self._data["proteomics"] = df_combined
                else:
                    self._data["proteomics"] = df_normal

            elif file_name == "proteomics_gene_level_MD_abundance_tumor.cct.gz":
                df_tumor = pd.read_csv(file_path, sep='\t', index_col=0)
                df_tumor = df_tumor.sort_index()
                df_tumor = df_tumor.transpose()

                # merge tumor and normal if normal data has already been read
                if "proteomics" in self._data:
                    df_normal = self._data["proteomics"]
                    df_combined = pd.concat([df_normal, df_tumor])
                    df_combined.index.name = "Patient_ID"
                    df_combined.columns.name = "Name"
                    self._data["proteomics"] = df_combined
                else:
                    self._data["proteomics"] = df_tumor

            elif file_name == "RNA_fusion_unfiltered_normal.tsv.gz":
                df_normal = pd.read_csv(file_path, sep='\t', index_col=0)
                df_normal = df_normal.rename(columns={"Sample": "Patient_ID"})
                df_normal = df_normal.set_index("Patient_ID")
                df_normal = df_normal.rename(index=mark_normal)

                if "gene_fusion" in self._data:
                    df_tumor = self._data ["gene_fusion"]
                    df_combined = pd.concat([df_normal, df_tumor])
                    df_combined.index.name = "Patient_ID"
                    df_combined.columns.name = "Name"
                    self._data["gene_fusion"] = df_combined
                else:
                    self._data["gene_fusion"] = df_normal

            elif file_name == "RNA_fusion_unfiltered_tumor.tsv.gz":
                df_tumor = pd.read_csv(file_path, sep='\t', index_col=0)
                df_tumor = df_tumor.rename(columns={"Sample": "Patient_ID"})
                df_tumor = df_tumor.set_index("Patient_ID")

                if "gene_fusion" in self._data:
                    df_normal = self._data ["gene_fusion"]
                    df_combined = pd.concat([df_normal, df_tumor])
                    df_combined.index.name = "Patient_ID"
                    df_combined.columns.name = "Name"
                    self._data["gene_fusion"] = df_combined
                else:
                    self._data["gene_fusion"] = df_tumor

            elif file_name == "SCNA_log2_gene_level.cct.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df.columns.name = "Name"
                df.index.name = "Patient_ID"
                self._data["CNV"] = df


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
        warnings.warn("The pdac data is currently strictly reserved for CPTAC investigators. "
            "Otherwise, you are not authorized to access these data. "
            "Additionally, even after these data become publicly available, "
            "they will be subject to a publication embargo "
            "(see https://proteomics.cancer.gov/data-portal/about/data-use-agreement "
            "or enter cptac.embargo() to open the webpage for more details).", 
            PublicationEmbargoWarning, stacklevel=2)
    def how_to_cite(self):
        return super().how_to_cite(cancer_type='pancreatic ductal adenocarcinoma', pmid='', unpublished=True)