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
from .dataset import DataSet
from .file_download import update_index
from .file_tools import validate_version, get_version_files_paths
from .dataframe_tools import *
from .exceptions import NoInternetError, FailedReindexWarning

class Gbm(DataSet):

    def __init__(self, version="latest"):
        """Load all of the gbm dataframes as values in the self._data dict variable, with names as keys, and format them properly."""

        # Call the parent DataSet __init__ function, which initializes self._data and other variables we need
        super().__init__("gbm")

        # Update the index, if possible. If there's no internet, that's fine.
        try:
            update_index(self._cancer_type)
        except NoInternetError:
            pass

        # Validate the version
        self._version = validate_version(version, self._cancer_type, use_context="init")

        # Get the paths to all the data files
        data_files = [
            "clinical_data_core.v1.0.20190802.tsv.gz",
            "mirnaseq_mirna_mature_tpm.v1.0.20190802.tsv.gz",
            "phosphoproteome_pnnl_d6.v1.0.20190802.tsv.gz",
            "proteome_pnnl_per_gene_d4.v1.0.20190802.tsv.gz",
            "proteome_tmt_design.v1.0.20190802.tsv.gz",
            "rnaseq_gdc_fpkm_uq.v1.0.20190802.tsv.gz",
            "tindaisy_all_cases_filtered.v1.0.20190802.maf.gz",
            "wgs_somatic_cnv_per_gene.v1.0.20190802.tsv.gz",
        ]
        data_files_paths = get_version_files_paths(self._cancer_type, self._version, data_files)

        # Load the data into dataframes in the self._data dict
        loading_msg = "Loading dataframes"
        for file_path in data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            if file_name == "clinical_data_core.v1.0.20190802.tsv.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                self._data["clinical"] = df

            if file_name == "mirnaseq_mirna_mature_tpm.v1.0.20190802.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')
                df = df.rename(columns={"name": "Name", "unique_id": "Database_ID"})
                df = df.set_index(["Name", "Database_ID"]) # We use a multiindex with database IDs, not just names, to avoid duplicate column headers
                df = df.drop(columns=["chromosome", "start", "end", "strand", "mirna_type", "mirbase_id", "precursor_id"])
                df = df.sort_index()
                df = df.transpose()
                self._data["miRNA"] = df

            if file_name == "phosphoproteome_pnnl_d6.v1.0.20190802.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')

                # Create our multiindex
                split_genes = df["site"].str.rsplit("-", n=1, expand=True) # Split the genes from the sites, splitting from the right since some genes have hyphens in their names, but the genes and sites are also separated by hyphens
                df = df.drop(columns="site")
                df = df.assign(Site=split_genes[1])
                df["Site"] = df["Site"].str.replace(r"[sty]", r"") # Get rid of all lowercase s, t, and y delimeters in the sites
                df = df.rename(columns={"gene": "Name", "peptide": "Peptide"})
                df = df.set_index(["Name", "Site", "Peptide"]) # Turn these columns into a multiindex

                df = df.sort_index()
                df = df.transpose()
                self._data["phosphoproteomics"] = df

            if file_name == "proteome_pnnl_per_gene_d4.v1.0.20190802.tsv.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.sort_index()
                df = df.transpose()
                self._data["proteomics"] = df

            if file_name == "proteome_tmt_design.v1.0.20190802.tsv.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df.index.name = "Patient_ID"
                self._data["experimental_design"] = df

            if file_name == "rnaseq_gdc_fpkm_uq.v1.0.20190802.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')
                df = df.rename(columns={"gene_name": "Name", "gene_id": "Database_ID"})
                df = df.set_index(["Name", "Database_ID"]) # We use a multiindex with Ensembl IDs, not just gene names, to avoid duplicate column headers
                df = df.drop(columns=["gene_type", "gene_status", "havana_gene", "full_length", "exon_length", "exon_num"])
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                self._data["transcriptomics"] = df

            if file_name == "tindaisy_all_cases_filtered.v1.0.20190802.maf.gz":
                df = pd.read_csv(file_path, sep='\t')
                split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n = 1, expand = True) # The first part of the barcode is the patient id, which we need want to make the index
                df["Tumor_Sample_Barcode"] = split_barcode[0]
                df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
                df = df.rename({"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
                df = df.sort_values(by=["Patient_ID", "Gene"])
                df = df.set_index("Patient_ID")
                self._data["somatic_mutation"] = df

            if file_name == "wgs_somatic_cnv_per_gene.v1.0.20190802.tsv.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.drop(columns=["gene_id", "gene_id_version", "original_symbol"])
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                self._data["CNV"] = df


        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data)

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
        clinical = self._data["clinical"]
        clinical = clinical.reindex(master_index)

        # Copy the sample status column from the experimental_design dataframe to the clinical dataframe
        sample_status_col = self._data["experimental_design"]["tumor_normal"].copy()
        sample_status_col = sample_status_col.str.title()
        clinical.insert(0, "Sample_Tumor_Normal", sample_status_col)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = clinical

        # Generate a sample ID for each patient ID
        sample_id_dict = generate_sample_id_map(master_index)

        # Give all the dataframes Sample_ID indices
        dfs_to_delete = []
        for name in self._data.keys(): # Only loop over keys, to avoid changing the structure of the object we're looping over
            df = self._data[name]
            df.index.name = "Patient_ID"
            keep_old = name in ["clinical", "experimental_design"] # Keep the old Patient_ID index as a column in the clinical and experimental_design dataframes, so we have a record of it.
            try:
                df = reindex_dataframe(df, sample_id_dict, "Sample_ID", keep_old)
            except ReindexMapError:
                warnings.warn(f"Error mapping sample ids in {name} dataframe. At least one Patient_ID did not have corresponding Sample_ID mapped in clinical dataframe. {name} dataframe not loaded.", FailedReindexWarning, stacklevel=2) # stacklevel=2 ensures that the warning is registered as originating from the file that called this __init__ function, instead of from here directly, because the former is more useful information.
                dfs_to_delete.append(name)
                continue

            self._data[name] = df

        for name in dfs_to_delete: # Delete any dataframes that had issues reindexing
            del self._data[name]

        # Set name of column axis to "Name" for all dataframes
        for name in self._data.keys():
            df = self._data[name]
            df.columns.name = "Name"
            self._data[name] = df

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
