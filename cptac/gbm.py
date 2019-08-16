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
                df = pd.read_csv(file_path, sep='\t')
                self._data[""] = df

            if file_name == "mirnaseq_mirna_mature_tpm.v1.0.20190802.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')
                self._data[""] = df

            if file_name == "phosphoproteome_pnnl_d6.v1.0.20190802.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')
                self._data[""] = df

            if file_name == "proteome_pnnl_per_gene_d4.v1.0.20190802.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')
                self._data[""] = df

            if file_name == "proteome_tmt_design.v1.0.20190802.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')
                self._data[""] = df

            if file_name == "rnaseq_gdc_fpkm_uq.v1.0.20190802.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')
                duplicates = df["gene_name"].duplicated(keep=False) # Find which gene names have duplicates
                needed_ids = df["gene_id"].where(duplicates, other="") # Get the unique Ensembl IDs for duplcated genes
                needed_ids = needed_ids.where(~duplicates, other=('|' + needed_ids)) # Prepend a "|" to the rows we're going to append Ensembl IDs to, to separate the gene name and the ID
                df["gene_name"] = df["gene_name"].str.cat(needed_ids) # Append the Ensembl IDs. Ta da!! No duplicates!!
                df = df.drop(columns=["gene_id", "gene_type", "gene_status", "havana_gene", "full_length", "exon_length", "exon_num"])
                df = df.set_index("gene_name")
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

        # FILL: Here, write code to format your dataframes properly. Requirements:
        # - All dataframes must be indexed by Sample_ID, not Patient_ID.
        #     - This means that two samples from the same patient will have the same Patient_ID, but different Sample_ID numbers.
        #     - Sample_ID numbers must be of the format S***, e.g. S001, S028, S144
        #     - clinical dataframe must contain a Patient_ID column that contains the Patient_ID for each sample
        #     - If the data did not come indexed with Sample_ID numbers, look at the Ovarian dataset for an example of generating Sample_ID numbers and mapping them to Patient_ID numbers.
        #     - Note that most datasets are originally indexed with a patient id, which we rename as the case id, and has a format like C3N-00352.
        #         - If one patient provided both a normal and a tumor sample, those samples will have the same patient/case id. Therefore, before any joining or reindexing, prepend an 'N' to all normal sample case ids, based on the column in the clinical dataframe indicating which samples are tumor or normal. See existing datasets for examples of how to do this.
        #
        # - Each dataframe's name must match the format for that type of dataframe in all the other datasets. 
        #     - E.g., if your binary mutations dataframe is named mutations_binary, you'd need to rename it to somatic_mutation_binary to match the other datasets' binary mutation dataframes.
        #
        # - If the new dataset has a dataframe not included in any other datasets, you must write a getter for it in the parent DataSet class, found in cptac/dataset.py, using the private method DataSet._get_dataframe
        #
        # - You'd also need to add the new dataframe's name to self._valid_omics_dfs if it's a valid omics df for the DataSet merge functions, or self._valid_metadata_dfs if it's a valid metadata df for DataSet.append_metadata_to_omics
        #     - Note that a dataframe with multiple rows for each sample, like the treatment dataframe in the Ovarian dataset, should not be a valid dataset for joining
        #
        # - If any dataframes are split between two files--such as one file for the tumor sample proteomics, and one file for the normal sample proteomics--they'll have been read into separate dataframes, and you need to merge those into one dataframe.
        #     - Make sure that samples coming from a normal file have an 'N' added to their Patient_ID numbers, to keep a record of which ones are normal samples.
        #
        # - If multiple dataframes are contained in one file--e.g. clinical and derived_molecular data are both in clinical.txt, as in Endometrial--separate them out here.
        #
        # - Make sure that column names are consistent--e.g., all Patient_ID columns should be labeled as such, not as Clinical_Patient_Key or Case_ID or something else. Rename columns as necessary to match this.
        #
        # - The clinical dataframe must contain a Sample_Tumor_Normal column, which contains either "Tumor" or "Normal" for each sample, according to its status.
        #
        # - Only the clinical dataframe should contain a Patient_ID column. The other dataframes should contain just a Sample_ID index, and the data.
        #
        # - The column axis of each dataframe should have None as the value of its .name attribute
        #
        # - The index of each dataframe should have "Sample_ID" as the value of its .name attribute, since that's what the index is.
        #
        # - Make sure to drop any excluded cases, as in Endometrial.
        #
        # - Make sure that in dataframes where each column header is the name of a gene, the columns are in alphabetical order.
        #
        # - Dataframes that require parsing sites: phosphoproteomics, acetylproteomics.
        #
        # - General workflow when you need to reindex a dataset to have Sample_IDs (S*** numbers). See cptac/renalccrcc.py for an example.
        #     1. Get all the dataframes in the dataset to have the same type of index (patient id, case id, etc.)
        #     2. Use dataframe_tools.unionize_indices to get a union of all indices in the dataset (because some dataframes might have samples not included in other dataframes, and we want to include all samples when we create our new index).
        #     3. Pass the master index to the DataFrame.reindex function, called on the clinical dataframe, to make sure all indices are stored in the clinical dataframe, e.g. clinical = clinical.reindex(master_index)
        #     4. Use dataframe_tools.generate_sample_id_map to generate a map from the current index to Sample_ID numbers
        #     5. Loop through all dataframes in the dataset, and use dataframe_tools.reindex_dataframe, passing the sample ID map, to reindex all dataframes with Sample IDs. Only keep the old index in the clinical dataframe, so we have a record of it.
        #     6. Make sure df.columns.name is None for all dataframes, and df.index.name is Sample_ID for all dataframes
        #
        # Below is an example of how to do this. You may need to tweak it for your dataset.


#*****************************#FILL: EXAMPLE CODE--MAKE SURE TO FIX SO IT ACTUALLY WORKS FOR YOUR DATASET!!!!******************************

        # Get a union of all dataframes' indices, with duplicates removed
#        master_index = unionize_indices(self._data)
#
#        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
#        clinical = self._data["clinical"]
#        clinical = clinical.reindex(master_index)
#
#        # Replace the clinical dataframe in the data dictionary with our new and improved version!
#        self._data['clinical'] = clinical
#
#        # Generate a sample ID for each patient ID
#        sample_id_dict = generate_sample_id_map(master_index)
#
#        # Give all the dataframes Sample_ID indices
#        dfs_to_delete = []
#        for name in self._data.keys(): # Only loop over keys, to avoid changing the structure of the object we're looping over
#            df = self._data[name]
#            df.index.name = "Patient_ID"
#            keep_old = (name == "clinical") # Keep the old Patient_ID index as a column in the clinical dataframe, so we have a record of it.
#            try:
#                df = reindex_dataframe(df, sample_id_dict, "Sample_ID", keep_old)
#            except ReindexMapError:
#                warnings.warn(f"Error mapping sample ids in {name} dataframe. At least one Patient_ID did not have corresponding Sample_ID mapped in clinical dataframe. {name} dataframe not loaded.", FailedReindexWarning, stacklevel=2) # stacklevel=2 ensures that the warning is registered as originating from the file that called this __init__ function, instead of from here directly, because the former is more useful information.
#                dfs_to_delete.append(name)
#                continue
#
#            self._data[name] = df
#
#        for name in dfs_to_delete: # Delete any dataframes that had issues reindexing
#            del self._data[name]
#
        # Drop name of column axis for all dataframes
        for name in self._data.keys():
            df = self._data[name]
            df.columns.name = None
            self._data[name] = df


#*****************#FILL: END OF EXAMPLE*********************
            

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message

        # FILL: remove the below message printing, if the dataset isn't under publiction embargo
        # Print data embargo warning, if the date hasn't passed yet.
#        today = datetime.date.today()
#        embargo_date = datetime.date() # FILL: Insert embargo date here.
#        if today < embargo_date:
#            warning = "WARNING: This data is under a publication embargo until # FILL: INSERT EMBARGO DATE HERE#. CPTAC is a community resource project and data are made available rapidly after generation for community research use. The embargo allows exploring and utilizing the data, but analysis may not be published until July 1, 2019. Please see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter cptac.embargo() to open the webpage for more details."
#            wrapped_list = textwrap.wrap(warning)
#            for line in wrapped_list:
#                print(line)
