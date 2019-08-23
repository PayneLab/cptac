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
import os
import warnings
from .dataset import DataSet
from .file_download import update_index
from .file_tools import validate_version, get_version_files_paths
from .dataframe_tools import *
from .exceptions import FailedReindexWarning, NoInternetError, ReindexMapError

class Ovarian(DataSet):

    def __init__(self, version="latest"):
        """Load all the ovarian dataframes as values in the self._data dict variable, with names as keys, and format them properly."""

        # Call the parent Dataset __init__() function, which initializes self._data and other variables we need
        super().__init__("ovarian")

        # Update the index, if possible. If there's no internet, that's fine.
        try:
            update_index(self._cancer_type)
        except NoInternetError:
            pass

        # Validate the index
        self._version = validate_version(version, self._cancer_type, use_context="init")

        # Get the paths to all the data files
        data_files = [
            "clinical.csv.gz",
            "cnv.tsv.gz",
            "definitions.txt",
            "phosphoproteomics.txt.gz",
            "proteomics.txt.gz",
            "somatic_38.maf.gz",
            "transcriptomics.tsv.gz",
            "treatment.csv.gz"]
        data_files_paths = get_version_files_paths(self._cancer_type, self._version, data_files)

        # Load the data files into dataframes in the self._data dict
        loading_msg = "Loading dataframes"
        for file_path in data_files_paths:

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of all the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            if file_name == "clinical.csv.gz" or file_name == "treatment.csv.gz":
                df = pd.read_csv(file_path, sep=",", index_col=0)
                df = df.rename(columns={"Participant_ID":"Patient_ID"})
                df = df.set_index("Patient_ID")
                self._data[df_name] = df #maps dataframe name to dataframe

            elif file_name == "cnv.tsv.gz":
                df = pd.read_csv(file_path, sep="\t", index_col=0)
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                self._data["CNV"] = df #maps dataframe name to dataframe

            elif file_name == "definitions.txt":
                with open(file_path, "r", errors="ignore") as definitions_file:
                    for line in definitions_file.readlines():
                        line = line.strip()
                        line = line.split("\t")
                        term = line[0]
                        definition = line[1]
                        self._definitions[term] = definition

            elif file_name == "phosphoproteomics.txt.gz" or file_name == "proteomics.txt.gz":
                df = pd.read_csv(file_path, sep='\t')
                   
                if file_name == "proteomics.txt.gz":
                    df = df[df["hgnc_symbol"].notnull()] # Drops all nan values in hgnc_symbol column

                    # Create our column multiindex
                    df = df.rename(columns={"hgnc_symbol": "Name", "refseq_peptide": "Database_ID"})
                    df = df.set_index(["Name", "Database_ID"])

                elif file_name == "phosphoproteomics.txt.gz":
                    df = df[df["site"].notnull()] # Drops all rows with nan values in site column

                    # Create our column multiindex
                    split_genes = df["site"].str.rsplit("-", n=1, expand=True) # Split the genes from the sites, splitting from the right since some genes have hyphens in their names, but the genes and sites are also separated by hyphens
                    df = df.drop(columns=["hgnc_symbol", "site"]) # hgnc_symbol is a duplicate of split_genes[0], and site is now in split_genes and will be re-inserted differently
                    df = df.assign(Name=split_genes[0], Site=split_genes[1])
                    df["Site"] = df["Site"].str.replace(r"[sty]", r"") # Get rid of all lowercase s, t, and y delimeters in the sites
                    df = df.rename(columns={"refseq_peptide": "Database_ID"})
                    df = df.set_index(["Name", "Site", "Peptide", "Database_ID"]) # Turn these columns into a multiindex

                df = df.sort_index()
                df = df.transpose()
                df.index = df.index.where(~df.index.str.startswith('C'), df.index.str[1:]) # Take C prefix off of indices for those samples that have them (tumor samples have C, normal have N)
                df = df.drop(index=df.index[df.index.str.startswith("OV_QC")]) # Drop all OV_QC samples--they're quality control samples not relevant for data analysis
                self._data[df_name] = df

            elif file_name == "somatic_38.maf.gz":
                df = pd.read_csv(file_path, sep = "\t", index_col=0)
                df = df.reset_index()
                split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n = 1, expand = True) # The first part of the barcode is the patient id, which we need to make a Patient_ID column
                df["Tumor_Sample_Barcode"] = split_barcode[0]
                df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]] # We only want these columns
                df = df.rename(columns={"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"})
                df = df.sort_values(by=["Patient_ID", "Gene"])
                df = df.set_index("Patient_ID")
                self._data['somatic_mutation'] = df

            elif file_name == "transcriptomics.tsv.gz":
                df = pd.read_csv(file_path, sep="\t", index_col=0)
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                date_cols = ['1-Dec', '1-Sep', '10-Mar', '10-Sep', '11-Sep', '12-Sep', '14-Sep', '15-Sep', '2-Mar', '2-Sep', '3-Mar', '3-Sep', '4-Mar', '4-Sep', '5-Mar', '6-Mar', '6-Sep', '7-Mar', '7-Sep', '8-Mar', '8-Sep', '9-Mar', '9-Sep']
                df = df.drop(columns=date_cols) # Drop all date values until new data is uploaded
                self._data[df_name] = df #maps dataframe name to dataframe

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data)

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN
        master_clinical = self._data['clinical'].reindex(master_index)

        # Add a column called Sample_Tumor_Normal to the clinical dataframe indicating whether each sample was a tumor or normal sample. Normal samples have a Patient_ID that begins with 'N'.
        clinical_status_col = generate_sample_status_col(master_clinical, normal_test=lambda sample: sample[0] == 'N')
        master_clinical.insert(0, "Sample_Tumor_Normal", clinical_status_col)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = master_clinical 

        # Generate a sample ID for each patient ID
        sample_id_dict = generate_sample_id_map(master_index)

        # Give every datafame a Sample_ID index
        dfs_to_delete = [] # If there's an issue reindexing a dataframe, we delete it. That shouldn't ever happen...
        for name in self._data.keys(): # Loop over the keys so we can alter the values without any issues
            df = self._data[name]
            df.index.name = "Patient_ID"
            keep_old = name in ("clinical", "treatment") # Keep the old Patient_ID index as a column in clinical and treatment, so we have a record of it.

            try:
                df = reindex_dataframe(df, sample_id_dict, "Sample_ID", keep_old)
            except ReindexMapError:
                warnings.warn(f"Error mapping sample ids in {name} dataframe. At least one Patient_ID did not have a corresponding Sample_ID mapped in clinical dataframe. {name} dataframe not loaded.", FailedReindexWarning, stacklevel=2)
                dfs_to_delete.append(name)
                continue
            self._data[name] = df

        for name in dfs_to_delete: # Delete any dataframes that had issues reindexing
            del self._data[name]

        # Set name of column axis to "Name" for all dataframes
        for name in self._data.keys(): # Loop over the keys so we can alter the values without any issues
            df = self._data[name]
            df.columns.name = "Name"
            self._data[name] = df

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
