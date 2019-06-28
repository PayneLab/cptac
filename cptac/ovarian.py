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
import glob
import textwrap
import datetime
from .dataset import DataSet
from .sync import get_version_files_paths
from .dataframe_tools import *

class Ovarian(DataSet):

    def __init__(self, version="latest"):
        """Load all the ovarian dataframes as values in the self._data dict variable, with names as keys, and format them properly."""

        # Call the parent Dataset __init__() function, which initializes self._data and other variables we need
        super().__init__()

        # Set the cancer_type instance variable 
        self._cancer_type = "ovarian"

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
        data_files_paths = get_version_files_paths(self._cancer_type, version, data_files)
        if data_files_paths is None: # Version validation error. get_version_files_paths already printed an error message.
            return None

        # Load the data files into dataframes in the self._data dict
        for file_path in data_files_paths:
            path_elements = file_path.split(os.sep) # Get a list of all the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            # Load the file, based on what it is
            print("Loading {} data...".format(df_name), end='\r') # Carriage return ending causes previous line to be erased.

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
                df = pd.read_csv(file_path, sep="\t", index_col = 0)
                if file_name == "proteomics.txt.gz":
                    df = df[df["hgnc_symbol"].notnull()] # Drops all nan values in hgnc_symbol column
                    df = df.set_index("hgnc_symbol")
                elif file_name == "phosphoproteomics.txt.gz":
                    df = df[df["site"].notnull()] # Drops all nan values in site column
                    df = df.drop(["refseq_peptide","Peptide"],axis=1)
                    df = df.set_index("site")

                df = df.sort_index()
                df = df.transpose()
                c_index = df.index.values.tolist()
                index_no_c = [id[1:] if id.startswith("C") else id for id in c_index]
                index_no_c = pd.Index(index_no_c)
                df = df.set_index(index_no_c) # Take C prefix off of indicies for those samples that have them (tumor samples have C, normal have N)
                full_index = df.index.values.tolist()
                ids_to_drop = [id for id in full_index if id.startswith('OV_QC')]
                df = df.drop(ids_to_drop) # Drop all OV_QC* samples--they're quality control samples not relevant for data analysis
                self._data[df_name] = df #maps dataframe name to dataframe

            elif file_name == "somatic_38.maf.gz":
                df = pd.read_csv(file_path, sep = "\t", index_col=0)
                df = df.reset_index()
                split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n = 1, expand = True) # The first part of the barcode is the patient id, which we need to make a Patient_ID column
                df["Tumor_Sample_Barcode"] = split_barcode[0]
                parsed_df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]] # We only want these columns
                parsed_df = parsed_df.rename(columns={"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"})
                parsed_df = parsed_df.set_index("Patient_ID")
                self._data['somatic_mutation'] = parsed_df #maps dataframe name to dataframe

            elif file_name == "transcriptomics.tsv.gz":
                df = pd.read_csv(file_path, sep="\t", index_col=0)
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                date_cols = ['1-Dec', '1-Sep', '10-Mar', '10-Sep', '11-Sep', '12-Sep', '14-Sep', '15-Sep', '2-Mar', '2-Sep', '3-Mar', '3-Sep', '4-Mar', '4-Sep', '5-Mar', '6-Mar', '6-Sep', '7-Mar', '7-Sep', '8-Mar', '8-Sep', '9-Mar', '9-Sep']
                df = df.drop(columns=date_cols) # Drop all date values until new data is uploaded
                self._data[df_name] = df #maps dataframe name to dataframe

            print("\033[K", end='\r') # Use ANSI escape sequence to clear previously printed line (cursor already reset to beginning of line with \r)

        print("Formatting dataframes...", end="\r")

        # Get a union of all dataframes' indicies, with duplicates removed
        master_index = unionize_indicies(self._data)

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
        for name in self._data.keys(): # Only loop over keys, to avoid changing the structure of the object we're looping over
            df = self._data[name]
            df.index.name = "Patient_ID"
            keep_old = name in ("clinical", "treatment") # Keep the old Patient_ID index as a column in clinical and treatment, so we have a record of it.

            df = reindex_dataframe(df, sample_id_dict, "Sample_ID", keep_old)
            if df is None:
                print(f"Error mapping sample ids in {name} dataframe. At least one Patient_ID did not have a corresponding Sample_ID mapped in clinical dataframe. {name} dataframe not loaded.")
                del self._data[name]
                continue

            self._data[name] = df

        # Drop name of column axis for all dataframes
        for name in self._data.keys():
            df_rename_col_axis = self._data[name]
            df_rename_col_axis.columns.name = None
            self._data[name] = df_rename_col_axis

        # Use ANSI escape sequence to clear previously printed line (cursor already reset to beginning of line with \r)
        print("\033[K", end='\r') 
