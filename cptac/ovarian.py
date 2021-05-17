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
from .dataset import Dataset
from .file_download import update_index
from .dataframe_tools import *
from .exceptions import FailedReindexWarning, ReindexMapError

class Ovarian(Dataset):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the ovarian dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        valid_versions = ["0.0", "0.0.1"] # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.

        data_files = {
            "0.0": [
                "clinical.csv.gz",
                "cnv.tsv.gz",
                "definitions.txt",
                "phosphoproteomics.txt.gz",
                "proteomics.txt.gz",
                "somatic_38.maf.gz",
                "transcriptomics.tsv.gz",
                "treatment.csv.gz"],
            "0.0.1": [
                "clinical.csv.gz",
                "cnv.tsv.gz",
                "definitions.txt",
                "Ovary_One_Year_Clinical_Data_20160927.xls",
                "phosphoproteomics.txt.gz",
                "proteomics.txt.gz",
                "somatic_38.maf.gz",
                "transcriptomics.tsv.gz",
                "treatment.csv.gz"],
        }

        super().__init__(cancer_type="ovarian", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data files into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths:

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
                    df["Site"] = df["Site"].str.replace(r"[sty]", r"", regex=True) # Get rid of all lowercase s, t, and y delimeters in the sites
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

            elif file_name == 'Ovary_One_Year_Clinical_Data_20160927.xls' and self._version == "0.0.1":
                df = pd.read_excel(file_path)

                # Replace redundant values for "not reported" with NaN
                nan_equivalents = ['Not Reported/ Unknown', 'Reported/ Unknown', 'Not Applicable',
                    'na', 'unknown', 'Not Performed', 'Unknown tumor status', 'Unknown', 
                    'Unknown Tumor Status', 'Not specified']

                df = df.replace(nan_equivalents, np.nan)

                # Rename PPID to Patient_ID and set as index
                df = df.rename(columns={'PPID': 'Patient_ID'})
                df = df.set_index("Patient_ID")
                df = df.sort_index()

                self._data["followup"] = df

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data, exclude="followup")

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN
        master_clinical = self._data['clinical'].reindex(master_index)

        # Add a column called Sample_Tumor_Normal to the clinical dataframe indicating whether each sample was a tumor or normal sample. Normal samples have a Patient_ID that begins with 'N'.
        clinical_status_col = generate_sample_status_col(master_clinical, normal_test=lambda sample: sample[0] == 'N')
        master_clinical.insert(0, "Sample_Tumor_Normal", clinical_status_col)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = master_clinical 

        # Edit the format of the Patient_IDs to have normal samples marked the same way as in other datasets. Currently, all the normal samples have an "N" prepended. We're going to erase that and put a ".N" at the end.
        self._data = reformat_normal_patient_ids(self._data, existing_identifier="N", existing_identifier_location="start")

        # Call function from dataframe_tools.py to sort all tables first by sample status, and then by the index
        self._data = sort_all_rows(self._data)

        # Call function from dataframe_tools.py to standardize the names of the index and column axes
        self._data = standardize_axes_names(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message

    def how_to_cite(self):
        return super().how_to_cite(cancer_type='high grade serous ovarian cancer', pmid=27372738)