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
from .dataframe_tools import *
from .exceptions import FailedReindexWarning, ReindexMapError

class Colon(Dataset):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the colon dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        valid_versions = ["0.0", "0.0.1"] # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.

        data_files = {
            "0.0": [
                "clinical.tsi.gz",
                "miRNA.cct.gz",
                "mutation_binary.cbt.gz",
                "mutation.txt.gz",
                "phosphoproteomics_normal.gz",
                "phosphoproteomics_tumor.gz",
                "proteomics_normal.cct.gz",
                "proteomics_tumor.cct.gz",
                "transcriptomics.gz"],
            "0.0.1": [
                "clinical.tsi.gz",
                "Colon_One_Year_Clinical_Data_20160927.xls",
                "Human__CPTAC_COAD__VU__SCNA__ExomeSeq__01_28_2016__BCM__Gene__BCM_CopyWriteR_GISTIC2.cct.gz",
                "miRNA.cct.gz",
                "mutation_binary.cbt.gz",
                "mutation.txt.gz",
                "phosphoproteomics_normal.gz",
                "phosphoproteomics_tumor.gz",
                "proteomics_normal.cct.gz",
                "proteomics_tumor.cct.gz",
                "transcriptomics.gz"],
        }

        super().__init__(cancer_type="colon", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            file_name_split = file_name.split(".")
            df_name = file_name_split[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            if file_name == 'Colon_One_Year_Clinical_Data_20160927.xls' and self._version == "0.0.1":
                df = pd.read_excel(file_path)

                # Replace redundant values for "not reported" with NaN
                nan_equivalents = ['Not Reported/ Unknown', 'Reported/ Unknown', 'Not Applicable',
                    'na', 'unknown', 'Not Performed', 'Unknown tumor status']

                df = df.replace(nan_equivalents, np.nan)

                # Rename and set index
                df = df.rename(columns={'PPID': 'Patient_ID'})
                df = df.set_index("Patient_ID")
                df = df.sort_index()

                self._data["followup"] = df

            elif file_name == "Human__CPTAC_COAD__VU__SCNA__ExomeSeq__01_28_2016__BCM__Gene__BCM_CopyWriteR_GISTIC2.cct.gz" and self._version == "0.0.1":
                df = pd.read_csv(file_path, sep="\t",index_col=0)
                df = df.sort_index()
                df = df.transpose()
                self._data["CNV"] = df

            else:
                df = pd.read_csv(file_path, sep="\t",index_col=0)
                df = df.sort_index()
                df = df.transpose()
                self._data[df_name] = df # Maps dataframe name to dataframe. self._data was initialized when we called the parent class __init__()

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Reformat and rename the somatic_mutation dataframe
        mut = self._data["mutation"]
        mut = mut.transpose() # Transpose it back to its original orientation
        mut = mut.sort_values(by="SampleID")
        mut = mut.reset_index()
        mut = mut[["SampleID","Gene","Variant_Type","Protein_Change"]]
        mut = mut.drop_duplicates(keep="first") # Get rid of rows that are now duplicates since we didn't keep the mRNA column. We do this before setting the index, because drop_duplicates doesn't consider the index.
        mut = mut.rename(columns={"SampleID":"Patient_ID", "Variant_Type":"Mutation", "Protein_Change":"Location"})
        mut = mut.sort_values(by=["Patient_ID", "Gene"])
        mut = mut.set_index("Patient_ID") # We only do this after the drop_duplicates call above because drop_duplicates doesn't consider the index, but we of course want the Patient_ID to be considered when identifying duplicate rows to drop.
        self._data["somatic_mutation"] = mut # Maps dataframe name to dataframe. self._data was initialized when we called the parent class __init__()
        del self._data["mutation"] # Delete the old version with the old name

        # Rename mutation_binary dataframe to somatic_mutation_binary
        self._data["somatic_mutation_binary"] = self._data["mutation_binary"]
        del self._data["mutation_binary"]

        # Separate clinical and derived molecular dataframes
        all_clinical_data = self._data.get("clinical")
        clinical_df = all_clinical_data.drop(columns=['StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'TumorPurity','immuneSubtype', 'CIN', 'Integrated.Phenotype', 'Transcriptomic_subtype', 'Proteomic_subtype', 'mutation_rate', 'Mutation_Phenotype'])
        derived_molecular_df = all_clinical_data[['StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'TumorPurity', 'immuneSubtype', 'CIN', 'Integrated.Phenotype', 'Transcriptomic_subtype', 'Proteomic_subtype', 'mutation_rate', 'Mutation_Phenotype']]

        # Format the dataframes
        clinical_df = clinical_df.apply(pd.to_numeric, errors="ignore")
        derived_molecular_df = derived_molecular_df.apply(pd.to_numeric, errors="ignore")
        derived_molecular_df = derived_molecular_df.sort_index(axis="columns")

        # Put them in our data dictionary
        self._data["clinical"] = clinical_df # Replaces original clinical dataframe
        self._data["derived_molecular"] = derived_molecular_df

        # Combine the two proteomics dataframes
        prot_tumor = self._data.get("proteomics_tumor")
        prot_normal = self._data.get("proteomics_normal") # Normal entries are already marked with 'N' on the end of the ID
        prot_combined = prot_tumor.append(prot_normal)
        self._data["proteomics"] = prot_combined
        del self._data["proteomics_tumor"]
        del self._data["proteomics_normal"]

        # Get phosphoproteomics dataframes, so we can process and combine them
        phos_tumor = self._data.get("phosphoproteomics_tumor")
        phos_normal = self._data.get("phosphoproteomics_normal") # Normal entries are not marked

        # Mark entries in phosphoproteomics_normal dataframe with an N at the end of the ID
        phos_normal = phos_normal.set_index(phos_normal.index + 'N')

        # Combine the two phosphoproteomics dataframes into one dataframe
        phos_combined = phos_tumor.append(phos_normal)

        # Create our phosphoproteomics columns multiindex
        multiindex = phos_combined.columns.str.split('[_:]', expand=True) # Split the column names into their constituent parts
        multiindex = multiindex.droplevel([2, 4]) # The third level is just empty strings, and the fifth is a duplicate of the second
        multiindex = multiindex.set_names(["Name", "Site", "Database_ID"])
        phos_combined.columns = multiindex
        phos_combined = phos_combined.sort_index(axis=1) # Put all the columns in alphabetical order
        self._data['phosphoproteomics'] = phos_combined
        del self._data["phosphoproteomics_tumor"]
        del self._data["phosphoproteomics_normal"]

        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data, exclude="followup")

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
        master_clinical = self._data['clinical'].reindex(master_index)

        # Add a column called Sample_Tumor_Normal to the clinical dataframe indicating whether each sample is a tumor or normal sample. Samples with a Patient_ID ending in N are normal.
        clinical_status_col = generate_sample_status_col(master_clinical, normal_test=lambda sample: sample[-1] == 'N')
        master_clinical.insert(0, "Sample_Tumor_Normal", clinical_status_col)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = master_clinical 

        # Edit the format of the Patient_IDs to have normal samples marked the same way as in other datasets
        # Currently, normal patient IDs have an "N" appended. We're going to make that a ".N"
        self._data = reformat_normal_patient_ids(self._data, existing_identifier="N", existing_identifier_location="end")

        # Call function from dataframe_tools.py to sort all tables first by sample status, and then by the index
        self._data = sort_all_rows(self._data)

        # Call function from dataframe_tools.py to standardize the names of the index and column axes
        self._data = standardize_axes_names(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message

    # Overload the default how_to_cite function, to provide the specific publication information for the Colon dataset
    def how_to_cite(self):
        """Print instructions for citing the data."""
        super().how_to_cite(cancer_type='colorectal cancer', pmid=31031003)