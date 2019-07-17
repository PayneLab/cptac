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
from .dataset import DataSet
from .file_download import update_index
from .file_tools import validate_version, get_version_files_paths
from .dataframe_tools import *

class Colon(DataSet):

    def __init__(self, version="latest"):
        """Load all of the colon dataframes as values in the self._data dict variable, with names as keys, and format them properly."""

        # Call the parent DataSet __init__ function, which initializes self._data and other variables we need
        super().__init__("colon")

        # Update the index, if possible. If there's no internet, update_index will return False, but we don't care in this context.
        update_index(self._cancer_type)

        # Validate the index
        self._version = validate_version(version, self._cancer_type, use_context="init")
        if self._version is None: # Validation error. validate_version already printed an error message.
            return

        # Overload the gene separator for column names in the phosphoproteomics dataframe. In the colon data, it's an underscore, not a dash like most datasets.
        self._gene_separator = "_"

        # Get the paths to all the data files
        data_files = [
            "clinical.tsi.gz",
            "miRNA.cct.gz",
            "mutation_binary.cbt.gz",
            "mutation.txt.gz",
            "phosphoproteomics_normal.gz",
            "phosphoproteomics_tumor.gz",
            "proteomics_normal.cct.gz",
            "proteomics_tumor.cct.gz",
            "transcriptomics.gz"]
        data_files_paths = get_version_files_paths(self._cancer_type, self._version, data_files)
        if data_files_paths is None: # Data error. get_version_files_paths already printed an error message.
            return

        # Load the data into dataframes in the self._data dict
        loading_msg = "Loading dataframes"
        for file_path in data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            file_name_split = file_name.split(".")
            df_name = file_name_split[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            df = pd.read_csv(file_path, sep="\t",index_col=0)
            df = df.transpose()
            self._data[df_name] = df # Maps dataframe name to dataframe. self._data was initialized when we called the parent class __init__()

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Reformat and rename the somatic_mutation dataframe
        mut = self._data["mutation"]
        mut = mut.transpose() # Transpose it back to its original orientation
        mut = mut.sort_values(by="SampleID")
        mut = mut[["SampleID","Gene","Variant_Type","Protein_Change"]]
        mut = mut.rename(columns={"SampleID":"Patient_ID", "Variant_Type":"Mutation", "Protein_Change":"Location"})
        mut = mut.sort_values(by=["Patient_ID", "Gene"])
        mut = mut.set_index("Patient_ID")
        self._data["somatic_mutation"] = mut # Maps dataframe name to dataframe. self._data was initialized when we called the parent class __init__()
        del self._data["mutation"] # Delete the old version with the old name

        # Rename mutation_binary dataframe to somatic_mutation_binary
        self._data["somatic_mutation_binary"] = self._data["mutation_binary"]
        del self._data["mutation_binary"]

        # Separate clinical and derived molecular dataframes
        all_clinical_data = self._data.get("clinical")
        clinical_df = all_clinical_data.drop(columns=['StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'TumorPurity','immuneSubtype', 'CIN', 'Integrated.Phenotype'])
        derived_molecular_df = all_clinical_data[['StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'TumorPurity', 'immuneSubtype', 'CIN', 'Integrated.Phenotype']]

        # Format clinical dataframe
        clinical_df = clinical_df.astype({"Age": float, "CEA": float, "mutation_rate": float}) # For one reason or another, these weren't automatically cast on loading.

        # Put them in our data dictionary
        self._data["clinical"] = clinical_df # Replaces original clinical dataframe
        self._data["derived_molecular"] = derived_molecular_df

        # Combine the two proteomics dataframes
        prot_tumor = self._data.get("proteomics_tumor")
        prot_normal = self._data.get("proteomics_normal") # Normal entries are marked with 'N' on the end of the ID
        prot_combined = prot_tumor.append(prot_normal)
        self._data["proteomics"] = prot_combined
        del self._data["proteomics_tumor"]
        del self._data["proteomics_normal"]

        # Get phosphoproteomics dataframes, so we can process and combine them
        phos_tumor = self._data.get("phosphoproteomics_tumor")
        phos_normal = self._data.get("phosphoproteomics_normal") # Normal entries are not marked

        # Mark entries in phosphoproteomics_normal dataframe with an N at the end of the ID, to match proteomics_normal
        phos_normal_indices = phos_normal.index.values.tolist()
        for i in range(len(phos_normal_indices)):
            index = phos_normal_indices[i]
            index_marked = index + 'N'
            phos_normal_indices[i] = index_marked
        new_phos_index = pd.Index(phos_normal_indices)
        phos_normal = phos_normal.set_index(new_phos_index)

        # Combine the two phosphoproteomics dataframes into one dataframe
        phos_combined = phos_tumor.append(phos_normal)
        phos_combined = phos_combined.rename(columns=lambda x: x.split(":")[0]) # Drop everything after ":" in column names--unneeded additional identifiers
        phos_combined = phos_combined.sort_index(axis=1) # Put all the columns in alphabetical order
        self._data['phosphoproteomics'] = phos_combined
        del self._data["phosphoproteomics_tumor"]
        del self._data["phosphoproteomics_normal"]

        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data)

        # Sort this master_index so all the samples with an N suffix are last. Because the N is a suffix, not a prefix, this is kind of messy.
        status_df = pd.DataFrame(master_index, columns=['Patient_ID']) # Create a new dataframe with the master_index as a column called "Patient_ID"
        status_col = []
        for index in master_index:
            if index[-1] == 'N':
                status_col.append("Normal")
            else:
                status_col.append("Tumor")
        status_df = status_df.assign(Status=status_col)
        status_df = status_df.sort_values(by=['Status', 'Patient_ID'], ascending=[False, True]) # Sorts first by status, and in descending order, so "Tumor" samples are first
        master_index = status_df["Patient_ID"].tolist()

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
        master_clinical = self._data['clinical'].reindex(master_index)

        # Add a column called Sample_Tumor_Normal to the clinical dataframe indicating whether each sample is a tumor or normal sample. Samples with a Patient_ID ending in N are normal.
        clinical_status_col = generate_sample_status_col(master_clinical, normal_test=lambda sample: sample[-1] == 'N')
        master_clinical.insert(0, "Sample_Tumor_Normal", clinical_status_col)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = master_clinical 

        # Generate a sample ID for each patient ID
        sample_id_dict = generate_sample_id_map(master_index)

        # Give all the dataframes Sample_ID indices
        dfs_to_delete = [] # If there's an issue reindexing a dataframe, we delete it. That shouldn't ever happen...
        for name in self._data.keys(): # Loop over the keys so we can alter the values without any issues
            df = self._data[name]
            df.index.name = "Patient_ID"
            keep_old = name == "clinical" # Keep the old Patient_ID index as a column in the clinical dataframe, so we have a record of it.

            df = reindex_dataframe(df, sample_id_dict, "Sample_ID", keep_old)
            if df is None:
                print(f"Error mapping sample ids in {name} dataframe. At least one Patient_ID did not have corresponding Sample_ID mapped in clinical dataframe. {name} dataframe not loaded.")
                dfs_to_delete.append(name)
                continue

            self._data[name] = df

        for name in dfs_to_delete: # Delete any dataframes that had issues reindexing
            del self._data[name]

        # Drop name of column axis for all dataframes
        for name in self._data.keys(): # Loop over the keys so we can alter the values without any issues
            df = self._data[name]
            df.columns.name = None
            self._data[name] = df

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message

    # Overload the default how_to_cite function, to provide the specific publication information for the Colon dataset
    def how_to_cite(self):
        """Print instructions for citing the data."""
        print('Please include the following statement in publications using colon cancer data accessed through this module:\n"Data used in this publication were generated by the Clinical Proteomic Tumor Analysis Consortium (NCI/NIH, <https://proteomics.cancer.gov/programs/cptac/>). Data were published at <https://www.ncbi.nlm.nih.gov/pubmed/31031003>. Data were accessed through the Python module cptac, available at <https://pypi.org/project/cptac/>."')
