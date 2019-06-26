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

import numpy as np
import pandas as pd
import os
import glob
from .dataset import DataSet
from .sync import get_version_files_paths

class RenalCcrcc(DataSet):

    def __init__(self, version="latest"):
        """Load all of the renalccrcc dataframes as values in the self._data dict variable, with names as keys, and format them properly."""

        # Call the parent DataSet __init__ function, which initializes self._data and other variables we need
        super().__init__()

        # Set the _cancer_type instance variable
        self._cancer_type = "renalccrcc"

        # FILL: The following overloading may or not be needed for your dataset.
        # Overload the gene separator for column names in the phosphoproteomics dataframe. In the renalccrcc data, it's an underscore, not a dash like most datasets.
        #self._gene_separator = "_"

        # FILL: If needed, overload the self._valid_omics_dfs and self._valid_metadata_dfs variables that were initialized in the parent DataSet init.

        # Get the paths to all the data files
        data_files = [] # FILL: Insert the new dataset's data file names into this list.
        data_files_paths = get_version_files_paths(self._cancer_type, version, data_files)
        if data_files_paths is None: # Version validation error. get_version_files_paths already printed an error message.
            return None

        # Load the data into dataframes in the self._data dict
        for file_path in data_files_paths: # Loops through files variable
            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            file_name_split = file_name.split(".")
            df_name = file_name_split[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            # Load the file, based on what it is
            print("Loading {} data...".format(df_name), end='\r') # Carriage return ending causes previous line to be erased.

            # FILL: Here, insert conditional statements to load all the data files as dataframes into the self._data dictionary. Consult existing datasets for examples.

            print("\033[K", end='\r') # Use ANSI escape sequence to clear previously printed line (cursor already reset to beginning of line with \r)

        print("Formatting dataframes...", end="\r")

        # FILL: Here, write code to format your dataframes properly. Requirements:
        # - All dataframes must be indexed by Sample_ID, not Patient_ID.
        #     - This means that two samples from the same patient will have the same Patient_ID, but different Sample_ID numbers.
        #     - Sample_ID numbers must be of the format S***, e.g. S001, S028, S144
        #     - clinical dataframe must contain a Patient_ID column that contains the Patient_ID for each sample
        #     - If the data did not come indexed with Sample_ID numbers, look at the Ovarian dataset for an example of generating Sample_ID numbers and mapping them to Patient_ID numbers.
        # - Each dataframe must have its name stored, all lowercase, in the DataFrame.name attribute.
        # - Each dataframe's name must match the format for that type of dataframe in all the other datasets. 
        #     - E.g., if your binary mutations dataframe is named mutations_binary, you'd need to rename it to somatic_mutation_binary to match the other datasets' binary mutation dataframes.
        # - If the new dataset has a dataframe not included in any other datasets, you must write a getter for it in the parent DataSet class, found in cptac/dataset.py
        # - You'd also need to add the new dataframe's name to self._valid_omics_dfs if it's a valid omics df for the DataSet merge functions, or self._valid_metadata_dfs if it's a valid metadata df for DataSet.append_metadata_to_omics
        # - If any dataframes are split between two files--such as one file for the tumor sample proteomics, and one file for the normal sample proteomics--they'll have been read into separate dataframes, and you need to merge those into one dataframe.
        #     - Make sure that samples coming from a normal file have an 'N' or some other indicator added to their Patient_ID numbers.
        # - If multiple dataframes are contained in one file--e.g. clinical and derived_molecular data are both in clinical.txt, as in Endometrial--separate them out here.
        # - Make sure that column names are consistent--e.g., all Patient_ID columns should be labeled as such, not as Clinical_Patient_Key or something else. Rename columns as necessary to match this.
        # - The clinical dataframe must contain a Sample_Tumor_Normal column, which contains either "Tumor" or "Normal" for each sample, according to its status.
        # - Only the clinical dataframe should contain a Patient_ID column. The other dataframes should contain just a Sample_ID index, and the data.
        # - The column axis of each dataframe should have None as the value of its .name attribute
        # - The index of each dataframe should have "Sample_ID" as the value of its .name attribute, since that's what the index is.
        # - Make sure to drop any excluded cases, as in Endometrial.
        # - Make sure that in dataframes where each column header is the name of a gene, the columns are in alphabetical order.
        # - If the dataset is still under publication embargo, print a warning after it's loaded.

        # Use ANSI escape sequence to clear previously printed line (cursor already reset to beginning of line with \r)
        print("\033[K", end='\r') 
