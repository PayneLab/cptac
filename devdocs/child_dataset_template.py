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

from .dataset import DataSet
from .dataframe_tools import *
from .exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError


################################################################################
# HOW TO USE THIS TEMPLATE FILE
#
# To make a class for a new dataset, copy this file and fill in the indicated
# sections, as described below.
#
# This file has sections marked with the word FILL, usually in triple quotes or
# preceded by three hashtags (###). To adapt this file for a new dataset,
# replace all of those marked fields with the proper values for your dataset.
# Additionally, there are some example code sections marked with START EXAMPLE 
# CODE and END EXAMPLE CODE. You need to replace the example code with the
# proper code for processing your dataset.
#
# This file uses dataframe processing functions imported from
# cptac/dataframe_tools.py. For more information on how to use those functions,
# you can read their docstrings in that file.
#
# If there's something confusing about this file, look at the files for existing
# datasets to provide examples of how this file would actually be implemented.
# If the new dataset you're adding has something weird that isn't addressed in
# this file, check the other datasets to see if any of them deal with a similar
# issue.
################################################################################

###FILL: Put in the actual name/acronym for the cancer type as the class name in the line below, in UpperCamelCase.
class NameOrAcronym(DataSet):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the gbm dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent DataSet class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        valid_versions = ["""FILL: Insert valid data versions here."""]

        ###FILL: Insert actual data files below
        data_files = {
            ###START EXAMPLE CODE###############################################
            "version_num": [
            "awesome_omics_data.tsv",
            "other_data_file.tsv"]
            ###END EXAMPLE CODE#################################################
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="""FILL: Insert cancer name or acronym here, in all lowercase""", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

            ###FILL: Insert if/elif statements to parse all data files. Example:
            ###START EXAMPLE CODE###############################################
            if file_name == "awesome_omics_data.tsv": # Note that we use the "file_name" variable to identify files. That way we don't have to use the whole path.
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.drop(columns=["columns", "we", "don't", "want"])
                df = df.do_some_formatting_thing()

                df = df.sort_index()
                df = df.transpose()
                self._data["awesomeomics"] = df

            elif file_name == "other_data_file.tsv":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.drop(columns=["columns", "we", "don't", "want"])
                df = df.super_crazy_dataframe_formatting_function()
                df = df.even_crazier()

                df = df.sort_index()
                df = df.transpose()
                self._data["lessawesomeomics"] = df

            ###END EXAMPLE CODE#################################################

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Get a union of all dataframes' indices, with duplicates removed
        ###FILL: If there are any tables whose index values you don't want
        ### included in the master index, pass them to the optional 'exclude'
        ### parameter of the unionize_indices function. This was useful, for
        ### example, when some datasets' followup data files included samples
        ### from cohorts that weren't in any data tables besides the followup
        ### table, so we excluded the followup table from the master index since
        ### there wasn't any point in creating empty representative rows for
        ### those samples just because they existed in the followup table.
        master_index = unionize_indices(self._data) 

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
        new_clinical = self._data["clinical"]
        new_clinical = new_clinical.reindex(master_index)

        # Add a column called Sample_Tumor_Normal to the clinical dataframe indicating whether each sample was a tumor or normal sample. Use a function from dataframe_tools to generate it.

        ###FILL: Your dataset should have some way that it marks the Patient IDs
        ### of normal samples. The example code below is for a dataset that
        ### marks them by putting an 'N' at the beginning of each one. You will
        ### need to write a lambda function that takes a given Patient_ID string
        ### and returns a bool indicating whether it corresponds to a normal
        ### sample. Pass that lambda function to the 'normal_test' parameter of
        ### the  generate_sample_status_col function when you call it. See 
        ### cptac/dataframe_tools.py for further function documentation.
        ###START EXAMPLE CODE###################################################
        sample_status_col = generate_sample_status_col(new_clinical, normal_test=lambda sample: sample[0] == 'N')
        ###END EXAMPLE CODE#####################################################

        new_clinical.insert(0, "Sample_Tumor_Normal", sample_status_col)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = new_clinical

        # Edit the format of the Patient_IDs to have normal samples marked the same way as in other datasets. 
        
        ###FILL: You will need to pass the proper parameters to correctly
        ### reformat the patient IDs in your dataset. The standard format is to
        ### have the string '.N' appended to the end of the normal patient IDs,
        ### e.g. the  normal patient ID corresponding to C3L-00378 would be
        ### C3L-00378.N (this way we can easily match two samples from the same
        ### patient). The example code below is for a dataset where all the
        ### normal samples have  an "N" prepended to the patient IDs. The
        ### reformat_normal_patient_ids function erases that and puts a ".N" at
        ### the end. See cptac/dataframe_tools.py for further function
        ### documentation.
        ###START EXAMPLE CODE###################################################
        self._data = reformat_normal_patient_ids(self._data, existing_identifier="N", existing_identifier_location="start")
        ###END EXAMPLE CODE#####################################################

        # Call function from dataframe_tools.py to sort all tables first by sample status, and then by the index
        self._data = sort_all_rows(self._data)

        # Call function from dataframe_tools.py to standardize the names of the index and column axes
        self._data = standardize_axes_names(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message

        ###FILL: If the dataset is not under publication embargo, you can remove
        ### the code block below. If it is password protected, still remove
        ### this warning, and instead keep the password protection warning
        ### below.
        # Print data embargo warning, if the date hasn't passed yet.
        today = datetime.date.today()
        embargo_date = datetime.date(year="""FILL: Insert embargo year""", month="""FILL: Insert embargo month""", day="""FILL: Insert embargo day""")
        if today < embargo_date:
            warnings.warn("The ###FILL: Insert dataset name### dataset is under publication embargo until ###FILL: Insert embargo date###. CPTAC is a community resource project and data are made available rapidly after generation for community research use. The embargo allows exploring and utilizing the data, but analysis may not be published until after the embargo date. Please see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter cptac.embargo() to open the webpage for more details.", PublicationEmbargoWarning, stacklevel=2)

        ###FILL: If the dataset is not password access only, remove the message
        ### below. If it's under publication embargo, still remove this
        ### warning, and keep the above warning about publication embargo.
        # Print password access only warning
        warnings.warn("The ###FILL: Insert dataset name### data is currently strictly reserved for CPTAC investigators. Otherwise, you are not authorized to access these data. Additionally, even after these data become publicly available, they will be subject to a publication embargo (see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter cptac.embargo() to open the webpage for more details).", PublicationEmbargoWarning, stacklevel=2)
