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
import datetime

from cptac.cancers.source import Source
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError


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
### For example, the endometrial dataset's class is called Endometrial; the BRCA dataset's class is called Brca; and the ccRCC dataset's class is called Ccrcc.
class NameOrAcronym(Source):

    def __init__(self, version="latest", no_internet=False):
        """Define load functions for all of the datatypes this source provides. in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        self.valid_versions = ["""FILL: Insert valid data versions here."""]

        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        ###FILL: Insert actual data files below
        self.data_files = {
            ###START EXAMPLE CODE###############################################
            "version_num": {
                "awesome_omics"  : "awesome_omics_data.tsv",
                "other_data"     : "other_data_file.tsv"
            }
            ###END EXAMPLE CODE#################################################
        }

        ###FILL: Again, insert actual data files below
        self.load_functions = {
            ### START EXAMPLE CODE ###############################################
            # Any load functions that are common to all versions can be added here
            # If a datatype is only available in a specific version, check for
            # that after initializing this dictionary and add those load functions
            'awesome_omics' : self.load_awesome_omics,
            'other_data' : self.load_other_data,
            # Basically any dataframe a user would want to access should have a load
            # function here, and nothing else. (e.g. load annotation does not belong)
            ### END EXAMPLE CODE #################################################
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="""FILL: Insert cancer acronym here, in all lowercase""", source="""FILL: Insert source here, in all lowercase""", version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

        ###FILL: Insert load functions to parse all data files. Example:
        ###START EXAMPLE CODE###############################################
    def load_awesome_omics(self):
        df_type = 'awesome_omics'
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.drop(columns=["columns", "we", "don't", "want"])
            df = df.do_some_formatting_thing()

            df = df.sort_index()
            df = df.transpose()
            # save df in self._data, basically 'self._data["awesomeomics"] = df' but with extra formatting
            self.save_df(df_type, df)

    def load_other_data(self):
        df_type = 'other_data'
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.super_crazy_dataframe_formatting_function()
            df = df.even_crazier()

            df = df.sort_index()
            df = df.transpose()
            self.save_df(df_type, df)

        # If you have any doubts, look through some of the other child sources
        # They have examples of many different situations, including combining
        # multiple files into one datatype, or splitting one file into
        # several datatypes.
        ###END EXAMPLE CODE#################################################


        ### NOTE: Any normal samples in your dataset will need to be marked.
        ### Throughout cptac, all normal sample patient ids end with .N
        ### For example: "C3N-01234.N"
        ### Sometimes the normal samples are already marked in the original data
        ### in some way, but not the way we want (e.g. they have an "N" at the
        ### beginning of the sample ID, instead of a ".N" at the end).
        ### Some datasets such as PDAC are different; instead of the normal
        ### samples already being marked, just not in the way we want, they're
        ### actually contained in a separate table, with no special marking on
        ### the sample ids. In those cases you would just mark the samples in
        ### the normal tables with the ".N" before appending them to the tumor tables.
        ### Depending on how the original data is formatted you could either do this
        ### in the individual load functions, or for every dataset in the save_df
        ### function, as will be demonstrated below


    ### The parent class (source.py) has a save_df function that each load function uses
    ### This function does some standard formatting and saves the parsed dataframe in a
    ### data dictionary. If there is special formatting that each load function in a
    ### particular source needs, the best way to do that is to override the save_df function.
    ### If that isn't needed, you can delete this
    # Override the save_df function from source.py
    def save_df(self, datatype, df):
        ### Here is an example of something you wouldn't want to do in each load function
        # Drop samples C3N.00545 and C3N.00545.N from the dataset.
        # They were excluded due to poor sample quality (see data freeze README; excluded in data freeze 3.0)
        cases_to_drop = ["C3N.00545", "C3N.00545.N"]
        df = df.drop(index=cases_to_drop, errors="ignore")

        # Replace periods with hyphens in all Patient_IDs
        df.index.name = "Patient_ID"
        df = df.reset_index()

        # Replace all '.' with '-'
        df["Patient_ID"] = df["Patient_ID"].str.replace(r"\.", "-", regex=True)
        # If there's a "-N" at the end, it's part of the normal identifier, which we want to actually be ".N"
        df["Patient_ID"] = df["Patient_ID"].str.replace(r"-N$", ".N", regex=True)

        # Set the index back to Patient_ID
        df = df.set_index("Patient_ID")

        ### We can check if no data has been loaded for this source yet, so we
        ### Don't display a warning until the user tries to use the data
        ### and we don't display the warning every time they use it.
        if self._data == {}:
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

        # Inherit the parent event
        super().save_df(datatype, df)
