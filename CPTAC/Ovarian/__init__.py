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
import glob
import os
import textwrap
import webbrowser
from .dataloader import get_dataframes
from .utilities import Utilities

message = "You have loaded the CPTAC Ovarian dataset. To view available dataframes, use CPTAC.Ovarian.list_data(). To view available functions for accessing and manipulating the dataframes, use CPTAC.Ovarian.list_api()."
wrapped_list = textwrap.wrap(message)
for line in wrapped_list:
    print(line)

data_version = "Most recent release"
print("Ovarian Data Version: {}\n".format(data_version))

data = dataloader.get_dataframes()

def list_data():
    """Print a list of loaded data frames and dimensions."""
    print("Below are the available ovarian data frames contained in this package:")
    for dataframe in data:
        print("\t", data[dataframe].name)
        print("\t", "\t", "Dimensions:", data[dataframe].shape)

def list_api():
    """Print docstrings for all accessible functions."""
    help(__name__)

def get_data():
    return data

def get_clinical():
    return data.get("clinical")

def get_treatment():
    return data.get("treatment")

def get_CNV():
    return data.get("CNV")

def get_phosphoproteomics():
    return data.get("phosphoproteomics")

def get_proteomics():
    return data.get("proteomics")

def get_mutations():
    """Get the somatic_mutation dataframe."""
    return data.get("somatic_mutation")

def get_transcriptomics():
    return data.get("transcriptomics")

def get_phosphosites(genes):
    """Returns dataframe with all phosphosites for a gene or list of genes.

    Parameters:
    genes (str, or list or array-like of str): The gene(s) to get phosphosites for. str if single, list or array-like of str if multiple.

    Returns:
    pandas DataFrame: The phosphosites for the specified gene(s).
    """
    phosphoproteomics = get_phosphoproteomics()
    return Utilities().get_omics_from_str_or_list(phosphoproteomics, genes)

def get_sample_status_map():
    """Get a pandas Series from the clinical dataframe, with sample ids as the index, and each sample's status (tumor or normal) as the values."""
    clinical = get_clinical()
    map = clinical["Sample_Tumor_Normal"]
    map.name = "Sample_Status"
    return map

def compare_omics(omics_df1, omics_df2, cols1=None, cols2=None):
    """Take specified column(s) from one omics dataframe, and merge with specified columns(s) from another omics dataframe. Intersection (inner join) of indicies is used.

    Parameters:
    omics_df1 (pandas DataFrame): First omics dataframe to select columns from.
    omics_df2 (pandas DataFrame): Second omics dataframe to select columns from.
    cols1 (str, or list or array-like of str, optional): Column(s) to select from omics_df1. str if one key, list or array-like of str if multiple. Defaults to None, in which case we'll select the entire dataframe.
    cols2 (str, or list or array-like of str, optional): Column(s) to select from omics_df2. str if one key, list or array-like of str if multiple. Defaults to None, in which case we'll select the entire dataframe.

    Returns:
    pandas DataFrame: The selected columns from omics_df1 and omics_df2, merged into one dataframe.
    """
    # Make sure it's the right kind of dataframe
    valid_dfs = [
        'phosphoproteomics',
        'proteomics',
        'CNV',
        'transcriptomics']
    invalid = False
    if (omics_df1.name not in valid_dfs):
        invalid = True
        print("{} is not a valid dataframe for this function.".format(omics_df1.name))
    if (omics_df2.name not in valid_dfs):
        invalid = True
        print("{} is not a valid dataframe for this function.".format(omics_df2.name))
    if invalid:
        print("Valid dataframe options:")
        for df_name in valid_dfs:
            print('\t' + df_name)
        return

    # Return the merge.
    return Utilities().compare_omics(omics_df1, omics_df2, cols1, cols2)

def append_metadata_to_omics(metadata_df, omics_df, metadata_cols=None, omics_genes=None):
    """Append columns from metadata dataframe to part or all of an omics dataframe. Intersection (inner join) of indicies is used.

    Parameters:
    omics_df (pandas DataFrame): Omics dataframe to append the metadata columns to.
    metadata_cols (str, or list or array-like of str, optional): Column(s) to select from the metadata dataframe. str if one gene, list or array-like of str if multiple. Default of None will select all columns in the dataframe.
    omics_genes (str, or list or array-like of str, optional): Gene(s) to select data for from the omics dataframe. str if one gene, list or array-like of str if multiple. Default will select entire dataframe.

    Returns:
    pandas DataFrame: The selected metadata columns, merged with all or part of the omics dataframe.
    """
    # Make sure metadata_df is the right kind of dataframe
    valid_metadata_dfs = [
        'clinical',
        'treatment']
    if (metadata_df.name not in valid_metadata_dfs):
        print("{} is not a valid dataframe for metadata_df parameter. Valid options:".format(metadata_df.name))
        for df_name in valid_metadata_dfs:
            print('\t' + df_name)
        return

    # Make sure omics_df is the right kind of dataframe
    valid_dfs = [
        'phosphoproteomics',
        'proteomics',
        'CNV',
        'transcriptomics']
    if (omics_df.name not in valid_dfs):
        print("{} is not a valid dataframe for omics_df parameter. Valid options:".format(omics_df.name))
        for df_name in valid_dfs:
            print('\t' + df_name)
        return

    # Return the merge.
    return Utilities().append_metadata_to_omics(metadata_df, omics_df, metadata_cols, omics_genes)

def append_mutations_to_omics(omics_df, mutation_genes, omics_genes=None, show_location=True):
    """Select all mutations for specified gene(s), and append to all or part of the given omics dataframe. Intersection (inner join) of indicies is used.

    Parameters:
    omics_df (pandas DataFrame): Omics dataframe to append the mutation data to.
    mutation_genes (str, or list or array-like of str): The gene(s) to get mutation data for. str if one gene, list or array-like of str if multiple.
    omics_genes (str, or list or array-like of str, optional): Gene(s) to select from the omics dataframe. str if one gene, list or array-like of str if multiple. Default will select entire dataframe.
    show_location (bool, optional): Whether to include the Location column from the mutation dataframe. Defaults to True.

    Returns:
    pandas DataFrame: The mutations for the specified gene, appended to all or part of the omics dataframe.
    """
    # Make sure omics_df is the right kind of dataframe
    valid_dfs = [
        'phosphoproteomics',
        'proteomics',
        'CNV',
        'transcriptomics']
    if (omics_df.name not in valid_dfs):
        print("{} is not a valid dataframe for omics_df parameter. Valid options:".format(omics_df.name))
        for df_name in valid_dfs:
            print('\t' + df_name)
        return

    # Return the merge.
    mutations = get_mutations()
    sample_status_map = get_sample_status_map()
    return Utilities().append_mutations_to_omics(mutations, omics_df, mutation_genes, omics_genes, show_location, sample_status_map)

def search(term):
    """Perform an online search of the provided term.

    Parameters:
    term (str): term to be searched

    Returns: None
    """
    url = "https://www.google.com/search?q=" + term
    print("Searching for", term, "in web browser...")
    webbrowser.open(url)

def embargo():
    """Open the CPTAC embargo details in a web browser."""
    print("Opening embargo details in web browser...")
    webbrowser.open("https://proteomics.cancer.gov/data-portal/about/data-use-agreement")

def version():
    """Print the version number of the CPTAC package."""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    version = {}
    with open(dir_path + os.sep + ".." + os.sep + "version.py") as fp: #.. required to navigate up to CPTAC folder from Ovarian folder
    	exec(fp.read(), version)
    return(version['__version__'])
