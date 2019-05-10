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
from .utilities import Utilities
from .dataframe import DataFrameLoader


def warning():
    print("\n","******PLEASE READ******")
    warning = "WARNING: This data is under a publication embargo until June 1, 2019. CPTAC is a community resource project and data are made available rapidly after generation for community research use. The embargo allows exploring and utilizing the data, but the data may not be in a publication until June 1, 2019. Please see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter embargo() to open the webpage for more details."
    wrapped_list = textwrap.wrap(warning)
    for line in wrapped_list:
        print(line)

dir_path = os.path.dirname(os.path.realpath(__file__)) #gets path to CPTAC package
data_directory = dir_path + os.sep + "Data" + os.sep #appends Data to path
path = data_directory + os.sep + "*.*" #appends "*.*" to path, which looks for all files
files = glob.glob(path) #puts all files into iterable variable
data = {}
print("Loading Ovarian CPTAC data:")
for file in files: #loops through files variable
    df = DataFrameLoader(file).createDataFrame()
    data[df.name] = df #maps dataframe name to dataframe
    
warning() #displays warning

def list_data():
    """
    Parameters
    None

    Prints list of loaded data frames and dimensions

    Returns
    None
    """

    print("Below are the available ovarian data frames contained in this package:")
    for dataframe in data:
        print("\t", data[dataframe].name)
        print("\t", "\t", "Dimensions:", data[dataframe].shape)
def list_api():
    """
    Parameters
    None

    Prints docstrings for all accessible functions

    Returns
    None
    """
    help(__name__)
def get_data():
    return data
def get_clinical():
    return data.get("clinical")

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

def embargo():
    """
    Parameters
    None

    Opens CPTAC embargo details in web browser

    Returns
    None
    """
    print("Opening embargo details in web browser...")
    webbrowser.open("https://proteomics.cancer.gov/data-portal/about/data-use-agreement")
def version():
    """
    Parameters
    None

    Prints version number of CPTAC package

    Returns
    Version number
    """
    version = {}
    with open(dir_path + os.sep + ".." + os.sep + "version.py") as fp: #.. required to navigate up to CPTAC folder from Endometrial folder
    	exec(fp.read(), version)
    return(version['__version__'])

# New merge functions
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

def append_clinical_to_omics(omics_df, clinical_cols=None, omics_cols=None):
    """Append columns from clinical dataframe to part or all of an omics dataframe. Intersection (inner join) of indicies is used.

    Parameters:
    omics_df (pandas DataFrame): Omics dataframe to append the clinical columns to.
    clinical_cols (str, or list or array-like of str, optional): Column(s) to select from the clinical dataframe. str if one gene, list or array-like of str if multiple. Default of None will select all columns in the dataframe.
    omics_cols (str, or list or array-like of str, optional): Column(s) to select from the omics dataframe. str if one gene, list or array-like of str if multiple. Default will select entire dataframe.

    Returns:
    pandas DataFrame: The selected clinical columns, merged with all or part of the omics dataframe.
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
    clinical = get_clinical()
    return Utilities().append_clinical_to_omics(clinical, omics_df, clinical_cols, omics_cols)

def append_mutations_to_omics(omics_df, mutation_genes, omics_genes=None, multiple_mutations=False, show_location=True):
    """Select all mutations for specified gene(s), and append to all or part of the given omics dataframe. Intersection (inner join) of indicies is used.

    Parameters:
    omics_df (pandas DataFrame): Omics dataframe to append the mutation data to.
    mutation_genes (str, or list or array-like of str): The gene(s) to get mutation data for. str if one gene, list or array-like of str if multiple.
    omics_genes (str, or list or array-like of str, optional): Gene(s) to select from the omics dataframe. str if one gene, list or array-like of str if multiple. Default will select entire dataframe.
    multiple_mutations (bool, optional): Whether to keep multiple mutations on the same gene for one sample, or only report the highest priority mutation.
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
    return Utilities().append_mutations_to_omics(mutations, omics_df, mutation_genes, omics_genes, multiple_mutations, show_location)
