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

import os
import sys
import webbrowser
import textwrap
import pandas as pd
import numpy as np
from .dataloader import get_dataframes
from .dataloader import get_dictionary
from .utilities import Utilities

message = "You have loaded the cptac endometrial dataset. To view available dataframes, use cptac.endometrial.list_data(). To view available functions for accessing and manipulating the dataframes, use cptac.endometrial.list_api()."
wrapped_list = textwrap.wrap(message)
for line in wrapped_list:
    print(line)

data_version = "2.1"
print("endometrial data version: {}\n".format(data_version))

dictionary = dataloader.get_dictionary()
data = dataloader.get_dataframes()

def list_data():
    """Print list of loaded data frames and dimensions."""
    print("Below are the available endometrial data frames contained in this package:")
    for dataframe in data.values():
        print("\t", dataframe.name)
        print("\t", "\t", "Dimensions:", dataframe.shape)

def list_api():
    """Print docstrings for all accessible functions."""
    help(__name__)

def get_clinical():
    """Get clinical dataframe."""
    return data["clinical"]

def get_derived_molecular():
    """Get derived_molecular dataframe."""
    return data["derived_molecular"]

def get_experimental_setup():
    """Get experimental_setup dataframe."""
    return data["experimental_setup"]

def get_acetylproteomics():
    """Get acetylproteomics dataframe."""
    return data["acetylproteomics"]

def get_proteomics():
    """Get proteomics dataframe."""
    return data["proteomics"]

def get_transcriptomics():
    """Gets transcriptomics dataframe."""
    return data["transcriptomics"]

def get_circular_RNA():
    """Gets circular_RNA dataframe."""
    return data["circular_RNA"]

def get_miRNA():
    """Gets miRNA dataframe."""
    return data["miRNA"]

def get_CNA():
    """Get the CNA dataframe."""
    return data["CNA"]

def get_phosphoproteomics():
    """Gets the phosphoproteomics dataframe."""
    return data["phosphoproteomics"]

def get_phosphoproteomics_gene():
    """Gets the phosphoproteomics_gene dataframe. The gene level phosphorylation measurement is an aggregate metric which potentially averages together individual measurements of different sites. Use get_phosphoproteomics() to view the data for individual sites."""
    return data["phosphoproteomics_gene"]

def get_phosphosites(genes):
    """Returns dataframe with all phosphosites of specified gene or list of genes.

    Parameters:
    genes (str, or list or array-like of str): gene or list of genes to use to select phosphosites. str if single, list or array-like of str if multiple.

    Returns:
    pandas DataFrame: The phosphoproteomics for the specified gene(s).
    """
    phosphoproteomics = get_phosphoproteomics()
    return Utilities().get_omics_cols(phosphoproteomics, genes)

def get_mutations():
    """Get the somatic_mutation dataframe."""
    return data["somatic_mutation"]

def get_mutations_binary():
    """Gets the somatic_mutation_binary dataframe, which has a binary value indicating, for each location on each gene, whether there was a mutation in that gene at that location, for each sample."""
    return data["somatic_mutation_binary"]

def get_sample_status_map():
    """Get a pandas Series from the clinical dataframe, with sample ids as the index, and each sample's status (tumor or normal) as the values."""
    clinical = get_clinical()
    raw_map = clinical["Proteomics_Tumor_Normal"]
    parsed_map = raw_map.where(raw_map == "Tumor", other="Normal") # Replace various types of normal (Adjacent_normal, Myometrium_normal, etc.) with just "Normal"
    parsed_map.name = "Sample_Status"
    return parsed_map

def compare_omics(omics_df1, omics_df2, cols1=None, cols2=None):
    """Take specified column(s) from one omics dataframe, and append to specified columns(s) from another omics dataframe. Intersection (inner join) of indicies is used.

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
        'acetylproteomics',
        'proteomics',
        'transcriptomics', # But not circular_RNA or miRNA--they have incompatible column names.
        'CNA',
        'phosphoproteomics',
        'phosphoproteomics_gene']
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
    """Joins columns from a metadata dataframe (clinical, derived_molecular, or experimental_setup) to part or all of an omics dataframe. Intersection (inner join) of indicies is used.

    Parameters:
    metadata_df (pandas DataFrame): Metadata dataframe to select columns from. Either clinical, derived_molecular, or experimental_setup.
    omics_df (pandas DataFrame): Omics dataframe to append the metadata columns to.
    metadata_cols (str, or list or array-like of str, optional): Column(s) to select from the metadata dataframe. str if one gene, list or array-like of str if multiple. Default is None, which will select the entire metadata dataframe.
    omics_genes (str, or list or array-like of str, optional): Gene(s) to select data for from the omics dataframe. str if one gene, list or array-like of str if multiple. Default is None, which will select entire dataframe.

    Returns:
    pandas DataFrame: The selected metadata columns, merged with all or part of the omics dataframe.
    """
    # Make sure metadata_df is the right kind of dataframe
    valid_metadata_dfs = [
        'clinical',
        'derived_molecular',
        'experimental_setup']
    if (metadata_df.name not in valid_metadata_dfs):
        print("{} is not a valid dataframe for metadata_df parameter. Valid options:".format(metadata_df.name))
        for df_name in valid_metadata_dfs:
            print('\t' + df_name)
        return

    # Make sure omics_df is the right kind of dataframe
    valid_omics_dfs = [
        'acetylproteomics',
        'proteomics',
        'transcriptomics', # But not circular_RNA or miRNA--they have incompatible column names.
        'CNA',
        'phosphoproteomics',
        'phosphoproteomics_gene']
    if (omics_df.name not in valid_omics_dfs):
        print("{} is not a valid dataframe for omics_df parameter. Valid options:".format(omics_df.name))
        for df_name in valid_omics_dfs:
            print('\t' + df_name)
        return

    # Return the merge.
    return Utilities().append_metadata_to_omics(metadata_df, omics_df, metadata_cols, omics_genes)

def append_mutations_to_omics(omics_df, mutation_genes, omics_genes=None, show_location=True):
    """Select all mutations for specified gene(s), and appends them to all or part of the given omics dataframe. Intersection (inner join) of indicies is used. Each location or mutation cell contains a list, which contains the one or more location or mutation values corresponding to that sample for that gene, or a value indicating that the sample didn't have a mutation in that gene.

    Parameters:
    omics_df (pandas DataFrame): Omics dataframe to append the mutation data to.
    mutation_genes (str, or list or array-like of str): The gene(s) to get mutation data for. str if one gene, list or array-like of str if multiple.
    omics_genes (str, or list or array-like of str, optional): Gene(s) to select from the omics dataframe. str if one gene, list or array-like of str if multiple. Default will select entire dataframe.
    show_location (bool, optional): Whether to include the Locations column from the mutation dataframe. Defaults to True.

    Returns:
    pandas DataFrame: The mutations for the specified gene, appended to all or part of the omics dataframe. Each location or mutation cell contains a list, which contains the one or more location or mutation values corresponding to that sample for that gene, or a value indicating that the sample didn't have a mutation in that gene.
    """
    # Make sure omics_df is the right kind of dataframe
    valid_dfs = [
        'acetylproteomics',
        'proteomics',
        'transcriptomics', # But not circular_RNA or miRNA--they have incompatible column names.
        'CNA',
        'phosphoproteomics',
        'phosphoproteomics_gene']
    if (omics_df.name not in valid_dfs):
        print("{} is not a valid dataframe for omics_df parameter. Valid options:".format(omics_df.name))
        for df_name in valid_dfs:
            print('\t' + df_name)
        return

    # Return the merge.
    somatic_mutation = get_mutations()
    sample_status_map = get_sample_status_map()
    return Utilities().append_mutations_to_omics(somatic_mutation, omics_df, mutation_genes, omics_genes, show_location, sample_status_map)

def define(term):
    """Define a term, if it is in the dataset's definitions dictionary.

    Parameters:
    term (str): term to be defined

    Returns:
    str: definition of provided term
    """
    if term in dictionary:
        print(dictionary[term])
    else:
        print(term, "not found in dictionary. Alternatively, cptac.define() can be used to perform a web search of the term provided.")

def search(term):
    """Search for a term in a web browser.

    Parameters:
    term (str): term to be searched

    Returns: None
    """
    url = "https://www.google.com/search?q=" + term
    print("Searching for", term, "in web browser...")
    webbrowser.open(url)

def embargo():
    """Open CPTAC embargo details in web browser."""
    print("Opening embargo details in web browser...")
    webbrowser.open("https://proteomics.cancer.gov/data-portal/about/data-use-agreement")

def version():
    """Print version number of cptac package."""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    version = {}
    with open(dir_path + os.sep + ".." + os.sep + "version.py") as fp: #.. required to navigate up to cptac folder from endometrial folder
    	exec(fp.read(), version)
    return(version['__version__'])
