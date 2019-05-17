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
import webbrowser
import pandas as pd
import numpy as np
from .dataloader import get_dataframes
from .utilities import Utilities

data = dataloader.get_dataframes()

def list_data():
	"""Print a list of available dataframes and their dimensions."""
	print("Below are the available colon data frames contained in this package:")
	for dataframe in data:
		print("\t", data[dataframe].name)
		print("\t", "\t", "Dimensions:", data[dataframe].shape)

def list_api():
    """Print docstrings for all accessible functions."""
    help(__name__)

def get_clinical():
	"""Get the clinical dataframe."""
	return data.get("clinical")

def get_derived_molecular():
    """Get the derived_molecular dataframe."""
    return data.get("derived_molecular")

def get_miRNA():
	"""Get the miRNA dataframe."""
	return data.get("miRNA")
	
def get_mutations():
	"""Get the somatic_mutation dataframe."""
	return data.get("somatic_mutation")

def get_mutations_binary():
    """Get the somatic_mutation_binary dataframe."""
    return data.get("somatic_mutation_binary")

def get_phosphoproteomics():
	"""Get the phosphoproteomics dataframe (both normal and tumor entries combined in one dataframe)."""
	return data.get("phosphoproteomics")

def get_proteomics():
	"""Get the proteomics dataframe (both normal and tumor entries combined in one dataframe)."""
	return data.get("proteomics")

def get_transcriptomics():
	"""Get the transcriptomics dataframe."""
	return data.get("transcriptomics")

def get_phosphosites(genes):
	"""Gets phosphosites for a gene or list or array-like of str of genes.

	Parameters:
	genes (str, or list or array-like of str): gene(s) to get the phosphosites for. str if single, list or array-like of str if multiple.

	Returns:
	pandas DataFrame: Phosphosites for the specified gene(s).
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

def append_metadata_to_omics(metadata_df, omics_df, metadata_cols=None, omics_cols=None):
    """Joins columns from a metadata dataframe (clinical or derived_molecular) to part or all of an omics dataframe. Intersection (inner join) of indicies is used.

    Parameters:
    metadata_df (pandas DataFrame): Metadata dataframe to select columns from. Either clinical or derived_molecular.
    omics_df (pandas DataFrame): Omics dataframe to append the metadata columns to.
    metadata_cols (str, or list or array-like of str, optional): Column(s) to select from the metadata dataframe. str if one gene, list or array-like of str if multiple. Default is None, which will select the entire metadata dataframe.
    omics_cols (str, or list or array-like of str, optional): Column(s) to select from the omics dataframe. str if one gene, list or array-like of str if multiple. Default is None, which will select entire dataframe.

    Returns:
    pandas DataFrame: The selected metadata columns, merged with all or part of the omics dataframe.
    """
    # Make sure metadata_df is the right kind of dataframe
    valid_metadata_dfs = [
        'clinical',
        'derived_molecular']
    if (metadata_df.name not in valid_metadata_dfs):
        print("{} is not a valid dataframe for metadata_df parameter. Valid options:".format(metadata_df.name))
        for df_name in valid_metadata_dfs:
            print('\t' + df_name)
        return

    # Make sure omics_df is the right kind of dataframe
    valid_dfs = [
        'phosphoproteomics',
        'proteomics',
        'transcriptomics']
    if (omics_df.name not in valid_dfs):
        print("{} is not a valid dataframe for omics_df parameter. Valid options:".format(omics_df.name))
        for df_name in valid_dfs:
            print('\t' + df_name)
        return

    # Return the merge.
    return Utilities().append_metadata_to_omics(metadata_df, omics_df, metadata_cols, omics_cols)

def append_mutations_to_omics(omics_df, mutation_genes, omics_genes=None, show_location=True):
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
    """
    Parameters
    term: string of term to be searched

    Performs online search of provided term

    Returns
    None
    """
    url = "https://www.google.com/search?q=" + term
    print("Searching for", term, "in web browser...")
    webbrowser.open(url)

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
    dir_path = os.path.dirname(os.path.realpath(__file__))
    version = {}
    with open(dir_path + os.sep + ".." + os.sep + "version.py") as fp: #.. required to navigate up to CPTAC folder from Endometrial folder, TODO: how to navigate from dataTest.py?
        exec(fp.read(), version)
    return(version['__version__'])
