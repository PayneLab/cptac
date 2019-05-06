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
import textwrap
import glob
import pandas as pd
from .dataframe import DataFrameLoader
from .utilities import Utilities

dir_path = os.path.dirname(os.path.realpath(__file__))
data_directory = dir_path + os.sep + "Data" + os.sep
path = data_directory + "*.*"
files = glob.glob(path) #puts all files into iterable variable
data = {}
print("Loading Colon CPTAC data:")
for file in files: #loops through files variable
    try:
        df = DataFrameLoader(file).createDataFrame()
        data[df.name] = df #maps dataframe name to dataframe
    except IOError:
        print("Error reading", file)
        print("Check that all file names coincide with DataFrameLoader specs")

def list_data():
	"""
	Parameters:
	None

	Prints a list of available dataframes and dimensions

	Returns:
	None
	"""
	print("Below are the available colon data frames contained in this package:")
	#dataframes = [clinical, miRNA, mutation, proteomics, transcriptomics, phosphoproteomics]
	for dataframe in data:
		print("\t", data[dataframe].name)
		print("\t", "\t", "Dimensions:", data[dataframe].shape)
	print("To access the data, use a get function with the data frame name, i.e. colon.get_proteomics()")
def list_api():
    """
    Parameters
    None

    Prints docstrings for all accessible functions

    Returns
    None
    """
    help(__name__)

def get_clinical():
	"""
	Parameters:
	None

	Returns:
	Clinical dataframe
	"""
	return data.get("clinical")
def get_miRNA():
	"""
	Parameters:
	None

	Returns:
	miRNA dataframe
	"""
	return data.get("miRNA")
def get_mutations(binary = False):
	"""
	Parameters:
	binary: boolean value to set whether the mutation data returned is normal or binary

	Returns:
	Mutation dataframe (either binary or nonbinary)
	"""
	if binary:
		return data.get("mutation_binary")
	return data.get("mutation")
def get_mutations_binary():
    """Get the mutation_binary dataframe."""
    return data.get("mutation_binary")
def get_phosphoproteomics():
	"""
	Parameters:
	None

	Returns:
	Phosphoproteomics dataframe (both normal and tumor entries combined in one dataframe)
	"""
	tumor = data.get("phosphoproteomics_tumor")
	normal = data.get("phosphoproteomics_normal") #normal entries are not marked
	combined = tumor.append(normal)
	combined.name = 'phosphoproteomics'
	return combined
def get_proteomics():
	"""
	Parameters:
	None

	Returns:
	Proteomics dataframe (both normal and tumor entries combined in one dataframe)
	"""
	tumor = data.get("proteomics_tumor")
	normal = data.get("proteomics_normal") #normal entries are marked with 'N' on the end of the ID
	combined = tumor.append(normal)
	combined.name = "proteomics"
	return combined
def get_transcriptomics():
	"""
	Parameters:
	None

	Returns:
	Transcriptomics dataframe
	"""
	return data.get("transcriptomics")
def get_phosphosites(genes):
	"""Gets phosphosites for a gene or list of genes.

	Parameters:
	genes (str or list): gene(s) to get the phosphosites for. str if single, list of strings if multiple.

	Returns:
	pandas.core.frame.DataFrame: Phosphosites for the specified gene(s).
	"""
	phosphoproteomics = get_phosphoproteomics()
	return Utilities().get_omics_from_str_or_list(phosphoproteomics, genes)

# New merge functions
def compare_omics(omics_df1, omics_df2, cols1=None, cols2=None):
    """Take specified column(s) from one omics dataframe, and merge with specified columns(s) from another omics dataframe.

    Parameters:
    omics_df1 (pandas.core.frame.DataFrame): First omics dataframe to select columns from.
    omics_df2 (pandas.core.frame.DataFrame): Second omics dataframe to select columns from.
    cols1 (str or list, optional): Column(s) to select from omics_df1. str if one key, list if multiple. Defaults to None, in which case we'll select the entire dataframe.
    cols2 (str or list, optional): Column(s) to select from omics_df2. str if one key, list if multiple. Defaults to None, in which case we'll select the entire dataframe.

    Returns:
    pandas.core.frame.DataFrame: The selected columns from omics_df1 and omics_df2, merged into one dataframe.
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

def append_clinical_to_omics(clinical_cols, omics_df, omics_cols=None):
    """Append columns from clinical dataframe to part or all of an omics dataframe.

    Parameters:
    clinical_cols (str or list): Column(s) to select from the clinical dataframe. str if one gene, list if multiple.
    omics_df (pandas.core.frame.DataFrame): Omics dataframe to append the clinical columns to.
    omics_cols (str or list, optional): Column(s) to select from the omics dataframe. str if one gene, list if multiple. Default will select entire dataframe.

    Returns:
    pandas.core.frame.DataFrame: The selected clinical columns, merged with all or part of the omics dataframe.
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
    clinical = get_clinical()
    return Utilities().append_clinical_to_omics(clinical, omics_df, clinical_cols, omics_cols)

def append_mutations_to_omics(mutation_genes, omics_df, omics_genes=None, multiple_mutations=False, show_location=True):
    """Select all mutations for specified gene(s), and append to all or part of the given omics dataframe.

    Parameters:
    mutation_genes (str or list): The gene(s) to get mutation data for. str if one gene, list if multiple.
    omics_df (pandas.core.frame.DataFrame): Omics dataframe to append the mutation data to.
    omics_genes (str or list, optional): Gene(s) to select from the omics dataframe. str if one gene, list if multiple. Default will select entire dataframe.
    multiple_mutations (bool, optional): Whether to keep multiple mutations on the same gene for one sample, or only report the highest priority mutation.
    show_location (bool, optional): Whether to include the Location column from the mutation dataframe. Defaults to True.

    Returns:
    pandas.core.frame.DataFrame: The mutations for the specified gene, appended to all or part of the omics dataframe.
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
    return Utilities().append_mutations_to_omics(mutations, omics_df, mutation_genes, omics_genes, multiple_mutations, show_location)