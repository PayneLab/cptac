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

# Combine the two proteomics dataframes
prot_tumor = data.get("proteomics_tumor")
prot_normal = data.get("proteomics_normal") #normal entries are marked with 'N' on the end of the ID
prot_combined = tumor.append(normal)
prot_combined.name = "proteomics"
data[prot_combined.name] = prot_combined
del data["proteomics_tumor"]
del data["proteomics_normal"]

# Get phosphoproteomics dataframes
phos_tumor = data.get("phosphoproteomics_tumor")
phos_normal = data.get("phosphoproteomics_normal") # Normal entries are not marked

# Mark entries in normal dataframe
phos_normal_indicies = phos_normal.index.values.tolist()
for i in range(len(phos_normal_indicies)):
    index = phos_normal_indicies[i]
    index_marked = index + 'N'
    phos_normal_indicies[i] = index_marked
phos_normal = phos_normal.set_index(phos_normal_indicies)

# Combine the two phosphoproteomics dataframes
phos_combined = phos_tumor.append(phos_normal)
phos_combined.name = 'phosphoproteomics'
data[phos_combined.name] = phos_combined
del data["phosphoproteomics_tumor"]
del data["phosphoproteomics_normal"]


def list_data():
	"""
	Parameters:
	None

	Prints a list of available dataframes and dimensions

	Returns:
	None
	"""
	print("Below are the available colon data frames contained in this package:")
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
	
def get_mutations():
	"""Get the somatic_mutation dataframe."""
	return data.get("somatic_mutation")

def get_mutations_binary():
    """Get the somatic_mutation_binary dataframe."""
    return data.get("somatic_mutation_binary")

def get_phosphoproteomics():
	"""
	Parameters:
	None

	Returns:
	Phosphoproteomics dataframe (both normal and tumor entries combined in one dataframe)
	"""
	return data.get("phosphoproteomics")

def get_proteomics():
	"""
	Parameters:
	None

	Returns:
	Proteomics dataframe (both normal and tumor entries combined in one dataframe)
	"""
	return data.get("proteomics")

def get_transcriptomics():
	"""
	Parameters:
	None

	Returns:
	Transcriptomics dataframe
	"""
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

def append_clinical_to_omics(omics_df, clinical_cols=None, omics_cols=None):
    """Append columns from clinical dataframe to part or all of an omics dataframe. Intersection (inner join) of indicies is used.

    Parameters:
    omics_df (pandas DataFrame): Omics dataframe to append the clinical columns to.
    clinical_cols (str, or list or array-like of str, optional): Column(s) to select from the clinical dataframe. str if one gene, list or array-like of str if multiple. Default of None will cause entire dataframe to be selected.
    omics_cols (str, or list or array-like of str, optional): Column(s) to select from the omics dataframe. str if one gene, list or array-like of str if multiple. Default will select entire dataframe.

    Returns:
    pandas DataFrame: The selected clinical columns, merged with all or part of the omics dataframe.
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
        'transcriptomics']
    if (omics_df.name not in valid_dfs):
        print("{} is not a valid dataframe for omics_df parameter. Valid options:".format(omics_df.name))
        for df_name in valid_dfs:
            print('\t' + df_name)
        return

    # Return the merge.
    mutations = get_mutations()
    return Utilities().append_mutations_to_omics(mutations, omics_df, mutation_genes, omics_genes, multiple_mutations, show_location)