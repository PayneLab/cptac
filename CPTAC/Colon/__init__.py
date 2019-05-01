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

def list():
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
def get_phosphosites(gene):
	"""
	Parameters
	gene: gene to get the phosphosites for

	Returns
	dataframe containing the phosphosites for the specified gene
	"""
	phosphoproteomics = get_phosphoproteomics()
	return Utilities().get_phosphosites(phosphoproteomics, gene)
def compare_gene(df1, df2, gene):
    """
    Parameters
    df1: omics dataframe (proteomics) to be selected from
    df2: other omics dataframe (transcriptomics) to be selected from
    gene: gene or list of genes to select from each of the dataframes

    Returns
    Dataframe containing common rows between provided dataframes and columns for the specified gene (or genes) from provided dataframes.
    """
    if isinstance(gene, str): #simple way to check for single gene string
        return Utilities().compare_gene(df1, df2, gene)
    else: #if not single gene string, then assuming an array was provided
        return Utilities().compare_genes(df1, df2, gene)
def compare_clinical(omics_data, clinical_col):
    """
    Parameters
    data: omics data for clinical data to be appended with
    clinical_col: column in clinical dataframe to be inserted into provided omics data

    Returns
    Dataframe with specified column from clinical dataframe added to specified dataframe (i.e., proteomics) for comparison and easy plotting
    """
    return Utilities().compare_clinical(get_clinical(), omics_data, clinical_col)
def compare_phosphosites(gene):
	"""
	Parameters
	proteomics: the proteomics dataframe
	phosphoproteomics: the phosphoproteomics dataframe
	gene: proteomics gene to query phosphoproteomics dataframe

	Searches for any phosphosites on the gene provided

	Returns
	Dataframe with a column from proteomics for the gene specified, as well as columns for all phosphoproteomics columns beginning with the specified gene
	"""
	proteomics = get_proteomics()
	phosphoproteomics = get_phosphoproteomics()
	return Utilities().compare_phosphosites(proteomics, phosphoproteomics, gene)
def compare_mutations(omics_data, omics_gene, mutations_gene = None):
    """
    Params
    omics_data: omics dataframe (i.e. proteomics, phosphoproteomics, transcriptomics)
    omics_gene: gene to select from omics data (used for mutation data if mutations_gene is left blank)
    mutations_gene: gene to select from somatic mutation data

    Returns
    Dataframe containing two columns, the omics data and the somatic mutation type for the gene(s) provided
    """
    if mutations_gene: #compare omics data of omics gene to mutations of mutations_gene
        return Utilities().merge_mutations_trans(omics_data, omics_gene, get_mutations(), mutations_gene)
    else: #compare omics data to mutations for same gene
        return Utilities().merge_mutations(omics_data, get_mutations(), omics_gene)
def compare_mutations_full(omics_data, omics_gene, mutations_gene = None):
    """
    Params
    omics_data: omics dataframe (i.e. proteomics, phosphoproteomics, transcriptomics)
    omics_gene: gene to select from omics data (used for somatic data if somatic_gene is left blank)
    mutations_gene: gene to select from somatic mutation data

    Returns
    Dataframe containing numeric omics data and categorical somatic data (including patient ID, mutation type, and mutation location)
    """
    if mutations_gene: #compare omics data of omics gene to mutations of mutations_gene
        return Utilities().merge_mutations_trans(omics_data, omics_gene, get_mutations(), mutations_gene, duplicates = True)
    else: #compare omics data to mutations for same gene
        return Utilities().merge_mutations(omics_data, get_mutations(), omics_gene, duplicates = True)

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
        'acetylproteomics',
        'proteomics',
        'transcriptomics_linear', # But not transcriptomics_circular or miRNA--they have incompatible column names.
        'CNA',
        'phosphoproteomics_site',
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
        'acetylproteomics',
        'proteomics',
        'transcriptomics_linear', # But not transcriptomics_circular or miRNA--they have incompatible column names.
        'CNA',
        'phosphoproteomics_site',
        'phosphoproteomics_gene']
    if (omics_df.name not in valid_dfs):
        print("{} is not a valid dataframe for omics_df parameter. Valid options:".format(omics_df.name))
        for df_name in valid_dfs:
            print('\t' + df_name)
        return

    # Return the merge.
    return Utilities().append_clinical_or_derived_molecular_to_omics(clinical, omics_df, clinical_cols, omics_cols)

def append_derived_molecular_to_omics(derived_molecular_cols, omics_df, omics_cols=None):
    """Append columns from derived_molecular dataframe to all or part of an omics dataframe.

    Parameters:
    derived_molecular_cols (str or list): Column(s) to select from the derived_molecular dataframe. str if one gene, list if multiple.
    omics_df (pandas.core.frame.DataFrame): Omics dataframe to append the derived_molecular columns to.
    omics_cols (str or list, optional): Column(s) to select from the omics dataframe. str if one gene, list if multiple. Default will select entire dataframe.

    Returns:
    pandas.core.frame.DataFrame: The selected derived_molecular columns, merged with all or part of the omics dataframe.
    """
    # Make sure omics_df is the right kind of dataframe
    valid_dfs = [
        'acetylproteomics',
        'proteomics',
        'transcriptomics_linear', # But not transcriptomics_circular or miRNA--they have incompatible column names.
        'CNA',
        'phosphoproteomics_site',
        'phosphoproteomics_gene']
    if (omics_df.name not in valid_dfs):
        print("{} is not a valid dataframe for omics_df parameter. Valid options:".format(omics_df.name))
        for df_name in valid_dfs:
            print('\t' + df_name)
        return

    # Return the merge.
    return Utilities().append_clinical_or_derived_molecular_to_omics(derived_molecular, omics_df, derived_molecular_cols, omics_cols)

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
        'acetylproteomics',
        'proteomics',
        'transcriptomics_linear', # But not transcriptomics_circular or miRNA--they have incompatible column names.
        'CNA',
        'phosphoproteomics_site',
        'phosphoproteomics_gene']
    if (omics_df.name not in valid_dfs):
        print("{} is not a valid dataframe for omics_df parameter. Valid options:".format(omics_df.name))
        for df_name in valid_dfs:
            print('\t' + df_name)
        return

    # Return the merge.
    return Utilities().append_mutations_to_omics(somatic_maf, omics_df, mutation_genes, omics_genes, multiple_mutations, show_location)