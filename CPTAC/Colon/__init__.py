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
def get_mutation(binary = False):
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
        return Utilities().merge_mutations_trans(omics_data, omics_gene, get_mutation(), mutations_gene)
    else: #compare omics data to mutations for same gene
        return Utilities().merge_mutations(omics_data, get_mutation(), omics_gene)
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
        return Utilities().merge_mutations_trans(omics_data, omics_gene, get_mutation(), mutations_gene, duplicates = True)
    else: #compare omics data to mutations for same gene
        return Utilities().merge_mutations(omics_data, get_mutation(), omics_gene, duplicates = True)
