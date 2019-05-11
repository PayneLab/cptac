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

def set_sample_id_index(df, sample_id_dict, drop_patient_ids): # private
    """Replaces a patient ID index with a sample ID index for a dataframe.

    Parameters:
    df (pandas DataFrame): The dataframe to change the index of.
    sample_id_dict (dict): A dictionary where the keys are the patient ids, and the values are the corresponding sample ids. Every value in the original dataframe's index must match a key in this dictionary.
    drop_patient_ids (bool): Whether to drop the patient id column after it is no longer the index.

    Returns:
    pandas DataFrame: The original dataframe, with the patient id index replaced with a sample id index.
    """
    sample_id_column = []
    for row in df.index:
        sample_id_column.append(sample_id_dict[row]) # This is why every value in the original dataframe's index has to match a key in sample_id_dict
    return_df = df.assign(Sample_ID=sample_id_column)
    if not drop_patient_ids:
        old_index_name = df.index.name
        if old_index_name is None:
            old_index_name = 'index'
        return_df = return_df.reset_index() # This gives the dataframe a default numerical index and makes the old index a column, which prevents it from being dropped when we set Sample_ID as the index.
        return_df = return_df.rename(columns={old_index_name:'Patient_ID'}) # Rename the old index as Patient_ID
    return_df = return_df.set_index('Sample_ID') # Make the Sample_ID column the index
    return_df.name = df.name
    return return_df

dir_path = os.path.dirname(os.path.realpath(__file__)) #gets path to CPTAC package
data_directory = dir_path + os.sep + "Data" + os.sep #appends Data to path
path = data_directory + os.sep + "*.*" #appends "*.*" to path, which looks for all files
files = glob.glob(path) #puts all files into iterable variable

# Load the dataframes
print("Loading Ovarian CPTAC data:")
data = {}
for file in files: #loops through files variable
    df = DataFrameLoader(file).createDataFrame()
    data[df.name] = df #maps dataframe name to dataframe

# Get a union of all dataframes' indicies, with duplicates removed
indicies = [df.index for df in data.values()]
master_index = pd.Index([])
for index in indicies:
    master_index = master_index.union(index)
master_index = master_index.drop_duplicates()

# Generate a sample ID for each patient ID
sample_id_dict = {}
for i in range(len(master_index)):
    patient_id = master_index[i]
    sample_id_dict[patient_id] = "S{:0>3}".format(i + 1) # Use string formatter to give each sample id the format S*** filled with zeroes, e.g. S001 or S023

# Put a mapping in the clinical dataframe of all patient ids to their sample ids, including patient ids for samples not originally in the clinical dataframe. Then, give the clinical dataframe a sample id index.
master_df = pd.DataFrame(index=master_index) 
master_clinical = data['clinical'].join(master_df, how='outer') # Do an outer join with the clinical dataframe, so that clinical has a row for every sample in the dataset
master_clinical = set_sample_id_index(master_clinical, sample_id_dict, drop_patient_ids=False) # Replace the patient id index with a sample id index in the clinical dataframe. Keep the patient ids so we can maps sample ids to their patient ids.
data['clinical'] = master_clinical # Replace the clinical dataframe in the data dictionary with our new and improved version!

# Give the other dataframes Sample_ID indicies, but don't keep the old index, since we have a mapping in the clinical dataframe of all sample ids to their patient ids.
for df in data.keys(): # Only loop over keys, to avoid changing the structure of the object we're looping over
    if df != 'clinical':
        data[df] = set_sample_id_index(data[df], sample_id_dict, drop_patient_ids=True)

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
