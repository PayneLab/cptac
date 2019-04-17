

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
def align_indices(data): #private
    all_index = data.get("clinical").index #grab the biggest index
    prot_index = data.get("proteomics").index
    c_prot = prot_index[0:83].str[1:] #drops letter off all indices with "C"
    all_index = all_index.append(c_prot).unique().sort_values() #sort all unique indices of all indices
    all_index = all_index.append(prot_index[83:]) #append all "N" indices
    patient_key = []
    for i in range(0,len(all_index)): #make S number for every index
        s = "S" + str(i + 1)
        patient_key.append(s)
    dict = {"patient_id":all_index,"patient_key":patient_key}
    key_id_map = pd.DataFrame(dict).set_index("patient_id") #map patient id to S number (called patient key)
    return key_id_map
def map_ids(data, key_id_map): #private, adds patient key (S number) to every data set for merging
    #transcriptomics, clinical, cnv, and somatic mutation data all match normally
    tran = data.get("transcriptomics")
    tran["patient_key"] = key_id_map.loc[tran.index]
    data["transcriptomics"] = tran

    clin = data.get("clinical")
    clin["patient_key"] = key_id_map.loc[clin.index]
    data["clinical"] = clin

    cnv = data.get("cnv")
    cnv["patient_key"] = key_id_map.loc[cnv.index]
    data["cnv"] = cnv

    som38 = data.get("somatic_38")
    som38["patient_key"] = key_id_map.loc[som38.index]
    som19 = data.get("somatic_19")
    som19["patient_key"] = key_id_map.loc[som19.index]
    data["somatic_38"] = som38
    data["somatic_19"] = som19

    #because of "C" and "N" prefix, proteomics and phosphoproteomics indices are mapped by dropping the "C" to find the propper patient key (S number)
    prot = data.get("proteomics")
    c = key_id_map.loc[prot.index.str[1:][0:83]]
    n = key_id_map.loc[prot.index[83:]]
    prot["patient_key"] = np.concatenate([c,n])
    data["proteomics"] = prot

    pho = data.get("phosphoproteomics")
    c = key_id_map.loc[pho.index.str[1:][0:83]]
    n = key_id_map.loc[pho.index[83:]]
    pho["patient_key"] = np.concatenate([c,n])
    data["phosphoproteomics"] = pho

    return data

dir_path = os.path.dirname(os.path.realpath(__file__)) #gets path to CPTAC package
data_directory = dir_path + os.sep + "Data" + os.sep #appends Data to path
path = data_directory + os.sep + "*.*" #appends "*.*" to path, which looks for all files
files = glob.glob(path) #puts all files into iterable variable
data = {}
print("Loading Ovarian CPTAC data:")
for file in files: #loops through files variable
    df = DataFrameLoader(file).createDataFrame()
    data[df.name] = df #maps dataframe name to dataframe
key_id_map = align_indices(data) #creates map for all unique indices to S number for merging
data = map_ids(data, key_id_map) #adds S number to all indicies for all data
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
    print("To access the data, use a get function with the data frame name, i.e. ovarian.get_proteomics()")
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

def get_cnv():
    return data.get("cnv")

def get_phosphoproteomics():
    return data.get("phosphoproteomics")

def get_proteomics():
    return data.get("proteomics")

def get_mutations(hg="38"):
    if hg == "38":
        return data.get("somatic_38") #Defaulting to somatic_38
    elif hg == "19":
        return data.get("somatic_19")
    else:
        print("Somatic",hg,"data not found, available options are somatic_19 and somatic_38")

def get_transcriptomics():
    return data.get("transcriptomics")

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
        return Utilities().compare_gene(df1, df2, gene, key_id_map)
    else: #if not single gene string, then assuming an array was provided
        return Utilities().compare_genes(df1, df2, gene, key_id_map)

def compare_clinical(omics_data, clinical_col):
    """
    Parameters
    data: omics data for clinical data to be appended with
    clinical_col: column in clinical dataframe to be inserted into provided omics data

    Returns
    Dataframe with specified column from clinical dataframe added to specified dataframe (i.e., proteomics) for comparison and easy plotting
    """
    return Utilities().compare_clinical(get_clinical(), omics_data, clinical_col, key_id_map)
def get_phosphosites(gene):
    """Returns dataframe with all phosphosites of specified gene name"""
    return Utilities().get_phosphosites(get_phosphoproteomics(), gene)
def compare_phosphosites(gene):
    """
    Parameters
    gene: proteomics gene to query phosphoproteomics dataframe

    Searches for any phosphosites on the gene provided

    Returns
    Dataframe with a column from proteomics for the gene specified, as well as columns for all phosphoproteomics columns beginning with the specified gene
    """
    return Utilities().compare_phosphosites(get_proteomics(), get_phosphoproteomics(), gene)
def compare_mutations(omics_data, omics_gene, mutations_gene = None):
    """
    Params
    omics_data: omics dataframe (i.e. proteomics, phosphoproteomics, transcriptomics)
    omics_gene: gene to select from omics data (used for somatic data if somatic_gene is left blank)
    mutations_gene: gene to select from somatic mutation data

    Returns
    Dataframe containing two columns, the omics data and the somatic mutation type for the gene(s) provided
    """
    if mutations_gene:
        return Utilities().merge_mutations_trans(omics_data, omics_gene, get_mutations(), mutations_gene, key_id_map)
    else:
        return Utilities().merge_mutations(omics_data, get_mutations(), omics_gene, key_id_map)
def compare_mutations_full(omics_data, omics_gene, mutations_gene = None):#doesn't work right now due to duplicate indices messing up the key_id_map
    """
    Params
    omics_data: omics dataframe (i.e. proteomics, phosphoproteomics, transcriptomics)
    omics_gene: gene to select from omics data (used for somatic data if somatic_gene is left blank)
    mutations_gene: gene to select from somatic mutation data

    Returns
    Dataframe containing numeric omics data and categorical somatic data (including patient ID, mutation type, and mutation location)
    """
    if mutations_gene:
        return Utilities().merge_mutations_trans(omics_data, omics_gene, get_mutations(), mutations_gene,  key_id_map, duplicates = True)
    else:
        return Utilities().merge_mutations(omics_data, get_mutations(), omics_gene, key_id_map, duplicates = True)

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
