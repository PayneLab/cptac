

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
    #TODO: What is the embargo date for the ovarian cancer data?
    warning = "WARNING: This data is under a publication embargo until July 1, 2019. CPTAC is a community resource project and data are made available rapidly after generation for community research use. The embargo allows exploring and utilizing the data, but the data may not be in a publication until June 1, 2019. Please see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter embargo() to open the webpage for more details."
    wrapped_list = textwrap.wrap(warning)
    for line in wrapped_list:
        print(line)

def read_data(fileName):
    df = None
    f = fileName.split(os.sep)
    f = f[len(f) - 1]
    name = f.split(".")[0]
    print("Loading",name,"data...")
    if name == "clinical":
        df = pd.read_csv(fileName, sep="\t")
    else:
        df = pd.read_csv(fileName, sep="\t", index_col=0)
    df.name = name
    return df

dir_path = os.path.dirname(os.path.realpath(__file__))
data_directory = dir_path + os.sep + "Data" + os.sep
path = data_directory + os.sep + "*.*"
files = glob.glob(path)
data = {}
print("Loading Ovarian CPTAC data:")
for file in files:
    df = DataFrameLoader(file).createDataFrame()
    data[df.name] = df
warning()


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

def get_somatic_mutations(hg="38"):
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
    gene: gene or array of genes to select from each of the dataframes

    Returns
    Dataframe containing two columns (or number of genes provided times two). Each column is the data for the specified gene from the two specified dataframes
    """
    if isinstance(gene, str): #simple way to check for single gene string
        return Utilities().compare_gene(df1, df2, gene)
    else: #if not single gene string, then assuming an array was provided
        return Utilities().compare_genes(df1, df2, gene)


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
    version = {}
    with open(dir_path + os.sep + ".." + os.sep + "version.py") as fp: #.. required to navigate up to CPTAC folder from Endometrial folder
    	exec(fp.read(), version)
    return(version['__version__'])
