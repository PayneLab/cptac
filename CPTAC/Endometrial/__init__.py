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
import pandas as pd
from .dataframe import DataFrameLoader
from .meta import MetaData
from .molecular import MolecularData
from .utilities import Utilities
from .queries import Queries

def warning():
    print("\n","******PLEASE READ******")
    warning = "WARNING: This data is under a publication embargo until July 1, 2019. CPTAC is a community resource project and data are made available rapidly after generation for community research use. The embargo allows exploring and utilizing the data, but the data may not be in a publication until July 1, 2019. Please see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter CPTAC.embargo() to open the webpage for more details."
    wrapped_list = textwrap.wrap(warning)
    for line in wrapped_list:
        print(line)

"""
Creates dictionary for linking Patient_Id with individual sample number (i.e. C3L-00006 with S001)
"""
def create_patient_ids(clinical):
    c = clinical[["Proteomics_Participant_ID"]][0:100] #S101 through S140 have no patient id
    s = c.index
    dictPrepDf = c.set_index('Proteomics_Participant_ID')
    dictPrepDf['idx'] = s
    patient_ids = dictPrepDf.to_dict()['idx']
    return patient_ids
def link_patient_ids(patient_ids, somatic):
    s = []
    for x in somatic["Patient_Id"]:
        if x in patient_ids.keys():
            s.append(patient_ids[x])
        else:
            s.append("NA")
    somatic["Clinical_Patient_Key"] = s
    return somatic
"""
Executes on import CPTAC statement. Selects files from docs folder in CPTAC package
utilizing DataFrameLoader from dataframe.py. Prints update as files are loaded into
dataframes.
"""
print("Loading CPTAC data:")

dir_path = os.path.dirname(os.path.realpath(__file__))
data_directory = dir_path + os.sep + "Data" + os.sep

print("Loading Dictionary...")
dict = {}
file = open(data_directory + "definitions.txt", "r")

for line in file:
    line = line.strip()
    line = line.split("\t")
    dict[line[0]] = line[1]
file.close()

print("Loading Clinical Data...")
clinical = DataFrameLoader(data_directory + "clinical.txt").createDataFrame()

print("Loading Proteomics Data...")
proteomics = DataFrameLoader(data_directory + "proteomics.cct.gz").createDataFrame()

print("Loading Transcriptomics Data...")
transcriptomics = DataFrameLoader(data_directory + "transcriptomics_linear.cct.gz").createDataFrame()
transcriptomics_circular = DataFrameLoader(data_directory + "transcriptomics_circular.cct.gz").createDataFrame()
miRNA = DataFrameLoader(data_directory + "miRNA.cct.gz").createDataFrame()

print("Loading CNA Data...")
cna = DataFrameLoader(data_directory + "CNA.cct.gz").createDataFrame()

print("Loading Phosphoproteomics Data...")
phosphoproteomics = DataFrameLoader(data_directory + "phosphoproteomics_site.cct.gz").createDataFrame()
phosphoproteomics_gene = DataFrameLoader(data_directory + "phosphoproteomics_gene.cct.gz").createDataFrame()

print("Loading Somatic Mutation Data...")
somatic_binary = DataFrameLoader(data_directory + "somatic.cbt.gz").createDataFrame()
somatic_binary.name = "somatic binary"
somatic_unparsed = pd.read_csv(data_directory + "somatic.maf.gz", sep = "\t")
somatic_unparsed.name = "somatic MAF unparsed"
somatic_maf = DataFrameLoader(data_directory + "somatic.maf.gz").createDataFrame()
patient_ids = create_patient_ids(clinical)
somatic_maf = link_patient_ids(patient_ids, somatic_maf)
somatic_maf.name = "somatic MAF"


#metaData = MetaData(clinical)
#molecularData = MolecularData(proteomics, transcriptome, cna, phosphoproteomics)
warning()
def list():
    """
    Parameters
    None

    Prints list of loaded data frames and dimensions

    Returns
    None
    """
    print("Below are the available data frames contained in this package:")
    data = [clinical, proteomics, transcriptomics, cna, phosphoproteomics, somatic_binary, somatic_maf]
    for dataframe in data:
        print("\t", dataframe.name)
        print("\t", "\t", "Dimensions:", dataframe.shape)
    print("To access the data, use a get function with the data frame name, i.e. CTPAC.get_proteomics()")
def define(term):
    """
    Parameters
    term: string of term to be defined

    Returns
    String definition of provided term
    """
    if term in dict:
        print(dict[term])
    else:
        print(term, "not found in dictionary. Alternatively, CPTAC.define() can be used to perform a web search of the term provided.")
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
def get_clinical():
    """
    Parameters
    None

    Returns
    Clincal dataframe
    """
    return clinical
def get_proteomics():
    """
    Parameters
    None

    Returns
    proteomics dataframe
    """
    return proteomics
def get_transcriptomics(circular = False, miRNA = False):
    """
    Parameters
    circular: boolean indicating whether to return circular RNA data
    miRNA: boolean indicating whether to return miRNA data

    Returns
    Transcriptomics dataframe
    """
    if circular:
        return transcriptomics_circular
    if miRNA:
        return miRNA
    return transcriptomics
def get_CNA():
    """
    Parameters
    None

    Returns
    CNA dataframe
    """
    return cna
def get_phosphoproteomics(gene_level = False):
    """
    Parameters
    gene_level: boolean indicating whether to return gene level phosphoproteomics (returns site level if false)

    Returns
    Phosphoproteomics dataframe
    """
    if gene_level:
        return phosphoproteomics_gene
    return phosphoproteomics
def get_phosphosites(gene):
    """Returns dataframe with all phosphosites of specified gene name"""
    return Utilities().get_phosphosites(phosphoproteomics, gene)
def get_somatic(binary=False, unparsed=False):
    """
    Parameters
    binary: boolean indicating whether to retrieve the somatic mutations binary data
    unparsed: boolean indicating whether to retrieve unparsed somatic mutations maf data

    Default behavior is to return parsed somatic mutations maf data

    Returns
    Somatic mutations dataframe corresponding with parameters provided
    """
    if binary:
        return somatic_binary
    if unparsed:
        return somatic_unparsed
    return somatic_maf
def get_clinical_cols():
    """
    Parameters
    None

    Returns
    List of clincal dataframe columns, aka data types (i.e. BMI, Diabetes)
    """
    return clinical.columns
def get_proteomics_cols():
    """
    Parameters
    None

    Returns
    List of columns of proteomics dataframe
    """
    return proteomics.columns
def get_transcriptomics_cols():
    """
    Parameters
    None

    Returns
    List of columns of transcriptomics dataframe
    """
    return transcriptomics.columns
def get_cohort_clinical(columns):
    """
    Parameters
    columns: single column name or array of column names to select for in the clinical dataframe

    Returns
    Dataframe of specified columns (or Series if one column) of clinical data
    """
    return clinical[columns]
def get_proteomics_quant(colon_ids):
    """
    Parameters
    colon_ids: string or list of string ids (i.e. S001, S068) to be selected from proteomics dataframe

    Returns
    Dataframe of specified rows (or Series if one row) of proteomics data
    """
    return proteomics.loc[colon_ids]
def get_cohort_proteomics(columns):
    """
    Parameters
    columns: single column name or array of column names to select for in the proteomics dataframe

    Returns
    Dataframe of specified columns (or Series if one column) of proteomics data
    """
    return proteomics[columns]
def get_cohort_transcriptomics(columns):
    """
    Parameters
    columns: single column name or array of column names to select for in the transcriptomics dataframe

    Returns
    Dataframe of specified columns (or Series if one column) of transcriptomics data
    """
    return transcriptomics[columns]
def get_cohort_cna(columns):
    """
    Parameters
    columns: single column name or array of column names to select for in the CNA dataframe

    Returns
    Dataframe of specified columns (or Series if one column) of CNA data
    """
    return cna[columns]
def get_cohort_phosphoproteomics(columns):
    """
    Parameters
    columns: single column name or array of column names to select for in the phosphoproteomics dataframe

    Returns
    Dataframe of specified columns (or Series if one column) of phosphoproteomics data
    """
    return phosphoproteomics[columns]
def get_patient_mutations(patient_id):
    """
    Parameters
    patient_id: Patient ID (i.e. C3L-00006, S018, etc.) to select from somatic mutation data

    Returns
    Dataframe containing data for provided patient ID
    """
    if len(patient_id) == 4: #S***
        return somatic_maf[somatic_maf["Patient_Id"] == patient_id]
    elif len(patient_id) > 0: #C3L-*****
        return somatic_maf[somatic_maf["Clinical_Patient_Key"] == patient_id]
    else:
        print("ERROR:", patient_id, "not a valid patient_id.")
def get_tumor_ids(tumor_type, query_type, value): #TODO: implement
    #"""
    #Parameters
    #tumor_type is the tumor type, e.g. colon
    #query_type is the type of tumor query, e.g. by SNP, mutated gene, outlier
    #value corresponds with the query type, e.g. TP53 for mutated gene or EGFR for outlier

    #Returns

    #"""
    dataframe = None #TODO what should the dataframe be?
    return Queries(dataframe).query(tumor_type, query_type, value)
def get_gene_mapping():
    #TODO implement
    return Utilities().get_gene_mapping()
def convert(snp_or_sap):
    #TODO implement
    return Utilities().convert(snp_or_sap)
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
        return Utilities().merge_mutations_trans(omics_data, omics_gene, somatic_maf, mutations_gene)
    else:
        return Utilities().merge_mutations(omics_data, somatic_maf, omics_gene)
def compare_mutations_full(omics_data, omics_gene, mutations_gene = None):
    """
    Params
    omics_data: omics dataframe (i.e. proteomics, phosphoproteomics, transcriptomics)
    omics_gene: gene to select from omics data (used for somatic data if somatic_gene is left blank)
    mutations_gene: gene to select from somatic mutation data

    Returns
    Dataframe containing numeric omics data and categorical somatic data (including patient ID, mutation type, and mutation location)
    """
    if mutations_gene:
        return Utilities().merge_mutations_trans(omics_data, omics_gene, somatic_maf, mutations_gene, duplicates = True)
    else:
        return Utilities().merge_mutations(omics_data, somatic_maf, omics_gene, duplicates = True)
def compare_clinical(omics_data, clinical_col):
    """
    Parameters
    data: omics data for clinical data to be appended with
    clinical_col: column in clinical dataframe to be inserted into provided omics data

    Returns
    Dataframe with specified column from clinical dataframe added to specified dataframe (i.e., proteomics) for comparison and easy plotting
    """
    #TODO: do we need clinical parameter? Could just grab it from loaded data?
    return Utilities().compare_clinical(clinical, omics_data, clinical_col)
def compare_phosphosites(gene):
    """
    Parameters
    gene: proteomics gene to query phosphoproteomics dataframe

    Searches for any phosphosites on the gene provided

    Returns
    Dataframe with a column from proteomics for the gene specified, as well as columns for all phosphoproteomics columns beginning with the specified gene
    """
    return Utilities().compare_phosphosites(proteomics, phosphoproteomics, gene)

def help():
    """
    Parameters
    None

    Opens github help page

    Returns
    None
    """
    print("Opening help.txt in web browser...")
    webbrowser.open("https://github.com/PayneLab/CPTAC/blob/master/doc/help.txt")
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
def start():
    #Might remove this function
    print("Welcome to our CPTAC data. Enter CPTAC.help() to open our Github help page.")
def version():
    version = {}
    with open(dir_path + '/version.py') as fp:
    	exec(fp.read(), version)
    return(version['__version__'])
