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
from .dataframe import DataFrameLoader
from .meta import MetaData
from .molecular import MolecularData
from .utilities import Utilities
from .queries import Queries

def warning():
    print("\n","******PLEASE READ******")
    warning = "WARNING: This data is under a publication embargo until July 1, 2019. CPTAC is a community resource project and data are made available rapidly after generation for community research use. The embargo allows exploring and utilizing the data, but the data may not be in a publication until July 1, 2019. Please see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter embargo() to open the webpage for more details."
    wrapped_list = textwrap.wrap(warning)
    for line in wrapped_list:
        print(line)

"""
Creates dictionary for linking Patient_Id with individual sample number (i.e. C3L-00006 with S001)
"""
def create_patient_ids(clinical): #private
    c = clinical[["Proteomics_Participant_ID"]][0:103] # S105 maps back to S001
    s = c.index
    dictPrepDf = c.set_index('Proteomics_Participant_ID')
    dictPrepDf['idx'] = s
    patient_ids = dictPrepDf.to_dict()['idx']
    return patient_ids
def link_patient_ids(patient_ids, somatic): #private
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
print("Loading Endometrial CPTAC data:")

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
clinical_file_data = DataFrameLoader(data_directory + "clinical.txt").createDataFrame()
casesToDrop = clinical_file_data[clinical_file_data["Case_excluded"] == "Yes"].index
clinical_unfiltered = clinical_file_data[[
    'Proteomics_Participant_ID', 'Case_excluded',  'Proteomics_Tumor_Normal',  'Country',
    'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity',
    'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN', 'Clin_Stage_Dist_Mets-cM', 'Path_Stage_Dist_Mets-pM',
    'tumor_Stage-Pathological', 'FIGO_stage', 'LVSI', 'BMI', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site',
    'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm',   'Num_full_term_pregnancies']]
clinical = clinical_unfiltered.drop(casesToDrop, errors = "ignore") #Drops all samples with Case_excluded == Yes
clinical = clinical.drop(['Case_excluded'], axis=1)
clinical_unfiltered.name = "clinical"
clinical.name = clinical_unfiltered.name
derived_molecular_u = clinical_file_data.drop(['Proteomics_Participant_ID', 'Case_excluded',  'Proteomics_Tumor_Normal',  'Country',
    'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity',
    'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN', 'Clin_Stage_Dist_Mets-cM', 'Path_Stage_Dist_Mets-pM',
    'tumor_Stage-Pathological', 'FIGO_stage', 'LVSI', 'BMI', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site',
    'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm',   'Num_full_term_pregnancies'], axis=1)
derived_molecular = derived_molecular_u.drop(casesToDrop, errors = "ignore")
derived_molecular_u.name = "derived_molecular"
derived_molecular.name = derived_molecular_u.name

print("Loading Acetylation Proteomics Data...")
acetylproteomics_u = DataFrameLoader(data_directory + "acetylproteomics.cct").createDataFrame()
acetylproteomics = acetylproteomics_u.drop(casesToDrop, errors = "ignore")
acetylproteomics.name = acetylproteomics_u.name

print("Loading Proteomics Data...")
proteomics_u = DataFrameLoader(data_directory + "proteomics.cct.gz").createDataFrame()
proteomics = proteomics_u.drop(casesToDrop, errors = "ignore")
proteomics.name = proteomics_u.name

print("Loading Transcriptomics Data...")
transcriptomics_u = DataFrameLoader(data_directory + "transcriptomics_linear.cct.gz").createDataFrame()
transcriptomics_circular_u = DataFrameLoader(data_directory + "transcriptomics_circular.cct.gz").createDataFrame()
miRNA_u = DataFrameLoader(data_directory + "miRNA.cct.gz").createDataFrame()

transcriptomics = transcriptomics_u.drop(casesToDrop, errors = "ignore")
transcriptomics_circular = transcriptomics_circular_u.drop(casesToDrop, errors = "ignore")
miRNA = miRNA_u.drop(casesToDrop, errors = "ignore")

transcriptomics.name = transcriptomics_u.name
transcriptomics_circular.name = transcriptomics_circular_u.name
miRNA.name = miRNA_u.name

print("Loading cna Data...")
cna_u = DataFrameLoader(data_directory + "cna.cct.gz").createDataFrame()
cna = cna_u.drop(casesToDrop, errors = "ignore")
cna.name = cna_u.name

print("Loading Phosphoproteomics Data...")
phosphoproteomics_u = DataFrameLoader(data_directory + "phosphoproteomics_site.cct.gz").createDataFrame()
phosphoproteomics_gene_u = DataFrameLoader(data_directory + "phosphoproteomics_gene.cct.gz").createDataFrame()

phosphoproteomics = phosphoproteomics_u.drop(casesToDrop, errors = "ignore")
phosphoproteomics_gene = phosphoproteomics_gene_u.drop(casesToDrop, errors = "ignore")
phosphoproteomics.name = phosphoproteomics_u.name
phosphoproteomics_gene.name = phosphoproteomics_gene_u.name

print("Loading Somatic Mutation Data...")
somatic_binary_u = DataFrameLoader(data_directory + "somatic.cbt.gz").createDataFrame()
somatic_binary = somatic_binary_u.drop(casesToDrop, errors = "ignore")
somatic_binary.name = "somatic binary"
somatic_unparsed_u = pd.read_csv(data_directory + "somatic.maf.gz", sep = "\t")
somatic_unparsed = somatic_unparsed_u.drop(casesToDrop, errors = "ignore")
somatic_unparsed.name = "somatic MAF unparsed"
somatic_maf_u = DataFrameLoader(data_directory + "somatic.maf.gz").createDataFrame()
patient_ids = create_patient_ids(clinical_unfiltered) #maps C3L-**** number to S*** number
somatic_maf_u = link_patient_ids(patient_ids, somatic_maf_u) #adds S*** number to somatic mutations dataframe
somatic_maf_u = somatic_maf_u.set_index("Clinical_Patient_Key")
somatic_maf = somatic_maf_u.drop(casesToDrop, errors = "ignore")
somatic_maf = somatic_maf.reset_index()
somatic_maf.name = "somatic MAF"


#metaData = MetaData(clinical)
#molecularData = MolecularData(proteomics, transcriptome, cna, phosphoproteomics)
warning()
def list_data():
    """
    Parameters
    None

    Prints list of loaded data frames and dimensions

    Returns
    None
    """
    print("Below are the available endometrial data frames contained in this package:")
    data = [clinical, derived_molecular, acetylproteomics, proteomics, transcriptomics, transcriptomics_circular, miRNA, cna, phosphoproteomics, phosphoproteomics_gene, somatic_binary, somatic_maf]
    for dataframe in data:
        print("\t", dataframe.name)
        print("\t", "\t", "Dimensions:", dataframe.shape)
    #print("To find how to access the data, view the documentation with either list_api() or visit the github page with help().")
def list_api():
    """
    Parameters
    None

    Prints docstrings for all accessible functions

    Returns
    None
    """
    help(__name__)
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
def unfiltered_warning():
    """
    Parameters
    None

    Prints warning to about the unfiltered data

    Returns
    None
    """

    message = "IMPORTANT! Data has been filtered due to quality check on samples. Inclusion of unfiltered samples in analyses is NOT recommended."
    print(message)
def get_clinical(unfiltered=False):
    """
    Parameters
    unfiltered: boolean indicating whether to return unfiltered clinical data, aka clinical["Case_excluded"] == "Yes"

    Returns
    Clinical dataframe
    """
    if unfiltered:
        unfiltered_warning()
        return clinical_unfiltered
    return clinical
def get_derived_molecular(unfiltered=False):
    """
    Parameters
    unfiltered: boolean indicating whether to return unfiltered derived molecular data

    Returns
    Derived Molecular dataframe
    """
    if unfiltered:
        unfiltered_warning()
        return derived_molecular_u
    return derived_molecular
def get_acetylproteomics(unfiltered=False):
    """
    Parameters
    unfiltered: boolean indicating whether to return unfiltered acetylproteomics data

    Returns
    Acetylproteomics dataframe
    """
    if unfiltered:
        unfiltered_warning()
        return acetylproteomics_u
    return acetylproteomics
def get_proteomics(unfiltered=False):
    """
    Parameters
    unfiltered: boolean indicating whether to return unfiltered proteomics data

    Returns
    Proteomics dataframe
    """
    if unfiltered:
        unfiltered_warning()
        return proteomics_u
    return proteomics
def get_transcriptomics(data_type="linear", unfiltered=False):
    """
    Parameters
    data_type: either "linear", "circular", or "miRNA". Indicates which transcriptomics dataframe you want.
    unfiltered: boolean indicating whether to return unfiltered transcriptomics data

    Returns
    Transcriptomics dataframe
    """
    if data_type == "linear":
        if unfiltered:
            unfiltered_warning()
            return transcriptomics_u
        return transcriptomics
    elif data_type == "circular":
        if unfiltered:
            unfiltered_warning()
            return transcriptomics_circular_u
        return transcriptomics_circular
    elif data_type == "miRNA":
        if unfiltered:
            unfiltered_warning()
            return miRNA_u
        return miRNA
    else:
        raise ValueError("Invalid value for get_transcriptomics() data_type parameter.\n\tYou passed: '{}'\n\tOptions: 'linear', 'circular', or 'miRNA'".format(data_type))

def get_transcriptomics_linear(unfiltered=False):
    """Gets transcriptomics_linear dataframe.

    Parameters:
    unfiltered (bool, optional): Whether to include unfiltered samples. Default is false.

    Returns:
    pandas.core.frame.DataFrame: The transcriptomics_linear dataframe.
    """
    if unfiltered:
        unfiltered_warning()
        return transcriptomics_u
    return transcriptomics

def get_transcriptomics_circular(unfiltered=False):
    """Gets transcriptomics_circular dataframe.

    Parameters:
    unfiltered (bool, optional): Whether to include unfiltered samples. Default is false.

    Returns:
    pandas.core.frame.DataFrame: The transcriptomics_circular dataframe.
    """
    if unfiltered:
        unfiltered_warning()
        return transcriptomics_circular_u
    return transcriptomics_circular

def get_miRNA(unfiltered=False):
    """Gets miRNA dataframe.

    Parameters:
    unfiltered (bool, optional): Whether to include unfiltered samples. Default is false.

    Returns:
    pandas.core.frame.DataFrame: The miRNA dataframe.
    """
    if unfiltered:
        unfiltered_warning()
        return miRNA_u
    return miRNA

def get_cna(unfiltered=False):
    """
    Parameters
    unfiltered: boolean indicating whether to return unfiltered cna data

    Returns
    cna dataframe
    """
    if unfiltered:
        unfiltered_warning()
        return cna_u
    return cna
def get_phosphoproteomics(gene_level=False, unfiltered=False):
    """
    Parameters
    gene_level: boolean indicating whether to return gene level phosphoproteomics (returns site level if false)
    unfiltered: boolean indicating whether to return unfiltered phosphoproteomics data

    Returns
    Phosphoproteomics dataframe
    """
    if gene_level:
        if unfiltered:
            unfiltered_warning()
            return phosphoproteomics_gene_u
        return phosphoproteomics_gene
    if unfiltered:
        unfiltered_warning()
        return phosphoproteomics_u
    return phosphoproteomics

def get_phosphoproteomics_site(unfiltered=False):
    """Gets the phosphoproteomics_site dataframe.

    Parameters:
    unfiltered (bool, optional): Whether to include unfiltered samples. Default is false.

    Returns:
    pandas.core.frame.DataFrame: The phosphoproteomics_site dataframe.
    """
    if unfiltered:
        unfiltered_warning()
        return phosphoproteomics_u
    return phosphoproteomics

def get_phosphoproteomics_gene(unfiltered=False):
    """Gets the phosphoproteomics_gene dataframe.

    Parameters:
    unfiltered (bool, optional): Whether to include unfiltered samples. Default is false.

    Returns:
    pandas.core.frame.DataFrame: The phosphoproteomics_gene dataframe.
    """
    if unfiltered:
        unfiltered_warning()
        return phosphoproteomics_gene_u
    return phosphoproteomics_gene

def get_phosphosites(genes):
    """Returns dataframe with all phosphosites of specified gene or list of genes.

    Parameters:
    genes (str or list): gene or list of genes to use to select phosphosites. str if single, list if multiple.

    Returns:
    pandas.core.frame.DataFrame: The phosphoproteomics for the specified genes.
    """
    return Utilities().get_omics_from_str_or_list(phosphoproteomics, genes)

def get_mutations(binary=False, unparsed=False, unfiltered=False):
    """
    Parameters
    binary: boolean indicating whether to retrieve the somatic mutations binary data
    unparsed: boolean indicating whether to retrieve unparsed somatic mutations maf data
    unfiltered: boolean indicating whether to return unfiltered somatic data

    Default behavior is to return parsed somatic mutations maf data

    Returns
    Somatic mutations dataframe corresponding with parameters provided
    """
    if binary:
        if unfiltered:
            unfiltered_warning()
            return somatic_binary_u
        return somatic_binary
    if unparsed:
        if unfiltered:
            unfiltered_warning()
            return somatic_unparsed_u
        return somatic_unparsed
    if unfiltered:
        unfiltered_warning()
        return somatic_maf_u
    return somatic_maf

def get_mutations_maf(unfiltered=False):
    """Gets the somatic_maf mutations dataframe.

    Parameters:
    unfiltered (bool, optional): Whether to include unfiltered samples. Default is false.

    Returns:
    pandas.core.frame.DataFrame: The somatic_maf mutations dataframe.
    """
    if unfiltered:
        unfiltered_warning()
        return somatic_maf_u
    return somatic_maf

def get_mutations_binary(unfiltered=False):
    """Gets the somatic_binary mutations dataframe.

    Parameters:
    unfiltered (bool, optional): Whether to include unfiltered samples. Default is false.

    Returns:
    pandas.core.frame.DataFrame: The somatic_binary mutations dataframe.
    """
    if unfiltered:
        unfiltered_warning()
        return somatic_binary_u
    return somatic_binary

def get_mutations_unparsed(unfiltered=False):
    """Gets the somatic_unparsed mutations dataframe.

    Parameters:
    unfiltered (bool, optional): Whether to include unfiltered samples. Default is false.

    Returns:
    pandas.core.frame.DataFrame: The somatic_unparsed mutations dataframe.
    """
    if unfiltered:
        unfiltered_warning()
        return somatic_unparsed_u
    return somatic_unparsed

def get_clinical_cols():
    """
    Parameters
    None

    Returns
    List of clincal dataframe columns, aka data types (i.e. BMI, Diabetes)
    """
    return clinical.columns
def get_derived_molecular_cols():
    """
    Parameters
    None

    Returns
    List of derived molecular dataframe columns
    """
    return derived_molecular.columns
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
    columns: single column name or array of column names to select for in the cna dataframe

    Returns
    Dataframe of specified columns (or Series if one column) of cna data
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
    """
    Under construction
    """
    #"""
    #Parameters
    #tumor_type is the tumor type, e.g. colon
    #query_type is the type of tumor query, e.g. by SNP, mutated gene, outlier
    #value corresponds with the query type, e.g. TP53 for mutated gene or EGFR for outlier

    #Returns

    #"""
    dataframe = None #TODO what should the dataframe be?
    return Queries(dataframe).query(tumor_type, query_type, value)

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
        'cna',
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
        'cna',
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
        'cna',
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
        'cna',
        'phosphoproteomics_site',
        'phosphoproteomics_gene']
    if (omics_df.name not in valid_dfs):
        print("{} is not a valid dataframe for omics_df parameter. Valid options:".format(omics_df.name))
        for df_name in valid_dfs:
            print('\t' + df_name)
        return

    # Return the merge.
    return Utilities().append_mutations_to_omics(somatic_maf, omics_df, mutation_genes, omics_genes, multiple_mutations, show_location)

def git_help():
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
def version():
    """
    Parameters
    None

    Prints version number of CPTAC package

    Returns
    Version number
    """
    version = {}
    with open(dir_path + os.sep + ".." + os.sep + "version.py") as fp: #.. required to navigate up to CPTAC folder from Endometrial folder, TODO: how to navigate from dataTest.py?
    	exec(fp.read(), version)
    return(version['__version__'])
