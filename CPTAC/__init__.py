import os
import webbrowser
from .dataframe import DataFrameLoader
from .meta import MetaData
from .molecular import MolecularData
from .utilities import Utilities
from .queries import Queries
#temporary fix for .N .T issue
def warning():
    print("\n","******PLEASE READ******")
    print("WARNING: This data is under a publication embargo until July 1, 2019. CPTAC is a community resource project and data are made available rapidly after generation for community research use. The embargo allows exploring and utilizing the data, but the data may not be in a publication until July 1, 2019. Please see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter CPTAC.embargo() to open the webpage for more details.")

"""
Executes on import CPTAC statement. Selects files from docs folder in CPTAC package
utilizing DataFrameLoader from dataframe.py. Prints update as files are loaded into
dataframes.
"""
dir_path = os.path.dirname(os.path.realpath(__file__))
data_directory = dir_path + os.sep + "Data" + os.sep

print("Loading Clinical Data...")
clinical = DataFrameLoader(data_directory + "clinical.txt.gz").createDataFrame()

print("Loading Proteomics Data...")
proteomics = DataFrameLoader(data_directory + "proteomics.cct.gz").createDataFrame()

print("Loading Transcriptomics Data...")
transcriptomics = DataFrameLoader(data_directory + "transcriptomics.cct.gz").createDataFrame()

print("Loading CNA Data...")
cna = DataFrameLoader(data_directory + "CNA.cct.gz").createDataFrame()

print("Loading Phosphoproteomics Data...")
phosphoproteomics = DataFrameLoader(data_directory + "phosphoproteomics.cct.gz").createDataFrame()

print("Loading Somatic Data...")
somatic = DataFrameLoader(data_directory + "somatic.cbt.gz").createDataFrame()

#metaData = MetaData(clinical)#, clinical_meta)
#molecularData = MolecularData(proteomics, transcriptome, cna, phosphoproteomics)
warning()
def list():
    """
    Prints list of loaded data frames and dimensions
    """
    print("Below are the available data frames contained in this package:")
    print("\t", clinical.name)
    print("\t","\t", "Dimensions:", clinical.shape)
    print("\t", proteomics.name)
    print("\t","\t", "Dimensions:", proteomics.shape)
    print("\t", transcriptomics.name)
    print("\t","\t", "Dimensions:", transcriptomics.shape)
    print("\t", cna.name)
    print("\t","\t", "Dimensions:", cna.shape)
    print("\t", phosphoproteomics.name)
    print("\t","\t", "Dimensions:", phosphoproteomics.shape)
    print("\t", somatic.name)
    print("\t","\t", "Dimensions:", somatic.shape)
    print("To access the data, use a get function with the data frame name, i.e. CTPAC.get_proteomics()")
def get_clinical():
    """Returns clincal dataframe"""
    return clinical
def get_proteomics():
    """Returns proteomics dataframe"""
    return proteomics
def get_transcriptomics():
    """Returns transcriptomics dataframe"""
    return transcriptomics
def get_CNA():
    """Returns CNA dataframe"""
    return cna
def get_phosphoproteomics():
    """Returns phosphoproteomics dataframe"""
    return phosphoproteomics
def get_somatic():
    """Returns somatic mutations dataframe"""
    return somatic
def get_meta_cols():
    """
    Returns list of clincal dataframe columns,
    aka data types (i.e. BMI, Diabetes)
    """
    return clinical.columns
def get_proteomics_cols():
    """Returns list of columns of proteomics dataframe"""
    return proteomics.columns
def get_transcriptomics_cols():
    """Returns list of columns of transcriptomics dataframe"""
    return transcriptomics.columns
def get_cohort_meta(cols):
    """Returns specified column or columns of clinical data"""
    return clinical[cols]
def get_proteomics_quant(colon_ids):
    """Returns specified row or rows of proteomics data"""
    return proteomics.loc(colon_ids)
def get_cohort_proteomics(cols):
    """Returns specified column or columns of proteomics data"""
    return proteomics[cols]
def get_cohort_transcriptomics(cols):
    """Returns specified column or columns of transcriptomics data"""
    return transcriptomics[cols]
def get_cohort_cna(cols):
    """Returns specified column or columns of CNA data"""
    return cna[cols]
def get_cohort_phosphoproteomics(cols):
    """Returns specified column or columns of phosphoproteomics data"""
    return phosphoproteomics[cols]

def get_tumor_ids(tumor_type, query_type, value):
    """
    tumor_type is the tumor type, e.g. colon
    query_type is the type of tumor query, e.g. by SNP, mutated gene, outlier
    value corresponds with the query type, e.g. TP53 for mutated gene or EGFR for outlier
    """
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
    Returns dataframe containing two columns. Each column is the data for the
    specified gene from the two specified dataframes
    """
    if isinstance(gene, str):
        return Utilities().compare_gene(df1, df2, gene)
    else:
        return Utilities().compare_genes(df1, df2, gene)
def compare_clinical(clinical, data, clinical_col):
    """
    Returns dataframe with specified column from clinical dataframe added to
    specified dataframe (i.e., proteomics) for comparison and easy plotting
    """
    #TODO: do we need clinical parameter? Could just grab it from loaded data?
    return Utilities().compare_clinical(clinical, data, clinical_col)
def compare_phosphosites(gene):
    return Utilities().compare_phosphosites(proteomics, phosphoproteomics, gene)

def help():
    """Opens github help page"""
    print("Opening help.txt in web browser...")
    webbrowser.open("https://github.com/PayneLab/CPTAC/blob/master/doc/help.txt")
def embargo():
    print("Opening embargo details in web browser...")
    webbrowser.open("https://proteomics.cancer.gov/data-portal/about/data-use-agreement")
def start():
    print("Welcome to our CPTAC data. Enter CPTAC.help() to open our Github help page.")
