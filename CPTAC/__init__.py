import os
import webbrowser
from .dataframe import DataFrameLoader
from .meta import MetaData
from .molecular import MolecularData
from .utilities import Utilities
from .queries import Queries

dir_path = os.path.dirname(os.path.realpath(__file__))
print("Loading Clinical Data...")
clinical = DataFrameLoader(dir_path + os.sep + "Data" + os.sep + "clinical.csv").createDataFrame()
clinical_meta = DataFrameLoader(dir_path + os.sep + "Data" + os.sep + "meta_clinical.csv").createDataFrame() #TODO isn't finished yet
print("Loading Proteomics Data...")
proteomics = DataFrameLoader(dir_path + os.sep + "Data" + os.sep + "proteomics.txt").createDataFrame()
print("Loading Transcriptome Data...")
transcriptome = DataFrameLoader(dir_path + os.sep + "Data" + os.sep + "transcriptome.txt").createDataFrame()
print("Loading CNA Data...")
cna = DataFrameLoader(dir_path + os.sep + "Data" + os.sep + "CNA.txt").createDataFrame()
print("Loading Phosphoproteomics Data...")
phosphoproteomics = DataFrameLoader(dir_path + os.sep + "Data" + os.sep + "phosphoproteomics.txt").createDataFrame()

metaData = MetaData(clinical, clinical_meta)
molecularData = MolecularData(proteomics, transcriptome, cna, phosphoproteomics)

def list():
    #TODO return list of data frames and dimensions
    print("Below are the available data frames contained in this package:")
    print("\t", clinical.name)
    print("\t","\t", "Dimensions:", clinical.shape)
    print("\t", proteomics.name)
    print("\t","\t", "Dimensions:", proteomics.shape)
    print("\t", transcriptome.name)
    print("\t","\t", "Dimensions:", transcriptome.shape)
    print("\t", cna.name)
    print("\t","\t", "Dimensions:", cna.shape)
    print("\t", phosphoproteomics.name)
    print("\t","\t", "Dimensions:", phosphoproteomics.shape)
    print("To access the data, use a get function with the data frame name, i.e. CTPAC.get_proteomics()")
def get_meta_completeness():
    return clinical
def get_proteomics():
    return proteomics
def get_transcriptome():
    return transcriptome
def get_CNA():
    return cna
def get_phosphoproteomics():
    return phosphoproteomics
def get_meta_cols():
    return clinical.columns
def get_cohort_meta(cols):
    return clinical[cols]
def get_proteomics_quant(colon_ids):
    return proteomics.loc(colon_ids)
def get_proteomics_cols():
    return proteomics.columns
def get_transcriptome_cols():
    return transcriptome.columns
def get_cohort_proteomics(cols):
    return proteomics[cols]
def get_cohort_transcriptome(cols):
    return transcriptome[cols]
def get_cohort_cna(cols):
    return cna[cols]
def get_cohort_phosphoproteomics(cols):
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
    return Utilities().get_gene_mapping()
def convert(snp_or_sap):
    return Utilities().convert(snp_or_sap)
def compare_gene(df1, df2, gene):
    return Utilities().compare_gene(df1, df2, gene)
def compare_clinical(clinical, data, clinical_col):
    return Utilities().compare_clinical(clinical, data, clinical_col)

def help():
    print("Opening help.txt in web browser...")
    webbrowser.open("https://github.com/PayneLab/CPTAC/blob/master/doc/help.txt")

def start():

    joke = u'Wenn ist das Nunst\u00fcck git und Slotermeyer? Ja! ... Beiherhund das Oder die Flipperwaldt gersput.'

    print("Welcome to our CPTAC data. Enter CPTAC.help() to open our Github help page.")
