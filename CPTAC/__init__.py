import os
from .dataframe import DataFrameLoader

clinical = DataFrameLoader("CPTAC" + os.sep + "Data" + os.sep + "clinical.csv").createDataFrame()
proteomics = DataFrameLoader("CPTAC" + os.sep + "Data" + os.sep + "proteomics.txt").createDataFrame()
transcriptome = DataFrameLoader("CPTAC" + os.sep + "Data" + os.sep + "transcriptome.txt").createDataFrame()
cna = DataFrameLoader("CPTAC" + os.sep + "Data" + os.sep + "CNA.txt").createDataFrame()
phosphoproteomics = DataFrameLoader("CPTAC" + os.sep + "Data" + os.sep + "phosphoproteomics.txt").createDataFrame()

def get_meta_completeness():
    return clinical

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

def start():

    joke = u'Wenn ist das Nunst\u00fcck git und Slotermeyer? Ja! ... Beiherhund das Oder die Flipperwaldt gersput.'

    print("Welcome to our CPTAC data. Below are a list of main functions for this package: \nfunction1()\nfunction2()\nfunction3()")
