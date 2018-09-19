import os
from .dataframe import DataFrameLoader

def get_meta_completeness():
    meta = DataFrameLoader("CPTAC" + os.sep + "Data" + os.sep + "clinical.csv").createDataFrame()
    return meta

def start():

    joke = u'Wenn ist das Nunst\u00fcck git und Slotermeyer? Ja! ... Beiherhund das Oder die Flipperwaldt gersput.'

    print("Welcome to our CPTAC data. Below are a list of main functions for this package: \nfunction1()\nfunction2()\nfunction3()")
