import numpy as np
import pandas as pd
import os
import re
import math
from .fileLoader import FileLoader
class DataFrameLoader:
    def __init__(self, fileName):
        self.fileName = fileName
    def createDataFrame(self):
        """
        Creates pandas DataFrame from provided file name. Parses file depending
        on type of file. Assigns name to dataframe based on name of file.
        """
        #checks if file ends with .csv followed by 0 to 7 dots or characters.
        #permits compressed files in various formats
        """if bool(re.search(r'\.csv[.|(a-z)]{,7}$', self.fileName)):
            df = pd.read_csv(self.fileName, index_col=0)
            df = df.iloc[1:]
            #TODO change implementation for excel file with all data in multiple sheets
            f = self.fileName.split(os.sep)
            f = f[len(f) - 1]
            if bool(re.search(r'^clinical\.csv[.|(a-z)]{,7}$', f)):
                df = df.apply(pd.to_numeric, errors='coerce')
            elif bool(re.search(r'^meta_clinical\.csv[.|(a-z)]{,7}$', f)):
                for num in range(0,len(df.index)):
                    if isinstance(df.index[num], str):
                        df = df.rename(index = {df.index[num]:df.index[num].replace(" ","_")})
            df.name = f.split(".")[0]
            return df"""
        if bool(re.search(r'\.txt[.|(a-z)]{,7}$', self.fileName)):
            df = pd.read_csv(self.fileName, sep="\t", index_col=0)
            #df = df.transpose() to put back if .cct doesn't work
            df = df.sort_index()

            f = self.fileName.split(os.sep)
            f = f[len(f) - 1]
            df.name = f.split(".")[0]
            return df
        elif bool(re.search(r'\.cct[.|(a-z)]{,7}$', self.fileName)):
            df = pd.read_csv(self.fileName, sep="\t", index_col=0)
            df = df.transpose()
            df = df.sort_index()
            f = self.fileName.split(os.sep)
            f = f[len(f) - 1]
            df.name = f.split(".")[0]
            return df

        else:
            print("Error reading", self.fileName)





# clinical = {'FIGO': [0,0,0,3],
#         'Diabetes': [0,0,1,0],
#         'BMI': [38.88, 39.76, 51.19, 21.57]}
# df = pd.DataFrame(clinical, index = ['C3L-06', 'C3L-08', 'C3L-32', 'C3L-139'])
# print(df)
# dictionary = {"iphone" : 2007,
# 		"iphone 3G" : 2008,
# 		"iphone 3GS" : 2009,
# 		"iphone 4" : 2010,
# 		"iphone 4S" : 2011,
# 		"iphone 5" : 2012}
# series = pd.Series(dictionary)
# print(series)
