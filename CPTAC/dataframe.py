import numpy as np
import pandas as pd
import os
from .fileLoader import FileLoader
class DataFrameLoader:
    def __init__(self, fileName):
        self.fileName = fileName
    def createDataFrame(self):
        file = FileLoader(self.fileName).readFile()
        if self.fileName.endswith('.csv'):
            df = pd.read_csv(file, index_col=0)
            df = df.iloc[1:]
            #TODO how to preserve type?
            df = df.apply(pd.to_numeric, errors='coerce')
            return df
        elif self.fileName.endswith('.txt'):
            line = file.readline()
            line = line.split()
            rows = line[1:] #C3L-00358 etc.
            line = file.readline()
            dict = {}
            while line:
                line = line.split()
                floats = []
                for num in line[1:]:
                    if num != 'NA':
                        floats.append(float(num))
                    else:
                        floats.append(None)

                dict.update({line[0]:floats})
                line = file.readline()
            df = pd.DataFrame(dict, rows)
            df = df.sort_index()
            f = self.fileName.split(os.sep)
            f = f[len(f) - 1]
            df.name = f.split(".")[0]
                    #print(df.head())
            return df
        else:
            print("Error reading file")





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
