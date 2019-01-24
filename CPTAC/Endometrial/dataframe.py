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
        Parameters
        None

        Returns
        Dataframe of parsed datafile depending on the data type
        """
        if bool(re.search(r'\.txt[.|(a-z)]{,7}$', self.fileName)):
            #temp fix for reading error on clinical_v2:
            file = open(self.fileName, "r", errors="ignore")
            df = pd.read_csv(file, sep="\t", index_col=0)
            df = df.sort_index()

            f = self.fileName.split(os.sep)
            f = f[len(f) - 1]
            df.name = f.split(".")[0]
            return df
        elif bool(re.search(r'\.(cct|cbt)[.|(a-z)]{,7}$', self.fileName)):
            df = pd.read_csv(self.fileName, sep="\t", index_col=0)
            df = df.transpose()
            df = df.sort_index()
            f = self.fileName.split(os.sep)
            f = f[len(f) - 1]
            df.name = f.split(".")[0]
            return df
        elif bool(re.search(r'\.maf[.|(a-z)]{,7}$', self.fileName)):
            df = pd.read_csv(self.fileName, sep = "\t")
            if "Tumor_Sample_Barcode" in df.columns:
                split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n = 1, expand = True)
                df["Tumor_Sample_Barcode"] = split_barcode[0]
            parsedDf = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
            parsedDf = parsedDf.rename({"Tumor_Sample_Barcode":"Patient_Id","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
            f = self.fileName.split(os.sep)
            f = f[len(f) - 1]
            parsedDf.name = f.split(".")[0] + " MAF"
            return parsedDf
        else:
            print("Error reading", self.fileName)
