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

class DataFrameLoader:
    def __init__(self, fileName):
        self.fileName = fileName
        f = fileName.split(os.sep)
        f = f[len(f) - 1]
        self.name = f.split(".")[0] #sets name as first section of file name (i.e. proteomics.txt.gz becomes proteomics)
    def createDataFrame(self):
        """
        Parameters
        None

        Returns
        Dataframe of parsed datafile depending on the data type
        """
        df = None
        print("Loading",self.name,"data...")
        if self.name == "clinical" or self.name.split("_")[0] == "proteomics" or self.name.split("_")[0] == "transcriptomics":
            df = pd.read_csv(self.fileName, sep="\t",index_col=0)
            df = df.transpose()
            df.name = self.name
        elif self.name.split("_")[0] == "phosphoproteomics": #column names are SLC12A8_T485_A0AV02:T485, should this be changed?
            df = pd.read_csv(self.fileName, sep="\t",index_col=0)
            df = df.transpose()
            df.name = self.name
        elif self.name == "miRNA": #column names are hsa-let-7a-2-3p, should this be changed?
            df = pd.read_csv(self.fileName, sep="\t",index_col=0)
            df = df.transpose()
            df.name = self.name
        elif self.name.split("_")[0] == "mutation":
            if self.fileName.split(os.sep)[-1].split(".")[1] == "cbt":
                df = pd.read_csv(self.fileName, sep="\t",index_col=0)
                df = df.transpose()
                df.name = self.name
            elif self.fileName.split(os.sep)[-1].split(".")[1] == "txt":
                df = pd.read_csv(self.fileName, sep="\t")#.set_index("SampleID").sort_index()
                df = df.sort_values(by="SampleID")
                df = df[["SampleID","Gene","Variant_Type","Protein_Change"]]
                df = df.rename({"Variant_Type":"Mutation","Protein_Change":"Location"},axis="columns")
                df.name = "somatic_" + self.name
        else:
            error_message = "Error reading " + self.fileName
            raise IOError(error_message)
            #print("Error reading", self.fileName)
        return df
