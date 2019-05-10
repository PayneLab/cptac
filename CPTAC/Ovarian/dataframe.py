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

import pandas as pd
import numpy as np
import os

class DataFrameLoader:
    def __init__(self,fileName):
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
        print("Loading",self.name,"data...")
        if self.name == "proteomics":
            df = pd.read_csv(self.fileName,sep="\t", index_col = 0)
            df = df[df["hgnc_symbol"].notnull()] #drops all nan values in hgnc_symbol column
            df = df.set_index("hgnc_symbol")
            df = df.sort_index()
            df = df.transpose()
            c_index = df.index[0:83].str[1:] #drops letter off all indices with "C"
            index = c_index.append(df.index[83:])
            df = df.set_index(index)
            idx = df.index.values.tolist()
            idx_to_drop = [id for id in idx if id.startswith('OV_QC')]
            df = df.drop(idx_to_drop) # Drop all OV_QC* samples--they're quality control samples not relevant for data analysis
            df.name = self.name
            return df
        elif self.name == "clinical":
            df = pd.read_csv(self.fileName, sep="\t")
            df = df.set_index("PPID")
            df = df[~df.index.duplicated(keep="first")]
            df.name = self.name
            return df
        elif self.name == "phosphoproteomics":
            df = pd.read_csv(self.fileName, sep = "\t",index_col = 0)
            df = df[df["site"].notnull()] #drops all nan values in site column
            df = df.drop(["refseq_peptide","Peptide"],axis=1)
            df = df.set_index("site")
            df = df.sort_index()
            df = df.transpose()
            c_index = df.index[0:83].str[1:] #drops letter off all indices with "C"
            index = c_index.append(df.index[83:])
            df = df.set_index(index)
            idx = df.index.values.tolist()
            idx_to_drop = [id for id in idx if id.startswith('OV_QC')]
            df = df.drop(idx_to_drop) # Drop all OV_QC* samples--they're quality control samples not relevant for data analysis
            df.name = self.name
            return df
        elif self.name == "transcriptomics":
            df = pd.read_csv(self.fileName, sep="\t", index_col=0)
            df = df.sort_index()
            df = df.transpose()
            df = df.sort_index()
            df = df.drop(columns = df.columns[0:23])#drop all date values until new data is uploaded
            df.name = self.name
            return df
        elif self.name == "cnv":
            df = pd.read_csv(self.fileName, sep="\t", index_col=0)
            df = df.sort_index()
            df = df.transpose()
            df = df.sort_index()
            df.name = "CNV"
            return df
        elif self.name == "somatic_38":
            df = pd.read_csv(self.fileName, sep = "\t")
            if "Tumor_Sample_Barcode" in df.columns:
                split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n = 1, expand = True)
                df["Tumor_Sample_Barcode"] = split_barcode[0]
            parsedDf = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
            parsedDf = parsedDf.rename({"Tumor_Sample_Barcode":"Patient_Id","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
            parsedDf = parsedDf.set_index("Patient_Id")
            parsedDf.name = 'somatic_mutation'
            return parsedDf
        else:
            print("Error reading", self.fileName)
