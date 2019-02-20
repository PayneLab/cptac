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
            df = df.transpose()
            df = df.sort_index()
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
            df.name = self.name
            return df
        elif self.name.split("_")[0] == "somatic":
            df = pd.read_csv(self.fileName, sep = "\t")
            if "Tumor_Sample_Barcode" in df.columns:
                split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n = 1, expand = True)
                df["Tumor_Sample_Barcode"] = split_barcode[0]
            parsedDf = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
            parsedDf = parsedDf.rename({"Tumor_Sample_Barcode":"Patient_Id","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
            parsedDf = parsedDf.set_index("Patient_Id")
            parsedDf.name = self.name
            return parsedDf
        else:
            print("Error reading", self.fileName)
