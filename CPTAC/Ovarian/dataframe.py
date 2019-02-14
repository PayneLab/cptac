import pandas as pd
import numpy as np
import os

class DataFrameLoader:
    def __init__(self,fileName):
        self.fileName = fileName
        f = fileName.split(os.sep)
        f = f[len(f) - 1]
        self.name = f.split(".")[0]
    def clip_indices(self, df, clip_number): #private
        """
        Parameters
        df: dataframe to be clipped
        clip_number: row index for first cutoff value. Results in df.iloc[0:clip_number], therefore leaving off the index number provided.

        Selects the first clip_number rows and clips the first character off of the indices, for dataframes with C***** and N*****.
        Possibly a better solution

        Returns
        Dataframe with amended and selected indices.
        """
        c_index = df.index[0:clip_number].str[1:]
        c_df_indices = df.iloc[0:clip_number].index #these next two lines of code eliminate SetWithCopy warning
        c_df = df.loc[c_df_indices,:] #The original c_df = df.iloc[0:c_number], then c_df["hgnc"] = c_index gave SetWithCopy warning
        c_df["hgnc"] = c_index
        c_df = c_df.set_index("hgnc")
        return c_df
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
            #c_df = self.clip_indices(df, clip_number=83) #83rd row is where the N***** indices start, is there a better solution??
            #c_df.name = self.name
            #return c_df

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
            #c_df = self.clip_indices(df, clip_number=83)
            #c_df.name = self.name
            #return c_df
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
