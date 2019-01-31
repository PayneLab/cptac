import pandas as pd
import numpy as np

class Utilities:
    def __init__(self):
        pass
    def compare_gene(self, df1, df2, gene):
        """
        Parameters
        df1: omics dataframe (proteomics) to be selected from
        df2: other omics dataframe (transcriptomics) to be selected from
        gene: gene to select from each of the dataframes

        Returns
        Dataframe containing two columns. Each column is the data for the specified gene from the two specified dataframes
        """
        if gene in df1.columns and gene in df2.columns:
            common = df1.index.intersection(df2.index)
            df1Matched = df1.loc[common]
            df1Matched = df1Matched.sort_index()
            df2Matched = df2.loc[common]
            df2Matched = df2Matched.sort_index()
            #TODO how to check if df1.name without crashing?
            dict = {df1.name:df1Matched[gene], df2.name:df2Matched[gene]}
            df = pd.DataFrame(dict, index = df1Matched.index)
            df.name = gene
            return df
        else:
            print(gene,"not found in provided dataframes. Please check that the specified gene is included in the provided dataframes")
    def compare_genes(self, df1, df2, genes):
        """
        Parameters
        df1: omics dataframe (proteomics) to be selected from
        df2: other omics dataframe (transcriptomics) to be selected from
        genes: gene or array of genes to select from each of the dataframes

        Returns
        Dataframe containing columns equal to the number of genes provided times two. Each two-column set is the data for each specified gene from the two specified dataframes
        """
        dfs = pd.DataFrame(index = df1.index.intersection(df2.index))
        for gene in genes:
            df = Utilities().compare_gene(df1, df2, gene)
            new_col1 = df1.name + "_" + gene
            new_col2 = df2.name + "_" + gene
            df = df.rename(columns = {df1.name:new_col1, df2.name:new_col2})
            dfs = dfs.add(df, fill_value=0)
        dfs.name = str(len(genes)) + " Genes Combined"
        return dfs
