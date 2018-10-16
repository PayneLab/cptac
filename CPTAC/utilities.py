import pandas as pd
class Utilities:
    def __init__(self):
        pass
    def get_gene_mapping(self):
        print("Under construction")
    def convert(self, snp_or_sap):
        print("Under construction")
    def compare_gene(self, df1, df2, gene):
        """
        Returns dataframe containing two columns. Each column is the data for the
        specified gene from the two specified dataframes
        """
        if gene in df1.columns and gene in df2.columns:
            gene = gene.upper()
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
        Returns dataframe of two column sets corresponding with the provided
        array of genes
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
    def compare_clinical(self, clinical, data, clinical_col):
        """
        Returns dataframe with specified column from clinical dataframe added to
        specified dataframe (i.e., proteomics) for comparison and easy plotting
        """
        if clinical_col in clinical:
            df = data[data.columns]
            df.insert(0, clinical_col, clinical[clinical_col])
            df.name = data.name + " with " + clinical_col
            return df
        else:
            print(clinical_col, "not found in clinical dataframe. You can check the available columns by entering CPTAC.get_clincal().columns")
