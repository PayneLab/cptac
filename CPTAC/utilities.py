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
    def compare_clinical(self, clinical, data, clinical_col):
        """
        Returns dataframe with specified column from clinical dataframe added to
        specified dataframe (i.e., proteomics) for comparison and easy plotting
        """
        common = clinical.index.intersection(data.index)
        clinicalMatched = clinical.loc[common]
        clinicalMatched = clinicalMatched.sort_index()
        dataMatched = data.loc[common]
        dataMatched = dataMatched.sort_index()
        dict = {clinical_col:clinical[clinical_col]}
        for num in range(0, len(dataMatched.columns)):
            column = dataMatched.columns[num]
            dict.update({column:dataMatched[column]})
        df = pd.DataFrame(dict, index = clinicalMatched.index)
        df.name = clinical_col + " with " + data.name
        return df
