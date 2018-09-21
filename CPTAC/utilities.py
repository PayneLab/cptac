import pandas as pd
class Utilities:
    def __init__(self):
        pass
    def get_gene_mapping(self):
        pass
    def convert(self, snp_or_sap):
        pass
    def compare_gene(self, df1, df2, gene):
        common = df1.index.intersection(df2.index)
        df1Matched = df1.loc[common]
        df1Matched.sort_index()
        df2Matched = df2.loc[common]
        df2Matched.sort_index()
        #TODO how to check if df1.name without crashing?
        dict = {df1.name:df1Matched[gene], df2.name:df2Matched[gene]}
        df = pd.DataFrame(dict, index = df1Matched.index)
        return df
