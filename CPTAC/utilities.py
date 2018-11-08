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
    def compare_mutations(self, df1, somatic, gene):
        if gene in df1.columns:
            df1_gene = df1[[gene]]
            if sum(somatic["Gene"] == gene) > 0:
                somatic_gene = somatic[somatic["Gene"] == gene]
                somatic_gene = somatic_gene.drop(columns = ["Gene"])
                somatic_gene = somatic_gene.set_index("Clinical_Patient_Key")
                merge = df1_gene.join(somatic_gene, how = "left")
                merge = merge.fillna(value = {'Mutation':"Wildtype"})
                return merge
            else:
                print("Gene", gene, "not found in somatic mutations.")
        else:
            print("Gene", gene, "not found in", df1.name, "data")
    def compare_mutations_trans(self, df1, df1Gene, somatic, somaticGene):
        if df1Gene in df1.columns:
            df1_gene = df1[[df1Gene]]
            if sum(somatic["Gene"] == somaticGene) > 0:
                somatic_gene = somatic[somatic["Gene"] == somaticGene]
                somatic_gene = somatic_gene.drop(columns = ["Gene"])
                somatic_gene = somatic_gene.set_index("Clinical_Patient_Key")
                merge = df1_gene.join(somatic_gene, how = "left")
                merge = merge.fillna(value = {'Mutation':"Wildtype"})
                return merge
            else:
                print("Gene", somaticGene, "not found in somatic mutations.")
        else:
            print("Gene", df1Gene, "not found in", df1.name,"data")
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
    def compare_phosphosites(self, proteomics, phosphoproteomics, gene):
        """
        Returns dataframe with a column from proteomics for the gene specified,
        as well as columns for all phosphoproteomics columns beginning with the specified gene
        """
        if gene in proteomics.columns:
            df = proteomics[[gene]]
            regex = gene + ".*"
            phosphosites = phosphoproteomics.filter(regex =(regex))
            df = df.add(phosphosites, fill_value=0)
            df.name = gene + " proteomics and phosphoproteomics"
            return df
        else:
            print(gene, "not found in proteomics dataframe. Available genes can be checked by entering CPTAC.get_proteomics().columns")
