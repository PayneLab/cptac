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
class Utilities:

    def __init__(self):
        pass
    def get_gene_mapping(self):
        print("Under construction")
    def convert(self, snp_or_sap):
        print("Under construction")
    def add_mutation_hierarchy(self, somatic):
        """
        Parameters
        somatic: somatic data to add mutation hierarchy to

        Retunrs
        Somatic mutation dataframe with added mutation hierarchy
        """
        mutation_hierarchy = {"Missense_Mutation":0,"In_Frame_Del":0,"In_Frame_Ins":0,"Splice_Site":1,"Frame_Shift_Ins":1,"Nonsense_Mutation":1,"Frame_Shift_Del":1,"Nonstop_Mutation":1}
        hierarchy = []
        for x in somatic["Mutation"]:
            if x in mutation_hierarchy.keys():
                hierarchy.append(mutation_hierarchy[x])
            else:
                hierarchy.append(float('NaN'))
        somatic["Mutation_Hierarchy"] = hierarchy
        return somatic
    def merge_somatic(self, somatic, gene, df_gene, multiple_mutations = False): #private
        """
        Parameters
        somatic: somatic mutations dataframe that will be parsed for specified gene data
        gene: string of gene to be selected for in somatic mutation data
        df_gene: selection of omics data for particular gene to be merged with somatic data
        multiple_mutations: boolean indicating whether to include multiple mutations for specified gene in an individual

        Returns
        Dataframe of merged somatic and omics dataframes based on gene provided
        """
        if sum(somatic["Gene"] == gene) > 0:
            somatic_gene = somatic[somatic["Gene"] == gene]
            somatic_gene = somatic_gene.drop(columns = ["Gene"])
            somatic_gene = somatic_gene.set_index("Clinical_Patient_Key")
            if not multiple_mutations:
                somatic_gene = self.add_mutation_hierarchy(somatic_gene)
                somatic_gene = somatic_gene.sort_values(by = ["Clinical_Patient_Key","Mutation_Hierarchy"], ascending = [True,False])
                somatic_gene = somatic_gene[~somatic_gene.index.duplicated(keep="first")]
            merge = df_gene.join(somatic_gene, how = "left")
            merge = merge.fillna(value = {'Mutation':"Wildtype"})
            merge["index"] = merge.index
            merge["Patient_Type"] = np.where(merge.index <= "S100", "Tumor", "Normal")
            merge.name = df_gene.columns[0] + " omics data with " + gene + " mutation data"
            return merge
        else:
            print("Gene", gene, "not found in somatic mutations.")
    def get_phosphosites(self, phosphoproteomics, gene):
        regex = gene + ".*"
        phosphosites = phosphoproteomics.filter(regex = (regex))
        if len(phosphosites.columns) == 0:
            print("Gene",gene, "not found in phosphoproteomics data")
        else:
            return phosphosites
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
        Dataframe containing columns equally to the number of genes provided times two. Each two-column set is the data for each specified gene from the two specified dataframes
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
    def merge_mutations(self, omics, somatic, gene, duplicates = False):
        """
        Parameters
        omics: dataframe containing specific omics data
        somatic: dataframe of somatic mutation data
        gene: string of specific gene to merge omics and somatic data on
        duplicates: boolean value indicating whether to include duplicate gene mutations for an individual

        Returns
        Dataframe of merged omics and somatic data based on gene provided
        """
        if gene in omics.columns:
            omics_gene_df = omics[[gene]]
            if duplicates:
                return self.merge_somatic(somatic, gene, omics_gene_df, multiple_mutations = True)
            else:
                return self.merge_somatic(somatic, gene, omics_gene_df)[[gene, "Mutation", "Patient_Type"]]
        elif omics.name.split("_")[0] == "phosphoproteomics":
            phosphosites = self.get_phosphosites(omics, gene)
            if len(phosphosites.columns) > 0:
                if duplicates:
                    return self.merge_somatic(somatic, gene, phosphosites, multiple_mutations = True)
                else:
                    columns = list(phosphosites.columns)
                    columns.append("Mutation")
                    columns.append("Patient_Type")
                    merged_somatic = self.merge_somatic(somatic, gene, phosphosites)
                    return merged_somatic[columns]

        else:
            print("Gene", gene, "not found in", omics.name, "data")
    def merge_mutations_trans(self, omics, omicsGene, somatic, somaticGene, duplicates = False):
        """
        Parameters
        omics: dataframe containing specific omics data (i.e. proteomics, transcriptomics)
        omicsGene: string of specific gene to merge from omics data
        somatic: dataframe of somatic mutation data
        somaticGene: string of specific gene to merge from somatic data
        duplicates: boolean value indicating whether to include duplicate gene mutations for an individual

        Returns
        Dataframe of merged omics data (based on specific omicsGene) with somatic data (based on specific somaticGene)
        """
        if omicsGene in omics.columns:
            omics_gene_df = omics[[omicsGene]]
            if duplicates:
                return self.merge_somatic(somatic, somaticGene, omics_gene_df, multiple_mutations = True)
            else:
                return self.merge_somatic(somatic, somaticGene, omics_gene_df)[[omicsGene, "Mutation", "Patient_Type"]]
        elif omics.name.split("_")[0] == "phosphoproteomics":
            phosphosites = self.get_phosphosites(omics, omicsGene)
            if len(phosphosites.columns) > 0:
                if duplicates:
                    return self.merge_somatic(somatic, somaticGene, phosphosites, multiple_mutations = True)
                else:
                    columns = list(phosphosites.columns)
                    columns.append("Mutation")
                    columns.append("Patient_Type")
                    merged_somatic = self.merge_somatic(somatic, somaticGene, phosphosites)
                    return merged_somatic[columns]
        else:
            print("Gene", omicsGene, "not found in", omics.name,"data")
    def compare_clinical(self, clinical, data, clinical_col):
        """
        Parameters
        clinical: clinical dataframe for omics data to be appended with
        data: omics data for clinical data to be appended with
        clinical_col: column in clinical dataframe to be inserted into provided omics data

        Returns
        Dataframe with specified column from clinical dataframe added to specified dataframe (i.e., proteomics) for comparison and easy plotting
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
        Parameters
        gene: proteomics gene to query phosphoproteomics dataframe

        Searches for any phosphosites on the gene provided

        Returns
        Dataframe with a column from proteomics for the gene specified, as well as columns for all phosphoproteomics columns beginning with the specified gene
        """
        if gene in proteomics.columns:
            df = proteomics[[gene]]
            phosphosites = self.get_phosphosites(phosphoproteomics, gene)
            if len(phosphosites.columns) > 0:
                df = df.add(phosphosites, fill_value=0)
                df.name = gene + " proteomics and phosphoproteomics"
                return df
        else:
            print(gene, "not found in proteomics dataframe. Available genes can be checked by entering CPTAC.get_proteomics().columns")
