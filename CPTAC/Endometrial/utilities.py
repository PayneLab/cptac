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
        """
        Under construction
        """
        print("Under construction")
    def convert(self, snp_or_sap):
        """
        Under construction
        """
        print("Under construction")
    def add_mutation_hierarchy(self, somatic): #private
        """
        Parameters
        somatic: somatic data to add mutation hierarchy to

        Returns
        Somatic mutation dataframe with added mutation hierarchy
        """
        mutation_hierarchy = {"Missense_Mutation":0,"In_Frame_Del":0,"In_Frame_Ins":0,"Splice_Site":1,"Frame_Shift_Ins":1,"Nonsense_Mutation":1,"Frame_Shift_Del":1,"Nonstop_Mutation":1}
        hierarchy = []
        for x in somatic["Mutation"]: #for every value in the Mutation column, append its value in the hard coded mutation hierarchy
            if x in mutation_hierarchy.keys():
                hierarchy.append(mutation_hierarchy[x])
            else:
                hierarchy.append(float('NaN'))
        somatic = somatic.assign(Mutation_Hierarchy =  hierarchy)
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
            somatic_gene = somatic[somatic["Gene"] == gene] #select for all mutations for specified gene
            somatic_gene = somatic_gene.drop(columns = ["Gene"]) #drop the gene column due to every value being the same
            somatic_gene = somatic_gene.set_index("Clinical_Patient_Key") #set index as S** number for merging
            if not multiple_mutations:
                somatic_gene = self.add_mutation_hierarchy(somatic_gene) #appends hierachy for sorting so correct duplicate can be kept
                somatic_gene = somatic_gene.sort_values(by = ["Clinical_Patient_Key","Mutation_Hierarchy"], ascending = [True,False]) #sorts by patient key, then by hierarchy so the duplicates will come with the higher number first
                somatic_gene = somatic_gene[~somatic_gene.index.duplicated(keep="first")] #keeps first duplicate row if indices are the same
            merge = df_gene.join(somatic_gene, how = "left") #left join omics data and mutation data (left being the omics data)
            merge = merge.fillna(value = {'Mutation':"Wildtype"}) #fill in all Mutation NA values (no mutation data) as Wildtype
            merge["index"] = merge.index #set index values as column ??Do we need this? Just duplicates the index column.
            merge["Sample_Status"] = np.where(merge.index <= "S104", "Tumor", "Normal") #add patient type, setting all samples up to S104 as Tumor, others as normal.
            merge.loc[merge.Sample_Status == "Normal","Mutation"] = "Wildtype_Normal" #change all Wildtype for Normal samples to Wildtype_Normal
            merge.loc[merge.Mutation == "Wildtype","Mutation"] = "Wildtype_Tumor" #change all other Wildtype (should be for Tumor samples with imputed Wildtype value) to Wildtype_Tumor
            merge.name = df_gene.columns[0] + " omics data with " + gene + " mutation data"
            return merge
        else:
            print("Gene", gene, "not found in somatic mutations.")
    def get_phosphosites(self, phosphoproteomics, gene):
        """
        Parameters
        phosphoproteomics: the phosphoproteomics dataframe
        gene: the gene we want to get the phosphosites for

        Returns
        dataframe containing the phosphosites for the specified gene
        """
        regex = gene + ".*" #set regular expression using specified gene
        phosphosites = phosphoproteomics.filter(regex = (regex)) #find all columns that match the regular expression, aka, all phosphosites for the specified gene
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
        if gene in df1.columns and gene in df2.columns: #check provided gene is in both provided dataframes
            common = df1.index.intersection(df2.index) #get rows common to df1 and df2
            df1Matched = df1.loc[common] #select all common rows in df1
            df1Matched = df1Matched.sort_index() #sort rows in ascending order
            df2Matched = df2.loc[common] #select all common rows in df2
            df2Matched = df2Matched.sort_index() #sort rows in ascending order
            assert(hasattr(df1,"name")); assert(hasattr(df2,"name")) #check that both dataframes have a name, which is assigned at
            dict = {df1.name:df1Matched[gene], df2.name:df2Matched[gene]} #create prep dictionary for dataframe mapping name to specified gene column
            df = pd.DataFrame(dict, index = df1Matched.index) #create dataframe with common rows as rows, and dataframe name to specified gene column as columns
            df.name = gene #dataframe is named as specified gene
            return df
        else:
            if gene not in df1.columns:
                if gene not in df2.columns:
                    print(gene,"not found in either of the provided dataframes. Please check that the specified gene is included in both of the provided dataframes.")
                else:
                    print(gene, "not found in", df1.name, "dataframe. Please check that the specified gene is included in both of the provided dataframes.")
            else:
                if gene not in df2.columns:
                    print(gene, "not found in", df2.name, "dataframe. Please check that the specified gene is included in both of the provided dataframes.")
                else: #Shouldn't reach this branch
                    print("Error asserting",gene,"in",df1.name,"and",df2.name,"dataframes.")
    def compare_genes(self, df1, df2, genes):
        """
        Parameters
        df1: omics dataframe (proteomics) to be selected from
        df2: other omics dataframe (transcriptomics) to be selected from
        genes: gene or array of genes to select from each of the dataframes

        Returns
        Dataframe containing columns equal to the number of genes provided times two. Each two-column set is the data for each specified gene from the two specified dataframes
        """
        dfs = pd.DataFrame(index = df1.index.intersection(df2.index)) #create empty returnable dataframe with common rows of df1 and df2 as rows
        for gene in genes: #loop through list of genes provided
            df = Utilities().compare_gene(df1, df2, gene) #create temp dataframe per gene in list
            new_col1 = df1.name + "_" + gene #create first new column using first dataframe name and gene
            new_col2 = df2.name + "_" + gene #create second new column using second dataframe name and gene
            df = df.rename(columns = {df1.name:new_col1, df2.name:new_col2}) #rename columns in returned dataframe
            dfs = dfs.add(df, fill_value=0) #append temp dataframe onto returnable dataframe
        dfs.name = str(len(genes)) + " Genes Combined" #Name returnable dataframe using number of genes provided
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
            omics_gene_df = omics[[gene]] #select omics data for specified gene
            if duplicates: #don't filter out duplicate sample mutations
                return self.merge_somatic(somatic, gene, omics_gene_df, multiple_mutations = True)
            else: #filter out duplicate sample mutations
                merged_with_duplicates = self.merge_somatic(somatic, gene, omics_gene_df)
                merged = merged_with_duplicates[[gene, "Mutation", "Sample_Status"]]
                merged.name = merged_with_duplicates.name
                return merged
        elif omics.name.split("_")[0] == "phosphoproteomics":
            phosphosites = self.get_phosphosites(omics, gene)
            if len(phosphosites.columns) > 0:
                if duplicates:#don't filter out duplicate sample mutations
                    return self.merge_somatic(somatic, gene, phosphosites, multiple_mutations = True)
                else:#filter out duplicate sample mutations
                    columns = list(phosphosites.columns)
                    columns.append("Mutation")
                    columns.append("Sample_Status")
                    merged_somatic_full = self.merge_somatic(somatic, gene, phosphosites)
                    merged_somatic = merged_somatic_full[columns] #select all phosphosites, mutation, and patient type columns
                    merged_somatic.name = merged_somatic_full.name
                    return merged_somatic
        else:
            print("Gene", gene, "not found in", omics.name, "data")
    def merge_mutations_trans(self, omics, omicsGene, somatic, somaticGene, duplicates = False): #same function as merge_mutations, except use somaticGene to select mutation data
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
                merged_somatic_full = self.merge_somatic(somatic, somaticGene, omics_gene_df)
                merged_somatic = merged_somatic_full[[omicsGene, "Mutation", "Sample_Status"]]
                merged_somatic.name = merged_somatic_full.name
                return merged_somatic
        elif omics.name.split("_")[0] == "phosphoproteomics":
            phosphosites = self.get_phosphosites(omics, omicsGene)
            if len(phosphosites.columns) > 0:
                if duplicates:
                    return self.merge_somatic(somatic, somaticGene, phosphosites, multiple_mutations = True)
                else:
                    columns = list(phosphosites.columns)
                    columns.append("Mutation")
                    columns.append("Sample_Status")
                    merged_somatic_full = self.merge_somatic(somatic, somaticGene, phosphosites)
                    merged_somatic = merged_somatic_full[columns]
                    merged_somatic.name = merged_somatic_full.name
                    return merged_somatic
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
            df = data[data.columns] #new df variable prevents insert function overwriting the original data
            df.insert(0, clinical_col, clinical[clinical_col]) #insert specified clinical column into specified omics data
            df.name = data.name + " with " + clinical_col
            return df
        else:
            print(clinical_col, "not found in clinical dataframe. You can check the available columns using get_clinical_cols()")
    def compare_derived_molecular(self, derived_molecular, data, dm_col):
        """
        Parameters
        derived_molecular: derived molecular dataframe for omics data to be appended with
        data: omics data for derived molecular data to be appended with
        dm_col: column in derived molecular dataframe to be inserted into provided omics data

        Returns
        Dataframe with specified column from derived molecular dataframe added to specified dataframe (i.e., proteomics) for comparison and easy plotting
        """
        if dm_col in derived_molecular:
            df = data[data.columns]
            df.insert(0, dm_col, derived_molecular[dm_col])
            df.name = data.name + " with " + dm_col
            return df
        else:
            print(dm_col, "not found in derived_molecular dataframe. You can check the available columns using get_derived_molecular_cols()")
    def compare_phosphosites(self, proteomics, phosphoproteomics, gene):
        """
        Parameters
        gene: proteomics gene to query phosphoproteomics dataframe

        Searches for any phosphosites on the gene provided

        Returns
        Dataframe with a column from proteomics for the gene specified, as well as columns for all phosphoproteomics columns beginning with the specified gene
        """
        if gene in proteomics.columns:
            df = proteomics[[gene]] #select proteomics data for specified gene
            phosphosites = self.get_phosphosites(phosphoproteomics, gene) #gets phosphosites for specified gene
            if len(phosphosites.columns) > 0:
                df = df.add(phosphosites, fill_value=0) #adds phosphosites columns to proteomics data for specified gene
                df.name = gene + " proteomics and phosphoproteomics"
                return df
        else:
            print(gene, "not found in proteomics dataframe. Available genes can be checked using get_proteomics().columns")
