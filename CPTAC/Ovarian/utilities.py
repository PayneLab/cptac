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
        if gene in df1.columns and gene in df2.columns: #check provided gene is in both provided dataframes
            common = df1.index.intersection(df2.index) #get rows common to df1 and df2
            df1Matched = df1.loc[common] #select all common rows in df1
            df1Matched = df1Matched.sort_index() #sort rows in ascending order
            df2Matched = df2.loc[common] #select all common rows in df2
            df2Matched = df2Matched.sort_index() #sort rows in ascending order
            assert(hasattr(df1,"name")); assert(hasattr(df2,"name")) #check that both dataframes have a name, which is assined at
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
        genes: gene or list of genes to select from each of the dataframes

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
            df = data[data.columns] #prep returnable dataframe due to DataFrame.insert() changing by reference. If only df = data, then insert() will change data as well
            if len(clinical.index) != len(df.index):
                clinical = clinical.reindex(df.index)
            values = clinical[clinical_col]
            df.insert(0, clinical_col, values) #inserts specified clinical column at the beginning of the returnable dataframe
            df.name = data.name + " with " + clinical_col #assigns dataframe name using data name and specified clinical column
            return df
        else:
            print(clinical_col, "not found in clinical dataframe. You can check the available columns by entering CPTAC.get_clincal().columns")

    def get_phosphosites(self, phosphoproteomics, gene):
        regex = gene + ".*"
        phosphosites = phosphoproteomics.filter(regex = (regex))
        if len(phosphosites.columns) == 0:
            print("Gene",gene, "not found in phosphoproteomics data")
        else:
            return phosphosites

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
