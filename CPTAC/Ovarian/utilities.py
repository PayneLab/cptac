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
    def add_mutation_hierarchy(self, somatic): #private
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
            somatic_gene = somatic_gene.set_index("Patient_Id")
            if not multiple_mutations:
                somatic_gene = self.add_mutation_hierarchy(somatic_gene) #appends hierachy for sorting so correct duplicate can be kept
                somatic_gene = somatic_gene.sort_values(by = ["Patient_Id","Mutation_Hierarchy"], ascending = [True,False]) #sorts by patient key, then by hierarchy so the duplicates will come with the lower number first
                somatic_gene = somatic_gene[~somatic_gene.index.duplicated(keep="first")] #keeps first duplicate row if indices are the same
            merge = df_gene.join(somatic_gene, how = "left")
            merge = merge.fillna(value = {'Mutation':"Wildtype"})
            #merge["index"] = merge.index
            #merge["Patient_Type"] = np.where(merge.index <= "S100", "Tumor", "Normal")
            merge.name = df_gene.columns[0] + " omics data with " + gene + " mutation data"
            return merge
        else:
            print("Gene", gene, "not found in somatic mutations.")
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
                return self.merge_somatic(somatic, gene, omics_gene_df)[[gene, "Mutation"]]#, "Patient_Type"]]
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
