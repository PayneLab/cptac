import pandas as pd
import numpy as np

class Utilities:
    def __init__(self):
        pass
    def compare_gene(self, df1, df2, gene, key_id_map):
        """
        Parameters
        df1: omics dataframe (proteomics) to be selected from
        df2: other omics dataframe (transcriptomics) to be selected from
        gene: gene to select from each of the dataframes

        Returns
        Dataframe containing two columns. Each column is the data for the specified gene from the two specified dataframes
        """
        if (type(df1) != pd.DataFrame) or (type(df2) != pd.DataFrame):
            print("Provided data not a dataframe, please check that both data inputs are dataframes")
            return
        if gene in df1.columns and gene in df2.columns: #check provided gene is in both provided dataframes
            common = df1.set_index("patient_key").index.intersection(df2.set_index("patient_key").index) #select for intersection of patient keys between two dataframes
            df1Matched = df1.set_index("patient_key").loc[common] #select for rows matching common patient keys in df1
            df2Matched = df2.set_index("patient_key").loc[common] #select for rows matching common patient keys in df2
            assert(hasattr(df1,"name")); assert(hasattr(df2,"name")) #check that both dataframes have a name, which is assigned in DataFrameLoader class
            dict = {df1.name:df1Matched[gene], df2.name:df2Matched[gene]} #create prep dictionary for dataframe mapping name to specified gene column
            df = pd.DataFrame(dict, index = df1Matched.index) #create dataframe with common rows as rows, and dataframe name to specified gene column as columns
            df["patient_id"] = key_id_map[key_id_map["patient_key"].isin(list(df.index))].index
            df["patient_key"] = df.index
            df = df.set_index("patient_id")
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
    def compare_genes(self, df1, df2, genes, key_id_map):
        """
        Parameters
        df1: omics dataframe (proteomics) to be selected from
        df2: other omics dataframe (transcriptomics) to be selected from
        genes: gene or list of genes to select from each of the dataframes

        Returns
        Dataframe containing columns equal to the number of genes provided times two. Each two-column set is the data for each specified gene from the two specified dataframes
        """
        if (type(df1) != pd.DataFrame) or (type(df2) != pd.DataFrame):
            print("Provided data not a dataframe, please check that both data inputs are dataframes")
            return
        common = df1.set_index("patient_key").index.intersection(df2.set_index("patient_key").index)
        common_index = key_id_map[key_id_map["patient_key"].isin(list(common))].index
        dfs = pd.DataFrame(index = common_index) #create empty returnable dataframe with common rows of df1 and df2 as rows
        for gene in genes: #loop through list of genes provided
            df = Utilities().compare_gene(df1, df2, gene, key_id_map) #create temp dataframe per gene in list (can Utilities().compare_gene be changed to self.compare_gene?)
            new_col1 = df1.name + "_" + gene #create first new column using first dataframe name and gene
            new_col2 = df2.name + "_" + gene #create second new column using second dataframe name and gene
            df = df.rename(columns = {df1.name:new_col1, df2.name:new_col2}) #rename columns in returned dataframe
            dfs = pd.concat([dfs,df[df.columns[0:2]]], axis=1) #append temp dataframe onto returnable dataframe, leaving off patient_key column until the end
        dfs["patient_key"] = key_id_map.loc[dfs.index] #add patient_key column
        dfs.name = str(len(genes)) + " Genes Combined" #Name returnable dataframe using number of genes provided
        return dfs

    def compare_clinical(self, clinical, data, clinical_col, key_id_map):
        """
        Parameters
        clinical: clinical dataframe for omics data to be appended with
        data: omics data for clinical data to be appended with
        clinical_col: column in clinical dataframe to be inserted into provided omics data

        Returns
        Dataframe with specified column from clinical dataframe added to specified dataframe (i.e., proteomics) for comparison and easy plotting
        """
        if clinical_col in clinical:
            df = data[data.columns].set_index("patient_key") #prep returnable dataframe due to DataFrame.insert() changing by reference. If only df = data, then insert() will change data as well
            clinical = clinical.set_index("patient_key") #set index as patient key for mapping
            clinical = clinical.reindex(df.index) #select clinical data with indices matching omics data
            values = clinical[clinical_col] #get values for clinical column
            df.insert(0, clinical_col, values) #inserts specified clinical column at the beginning of the returnable dataframe
            df = df.assign(patient_id = key_id_map[key_id_map["patient_key"].isin(list(df.index))].index) #make patient id (non-S number)
            df = df.assign(patient_key = df.index) #move index to column (S number)
            df = df.set_index("patient_id") #set index as patient id (non-S number)
            df.name = data.name + " with " + clinical_col #assigns dataframe name using data name and specified clinical column
            return df
        else:
            print(clinical_col, "not found in clinical dataframe. You can check the available columns using get_clincal().columns")

    def get_phosphosites(self, phosphoproteomics, gene):
        """
        Parameters
        phosphoproteomics: the phosphoproteomics dataframe
        gene: the gene we want to get the phosphosites for

        Returns
        dataframe containing the phosphosites for the specified gene
        """
        regex = gene + "-.*" #set regular expression using specified gene
        phosphosites = phosphoproteomics.filter(regex = (regex)) #find all columns that match the regular expression, aka, all phosphosites for the specified gene
        if len(phosphosites.columns) == 0:
            print("Gene",gene, "not found in phosphoproteomics data")
        phosphosites.name = 'phosphosites_{}'.format(gene)
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
            df = proteomics[[gene]] #select proteomics data for specified gene
            phosphosites = self.get_phosphosites(phosphoproteomics, gene) #gets phosphosites for specified gene
            if len(phosphosites.columns) > 0:
                df = df.add(phosphosites, fill_value=0) #adds phosphosites columns to proteomics data for specified gene
                df.name = gene + " proteomics and phosphoproteomics"
                return df
        else:
            print(gene, "not found in proteomics dataframe. Available genes can be checked using get_proteomics().columns")
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
    def merge_somatic(self, somatic, gene, df_gene, key_id_map, multiple_mutations = False): #private
        """
        Parameters
        somatic: somatic mutations dataframe that will be parsed for specified gene data
        gene: string of gene to be selected for in somatic mutation data
        df_gene: selection of omics data for particular gene to be merged with somatic data
        multiple_mutations: boolean indicating whether to include multiple mutations for specified gene in an individual

        Returns
        Dataframe of merged somatic and omics dataframes based on gene provided
        """
        #TODO: use patient_key instead of patient_id, proteomics currently returning all na, therefore 155 gives all wildtypeov
        if sum(somatic["Gene"] == gene) > 0:
            somatic_gene = somatic[somatic["Gene"] == gene] #select for all mutations for specified gene
            somatic_gene = somatic_gene.drop(columns = ["Gene"]) #drop the gene column due to every value being the same
            somatic_gene = somatic_gene.set_index("patient_key") #set index as patient key (S number)
            if not multiple_mutations: #if you want to remove duplicate indices
                somatic_gene = self.add_mutation_hierarchy(somatic_gene) #appends hierachy for sorting so correct duplicate can be kept
                somatic_gene["forSort"] = somatic_gene.index.str[1:] #creates separate column for sorting
                somatic_gene[["forSort"]] = somatic_gene[["forSort"]].apply(pd.to_numeric) #converts string column of patient key numbers to floats for sorting
                somatic_gene = somatic_gene.sort_values(by = ["forSort","Mutation_Hierarchy"], ascending = [True,False]) #sorts by patient key, then by hierarchy so the duplicates will come with the lower number first
                somatic_gene = somatic_gene.drop(columns=["forSort"]) #drops sorting column
                somatic_gene = somatic_gene[~somatic_gene.index.duplicated(keep="first")] #keeps first duplicate row if indices are the same
            merge = df_gene.join(somatic_gene, how = "left") #merges dataframes based on indices, how = "left" defaulting to the df_gene indices. If indices don't match, then mutation column will be NaN
            merge[["Mutation"]] = merge[["Mutation"]].fillna(value = "Wildtype") #fill in all Mutation NA values (no mutation data) as Wildtype
            if multiple_mutations:
                patient_ids = []
                patient_keys = list(merge.index)
                for key in patient_keys:
                    patient_ids.append(key_id_map[key_id_map["patient_key"] == key].index.values[0])
                assert(len(patient_ids) == len(patient_keys))
                merge["patient_id"] = patient_ids
                merge = merge.sort_values(by = ["patient_id"])
            else:
                merge["patient_id"] = key_id_map[key_id_map["patient_key"].isin(list(merge.index))].index #reverse lookup for patient key (S number) to patient id (non-S number)
            merge["patient_key"] = merge.index #move index to column
            merge = merge.set_index("patient_id") #set index to patient id (non-S number)
            merge["Sample_Status"] = np.where(merge.index <= "26OV013", "Tumor", "Normal") #26OV013 is the last patient id before the "N******" ids
            merge.loc[merge.Sample_Status == "Normal","Mutation"] = "Wildtype_Normal" #change all Wildtype for Normal samples to Wildtype_Normal
            merge.loc[merge.Mutation == "Wildtype","Mutation"] = "Wildtype_Tumor" #change all other Wildtype (should be for Tumor samples with imputed Wildtype value) to Wildtype_Tumor
            merge = merge.fillna(value={'Location':'No_mutation'}) # If there's no location, there wasn't a mutation--make it easier for people to understand what that means.
            merge.name = df_gene.columns[0] + " omics data with " + gene + " mutation data"
            return merge
        else:
            print("Gene", gene, "not found in somatic mutations.")
    def merge_mutations(self, omics, somatic, gene, key_id_map, duplicates = False):
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
            omics_gene_df = omics[[gene,"patient_key"]].set_index("patient_key") #select omics data for specified gene, setting index to patient key (S number) for merging
            if duplicates: #TODO: this breaks right now, merge_somatic can't handle duplicate samples
                return self.merge_somatic(somatic, gene, omics_gene_df, key_id_map, multiple_mutations = True)
            else: #filter out duplicate sample mutations
                return self.merge_somatic(somatic, gene, omics_gene_df, key_id_map)[[gene, "Mutation", "patient_key", "Sample_Status"]]
        elif omics.name.split("_")[0] == "phosphoproteomics":
            phosphosites = self.get_phosphosites(omics, gene)
            if len(phosphosites.columns) > 0:
                phosphosites = phosphosites.assign(patient_key = omics["patient_key"])
                phosphosites = phosphosites.set_index("patient_key")
                if duplicates:
                    return self.merge_somatic(somatic, gene, phosphosites, key_id_map, multiple_mutations = True)
                else:
                    columns = list(phosphosites.columns)
                    columns.append("Mutation")
                    columns.append("patient_key")
                    columns.append("Sample_Status")
                    merged_somatic = self.merge_somatic(somatic, gene, phosphosites, key_id_map)
                    return merged_somatic[columns]

        else:
            print("Gene", gene, "not found in", omics.name, "data")
    def merge_mutations_trans(self, omics, omicsGene, somatic, somaticGene, key_id_map, duplicates = False):  #same functonality as merge_mutations, but uses somaticGene for merge_somatic
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
        merged_somatic = None
        if omicsGene in omics.columns:
            omics_gene_df = omics[[omicsGene,"patient_key"]].set_index("patient_key")
            if duplicates:
                merged_somatic = self.merge_somatic(somatic, somaticGene, omics_gene_df, key_id_map, multiple_mutations = True)
            else:
                merged_somatic = self.merge_somatic(somatic, somaticGene, omics_gene_df, key_id_map)[[omicsGene, "Mutation", "patient_key", "Sample_Status"]]
        elif omics.name.split("_")[0] == "phosphoproteomics":
            phosphosites = self.get_phosphosites(omics, omicsGene)
            if len(phosphosites.columns) > 0:
                phosphosites = phosphosites.assign(patient_key = omics["patient_key"])
                phosphosites = phosphosites.set_index("patient_key")
                if duplicates:
                    merged_somatic = self.merge_somatic(somatic, somaticGene, phosphosites, key_id_map, multiple_mutations = True)
                else:
                    columns = list(phosphosites.columns)
                    columns.append("Mutation")
                    columns.append("patient_key")
                    columns.append("Sample_Status")
                    merged_somatic = self.merge_somatic(somatic, somaticGene, phosphosites, key_id_map)
                    merged_somatic =  merged_somatic[columns]
        else:
            print("Gene", omicsGene, "not found in", omics.name,"data")
            return
        if merged_somatic is None:
            return
        merged_somatic = merged_somatic.rename(columns={omicsGene:omicsGene + '_omics', 'Mutation':somaticGene + '_Mutation', 'Location':somaticGene + '_Location'}) # Add the gene name to the column headers, so that it's clear which gene the data is for.
        return merged_somatic
