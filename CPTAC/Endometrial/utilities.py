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
            #merge["index"] = merge.index #set index values as column
            merge["Sample_Status"] = np.where(merge.index <= "S104", "Tumor", "Normal") #add patient type, setting all samples up to S104 as Tumor, others as normal.
            merge.loc[merge.Sample_Status == "Normal","Mutation"] = "Wildtype_Normal" #change all Wildtype for Normal samples to Wildtype_Normal
            merge.loc[merge.Mutation == "Wildtype","Mutation"] = "Wildtype_Tumor" #change all other Wildtype (should be for Tumor samples with imputed Wildtype value) to Wildtype_Tumor
            merge = merge.drop(columns=['Patient_Id']) # We don't need this column
            merge = merge.fillna(value={'Location':'No_mutation'}) # If there's no location, there wasn't a mutation--make it easier for people to understand what that means.
            merge.name = df_gene.columns[0] + " omics data with " + gene + " mutation data"
            return merge
        else:
            print("Gene", gene, "not found in somatic mutations.")
# Obsolete. Replaced by get_col_from_omics.
#    def get_phosphosites(self, phosphoproteomics, gene):
#        """
#        Parameters
#        phosphoproteomics: the phosphoproteomics dataframe
#        gene: the gene we want to get the phosphosites for
#
#        Returns
#        dataframe containing the phosphosites for the specified gene
#        """
#        regex = gene + "-.*" #set regular expression using specified gene
#        phosphosites = phosphoproteomics.filter(regex = (regex)) #find all columns that match the regular expression, aka, all phosphosites for the specified gene
#        if len(phosphosites.columns) == 0:
#            print("Gene",gene, "not found in phosphoproteomics data")
#        phosphosites.name = 'phosphosites_{}'.format(gene)
#        return phosphosites

    def get_col_from_omics(omics_df, key): # private
        """Based on a single key, select a column or columns from an omics dataframe. Handles key parsing if needed, based on dataframe name.

        Parameters:
        omics_df (pandas.core.frame.DataFrame): omics dataframe to select colum(s) from.
        key (str): key to use to select column(s). 

        Returns: 
        pandas.core.frame.DataFrame: The selected column(s) from the dataframe.
        """
        if omics_df.name == ('phosphoproteomics_site' or 'acetylproteomics'): # These are the dataframes that require us to parse the key
            regex = key + "-.*" # Build a regex to get all columns that match the key
        else:
            regex = key

        selected = omics_df.filter(regex = (regex)) # Find all columns that match the key. If only one column matches, DataFrame.filter will still return a dataframe, not a series :)
        if len(selected.columns) == 0: # If none of the columns matched the key, print an error message and return None.
            print('{} did not match any columns in {} dataframe. Please double check that it is included in the dataframe.'.format(key, omics_df.name))
            return
        selected.name = "{} for {}".format(omics_df.name, key) # Give it a name!
        return selected

    def get_cols_from_omics(omics_df, keys): # private
        """Based on a list of keys, select multiple columns from an omics dataframe, and return the selected columns as one dataframe.

        Parameters:
        omics_df (pandas.core.frame.DataFrame): omics dataframe to select column(s) from.
        keys (list): list of keys, as strings, to use to select columns.

        Returns:
        pandas.core.frame.DataFrame: The selected columns from the dataframe.
        """
        df = pd.DataFrame(index=omics_df.index) # Create an empty dataframe, which we'll fill with the columns we select using our keys, and then return.
        for key in keys:
            selected = self.get_col_from_omics(omics_df, key) # Get the columns for that key from the dataframe
            if selected is None: # If it didn't match any columns, get_col_from_omics will have printed an error message. Return None.
                return
            df = df.add(selected, fill_value=0) # Otherwise, append the columns to our dataframe we'll return.
        df.name = "{} for {} genes".format(omics_df.name, len(keys)) # NAME the dataFRAME!
        return df

    def select_omics_from_str_or_list(omics_df, keys): # private
        """Determines whether you passed it a single key or a list of keys, selects the corresponding columns from the omics dataframe, and returns them.

        Parameters:
        omics_df (pandas.core.frame.DataFrame): omics dataframe to select columns from
        keys (str or list): key(s) to use to select columns from omics_df. str if one key, list if multiple.

        Returns:
        pandas.core.frame.DataFrame: the columns from omics_df corresponding to the key(s), as a dataframe.
        """
        if isinstance(keys, str): # If it's a single key, feed it to the proper function
            return self.get_col_from_omics(omics_df, keys) 
        elif isinstance(keys, list): # If it's a list of keys, feed it to the proper function
            return self.get_cols_from_omics(omics_df, keys)
        else: # If it's neither, they done messed up. Tell 'em.
            print("Keys parameter {} is of invalid type {}. Valid types: str or list.".format(keys, type(keys)))

    def compare_omics(self, df1, keys1, df2, keys2):
        """Select columns for one key or a list of keys from one omics dataframe, and columns for another key or list of keys from another omics dataframe, and join them into one dataframe.

        Parameters:
        df1 (pandas.core.frame.DataFrame): first omics dataframe to select from.
        keys1 (str or list): key(s) for column(s) to select from the first omics dataframe. str if one key, list of strings if multiple keys.
        df2 (pandas.core.frame.DataFrame): second omics dataframe to select from.
        keys2 (str or list): key(s) for column(s) to select from the second omics dataframe. str if one key, list of strings if multiple keys.

        Returns:
        pandas.core.frame.DataFrame: The data from the selected columns from the two dataframes, joined into one dataframe.
        """
        selected1 = select_omics_from_str_or_list(df1, keys1)
        selected2 = select_omics_from_str_or_list(df2, keys2)

        if (selected1 is not None) and (selected2 is not None): # If either selector returned None, the key(s) didn't match any columns, and it printed an informative error message already. We'll return None.
            df = selected1.join(selected2, how='inner') # Join the rows common to both dataframes
            df = df.sort_index() # Sort rows in ascending order
            df.name = "{} with {}".format(selected1.name, selected2.name) # Give it a nice name identifying the data in it.
            return df

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
            phosphosites = self.get_col_from_omics(omics, gene)
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
        merged_somatic = None
        if omicsGene in omics.columns:
            omics_gene_df = omics[[omicsGene]]
            if duplicates:
                merged_somatic = self.merge_somatic(somatic, somaticGene, omics_gene_df, multiple_mutations = True)
            else:
                merged_somatic_full = self.merge_somatic(somatic, somaticGene, omics_gene_df)
                merged_somatic = merged_somatic_full[[omicsGene, "Mutation", "Sample_Status"]]
                merged_somatic.name = merged_somatic_full.name
        elif omics.name.split("_")[0] == "phosphoproteomics":
            phosphosites = self.get_col_from_omics(omics, omicsGene)
            if len(phosphosites.columns) > 0:
                if duplicates:
                    merged_somatic = self.merge_somatic(somatic, somaticGene, phosphosites, multiple_mutations = True)
                else:
                    columns = list(phosphosites.columns)
                    columns.append("Mutation")
                    columns.append("Sample_Status")
                    merged_somatic_full = self.merge_somatic(somatic, somaticGene, phosphosites)
                    merged_somatic = merged_somatic_full[columns]
                    merged_somatic.name = merged_somatic_full.name
        else:
            print("Gene", omicsGene, "not found in", omics.name,"data")
            return
        if merged_somatic is None:
            return
        merged_somatic = merged_somatic.rename(columns={omicsGene:omicsGene + '_omics', 'Mutation':somaticGene + '_Mutation', 'Location':somaticGene + '_Location'}) # Add the gene name to the column headers, so that it's clear which gene the data is for.
        return merged_somatic

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
# Obsolete. Replaced by compare_omics.
#    def compare_phosphosites(self, proteomics, phosphoproteomics, gene):
#        """
#        Parameters
#        gene: proteomics gene to query phosphoproteomics dataframe
#
#        Searches for any phosphosites on the gene provided
#
#        Returns
#        Dataframe with a column from proteomics for the gene specified, as well as columns for all phosphoproteomics columns beginning with the specified gene
#        """
#        if gene in proteomics.columns:
#            df = proteomics[[gene]] #select proteomics data for specified gene
#            phosphosites = self.get_col_from_omics(phosphoproteomics, gene) #gets phosphosites for specified gene
#            if len(phosphosites.columns) > 0:
#                df = df.add(phosphosites, fill_value=0) #adds phosphosites columns to proteomics data for specified gene
#                df.name = gene + " proteomics and phosphoproteomics"
#                return df
#        else:
#            print(gene, "not found in proteomics dataframe. Available genes can be checked using get_proteomics().columns")
