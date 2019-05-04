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
import re

class Utilities:

    def __init__(self):
        pass
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

# Next 4 functions are for working with omics data
    def get_col_from_omics(self, omics_df, gene): # private
        """Based on a single gene, select a column or columns from an omics dataframe. If dataframe is phospho- or acetylproteomics, grabs all columns that match the gene.

        Parameters:
        omics_df (pandas.core.frame.DataFrame): omics dataframe to select colum(s) from.
        gene (str): gene to use to select column(s). 

        Returns: 
        pandas.core.frame.DataFrame: The selected column(s) from the dataframe.
        """
        if omics_df.name == ('phosphoproteomics_site' or 'acetylproteomics'):
            col_regex = gene + "-.*" # Build a regex to get all columns that match the gene
        else:
            col_regex = '^{}$'.format(gene)

        selected = omics_df.filter(regex=col_regex) # Find all columns that match the gene. If only one column matches, DataFrame.filter will still return a dataframe, not a series :)
        if len(selected.columns) == 0: # If none of the columns matched the gene, print an error message and return None.
            print('{} did not match any columns in {} dataframe. Please double check that it is included in the dataframe.'.format(gene, omics_df.name))
            return
        selected = selected.rename(columns=lambda x:'{}_{}'.format(x, omics_df.name)) # Append dataframe name to end of each column header, to preserve info when we merge dataframes
        selected.name = "{} for {}".format(omics_df.name, gene) # Give the dataframe a name!
        return selected

    def get_cols_from_omics(self, omics_df, genes): # private
        """Based on a list of genes, select multiple columns from an omics dataframe, and return the selected columns as one dataframe.

        Parameters:
        omics_df (pandas.core.frame.DataFrame): omics dataframe to select column(s) from.
        genes (list): list of genes, as strings, to use to select columns.

        Returns:
        pandas.core.frame.DataFrame: The selected columns from the dataframe.
        """
        df = pd.DataFrame(index=omics_df.index) # Create an empty dataframe, which we'll fill with the columns we select using our genes, and then return.
        for gene in genes:
            selected = self.get_col_from_omics(omics_df, gene) # Get the columns for that gene from the dataframe
            if selected is None: # If it didn't match any columns, get_col_from_omics will have printed an error message. Return None.
                return
            df = df.join(selected, how='left') # Otherwise, append the columns to our dataframe we'll return.
        df.name = "{} for {} genes".format(omics_df.name, len(genes)) # Name the dataframe!
        return df

    def get_omics_from_str_or_list(self, omics_df, genes):
        """Determines whether you passed it a single gene or a list of genes, selects the corresponding columns from the omics dataframe, and returns them.

        Parameters:
        omics_df (pandas.core.frame.DataFrame): omics dataframe to select columns from
        genes (str or list): gene(s) to use to select columns from omics_df. str if one gene, list if multiple. Passing None will select the entire omics dataframe.

        Returns:
        pandas.core.frame.DataFrame: the columns from omics_df corresponding to the gene(s), as a dataframe.
        """
        if isinstance(genes, str): # If it's a single gene, feed it to the proper function
            return self.get_col_from_omics(omics_df, genes) 
        elif isinstance(genes, list): # If it's a list of genes, feed it to the proper function
            return self.get_cols_from_omics(omics_df, genes)
        elif genes is None: # If it's the default of None, rename columns and return the entire dataframe
            return_df = omics_df.rename(columns=lambda x:'{}_{}'.format(x, omics_df.name)) # Append dataframe name to end of each column header, to preserve info when we merge dataframes
            return_df.name = omics_df.name # Name the return dataframe
            return return_df
        else: # If it's none of those, they done messed up. Tell 'em.
            print("Genes parameter {} is of invalid type {}. Valid types: str, list, or NoneType.".format(genes, type(genes)))

    def compare_omics(self, df1, df2, genes1, genes2):
        """Select columns for one gene or a list of genes from one omics dataframe, and columns for another gene or list of genes from another omics dataframe, and join them into one dataframe.

        Parameters:
        df1 (pandas.core.frame.DataFrame): first omics dataframe to select from.
        genes1 (str or list): gene(s) for column(s) to select from the first omics dataframe. str if one gene, list of strings if multiple genes. Passing None will select the entire omics dataframe.
        df2 (pandas.core.frame.DataFrame): second omics dataframe to select from.
        genes2 (str or list): gene(s) for column(s) to select from the second omics dataframe. str if one gene, list of strings if multiple genes. Passing None will select the entire omics dataframe.

        Returns:
        pandas.core.frame.DataFrame: The data from the selected columns from the two dataframes, joined into one dataframe.
        """
        selected1 = self.get_omics_from_str_or_list(df1, genes1)
        selected2 = self.get_omics_from_str_or_list(df2, genes2)

        if (selected1 is not None) and (selected2 is not None): # If either selector returned None, the gene(s) didn't match any columns, and it printed an informative error message already. We'll return None.
            df = selected1.join(selected2, how='inner') # Join the rows common to both dataframes
            df = df.sort_index() # Sort rows in ascending order
            df.name = "{}, with {}".format(selected1.name, selected2.name) # Give it a nice name identifying the data in it.
            return df

# Next 4 functions are for working with clinical and derived molecular data
    def get_col_from_clinical_or_derived_molecular(self, df, col):
        """Get a single column from the clinical or derived_molecular dataframe.

        Parameters:
        df (pandas.core.frame.DataFrame): The dataframe to select the column from. Either clinical or derived_molecular.
        col (str): The column to select from the dataframe.

        Returns:
        pandas.core.frame.DataFrame: The specified column from the given dataframe. We return it as a dataframe for easy merging later.
        """
        if col not in df.columns.values: # If they didn't give us one of the actual columns, tell them and return None.
            print('{} column not found in {} dataframe. Please double check that it is included in the dataframe.'.format(col, df.name))
            return
        
        selected = df.loc[:, [col]] # Select the column from the dataframe, keeping it as a dataframe
        selected.name = '{} from {}'.format(col, df.name) # Give it a name, identifying which dataframe it came from
        return selected

    def get_cols_from_clinical_or_derived_molecular(self, df, cols):
        """Select several columns, given as a list, from either the clinical or derived_molecular dataframe.

        Parameters:
        df (pandas.core.frame.DataFrame): The dataframe to select the column from. Either clinical derived_molecular.
        cols (list): A list of the columns, as strings, to select from the dataframe.

        Returns:
        pandas.core.frame.DataFrame: The specified columns from the given dataframe.
        """
        return_df = pd.DataFrame(index=df.index) # Create an empty dataframe, which we'll fill with the columns we select, and then return.
        for col in cols:
            selected = self.get_col_from_clinical_or_derived_molecular(df, col) # Get the columns from the given dataframe
            if selected is None: # If it didn't match any columns, get_col_from_omics will have printed an error message. Return None.
                return
            return_df = return_df.join(selected, how='left') # Otherwise, append the columns to our dataframe we'll return.
        return_df.name = "{} columns from {}".format(len(cols), df.name) # Name the dataframe!
        return return_df

    def get_clinical_or_derived_molecular_from_str_or_list(self, df, cols):
        """Determines whether you passed it a single column or a list of columns, then selects them from the given dataframe, and returns them.

        Parameters:
        df (pandas.core.frame.DataFrame): The dataframe to select column(s) from. Either clinical or derived_molecular.
        cols (str or list): Column(s) to select from the dataframe. str if one, list if multiple.

        Returns:
        pandas.core.frame.DataFrame: The specificed columns from the given dataframe.
        """
        if isinstance(cols, str): # If it's a single column, feed it to the proper function
            return self.get_col_from_clinical_or_derived_molecular(df, cols)
        elif isinstance(cols, list): # If it's a list of columns, feed it to the proper function
            return self.get_cols_from_clinical_or_derived_molecular(df, cols)
        else: # If it's neither of those, they done messed up. Tell 'em.
            print("Columns parameter {} is of invalid type {}. Valid types: str or list.".format(cols, type(cols)))

    def append_clinical_or_derived_molecular_to_omics(self, df, omics_df, df_cols, omics_cols):
        """Selected the specified columns from either the clinical or derived_molecular dataframe, and append to all or part of the given omics dataframe.

        Parameters:
        df (pandas.core.frame.DataFrame): Either the clinical or derived_molecular dataframe, from which we'll select our columns to append.
        omics_df (pandas.core.frame.DataFrame): The omics dataframe to append to.
        df_cols (str or list): Column(s) to select from the clinical or derived_molecular dataframe. str if one column, list of strings if multiple.
        omics_cols (str or list): Gene(s) to select data for from the omics dataframe. Passing None will select the entire omics dataframe.

        Returns:
        pandas.core.frame.DataFrame: The selected columns from the clinical or derived_molecular dataframe, appended to the selected columns from the omics dataframe.
        """  
        df_selected = self.get_clinical_or_derived_molecular_from_str_or_list(df, df_cols)
        omics_selected = self.get_omics_from_str_or_list(omics_df, omics_cols)

        if (df_selected is not None) and (omics_selected is not None): # If either selector returned None, the key(s) didn't match any columns, and it printed an informative error message already. We'll return None.
            df = df_selected.join(omics_selected, how='inner') # Join the rows common to both dataframes
            df = df.sort_index() # Sort rows in ascending order
            df.name = "{}, with {}".format(df_selected.name, omics_selected.name) # Give it a nice name identifying the data in it.
            return df

# Next 4 functions are for working with mutation data
    def get_mutations_for_gene(self, somatic, gene, multiple_mutations):
        """Gets all the mutations for a specific gene, for all patients.

        Parameters:
        somatic (pandas.core.frame.DataFrame): The somatic mutation dataframe that we'll grab the mutation data from.
        gene (str): The gene to grab mutations for.
        multiple_mutations (bool): Whether to keep multiple mutations on the gene for each patient, or only report the highest priority mutation per patient. 

        Returns:
        pandas.core.frame.DataFrame: The mutations in each patient for the specified gene.
        """
        mutations = somatic[somatic["Gene"] == gene]

        if len(mutations) == 0: # If the gene doesn't match any in the dataframe, tell them, and return None.
            print("{} gene not found in somatic mutations data.".format(gene))
            return
        
        mutations = mutations.drop(columns=["Gene"]) # Drop the gene column due to every value being the same
        mutations = mutations.set_index("Clinical_Patient_Key") # Set index as S*** number for merging
        mutations = mutations.drop(columns=['Patient_Id']) # We don't need this column
        if not multiple_mutations: # Filter out multiple mutations for a single sample
            mutations = self.add_mutation_hierarchy(mutations) # Appends hierachy for sorting so correct duplicate can be kept
            mutations = mutations.sort_values(by = ["Clinical_Patient_Key","Mutation_Hierarchy"], ascending = [True,False]) # Sorts by patient key, then by hierarchy so the duplicates will come with the higher number first
            mutations = mutations[~mutations.index.duplicated(keep="first")] # Keeps first duplicate row if indices are the same
            mutations = mutations.drop(columns=['Mutation_Hierarchy']) # Drop the Mutation_Hierarchy now that we've dropped the duplicates

        mutations = mutations.rename(columns=lambda x:'{}_{}'.format(gene, x)) # Add the gene name to end beginning of each column header, to preserve info when we merge dataframes
        mutations.name = 'Somatic mutation data for {} gene'.format(gene)
        return mutations

    def get_mutations_for_genes(self, somatic, genes, multiple_mutations):
        """Gets all the mutations for a list of genes, for all patients.

        Parameters:
        somatic (pandas.core.frame.DataFrame): The somatic mutation dataframe that we'll grab the mutation data from.
        genes (list): The genes, as strings, to grab mutations for.
        multiple_mutations (bool): Whether to keep multiple mutations on a single gene for each patient, or only report the highest priority mutation per patient. 

        Returns:
        pandas.core.frame.DataFrame: The mutations in each patient for the specified genes.
        """
        df = pd.DataFrame(index=somatic['Clinical_Patient_Key'].drop_duplicates()) # Create an empty dataframe, which we'll fill with the columns we select using our genes, and then return.
        for gene in genes:
            selected = self.get_mutations_for_gene(somatic, gene, multiple_mutations) # Get the mutations for our gene
            if selected is None: # If there's no mutation data for that gene, get_mutations_for_gene will have printed an error message. Return None.
                return
            df = df.join(selected, how='left') # Otherwise, append the columns to our dataframe we'll return.
        df.name = "Somatic mutation data for {} genes".format(len(genes)) # Name the dataframe!
        return df

    def get_mutations_from_str_or_list(self, somatic, genes, multiple_mutations):
        """Determines whether you passed it a single gene or a list of genes, selects the corresponding mutations from the somatic mutation dataframe, and returns them.

        Parameters:
        somatic (pandas.core.frame.DataFrame): The somatic mutation dataframe we'll grab the mutation data from.
        genes (str or list): gene(s) to select mutations for. str if one gene, list if multiple.
        multiple_mutations (bool): Whether to keep multiple mutations on the gene for each patient, or only report the highest priority mutation per patient. 

        Returns:
        pandas.core.frame.DataFrame: the mutation data corresponding to the gene(s), as a dataframe.
        """
        if isinstance(genes, str): # If it's a single gene, feed it to the proper function
            return self.get_mutations_for_gene(somatic, genes, multiple_mutations) 
        elif isinstance(genes, list): # If it's a list of genes, feed it to the proper function
            return self.get_mutations_for_genes(somatic, genes, multiple_mutations)
        else: # If it's neither of those, they done messed up. Tell 'em.
            print("Genes parameter {} is of invalid type {}. Valid types: str or list.".format(genes, type(genes)))

    def append_mutations_to_omics(self, somatic, omics_df, mutation_genes, omics_genes, multiple_mutations, show_location):
        """Select all mutations for specified gene(s), and append to all or part of the given omics dataframe.

        Parameters:
        somatic (pandas.core.frame.DataFrame): Somatic mutation dataframe we'll get the dataframe.
        omics_df (pandas.core.frame.DataFrame): Omics dataframe to append the mutation data to.
        mutation_genes (str or list): The gene(s) to get mutation data for. str if one gene, list if multiple.
        omics_genes (str or list): Gene(s) to select from the omics dataframe. str if one gene, list if multiple. Passing None will select the entire omics dataframe.
        multiple_mutations (bool): Whether to keep multiple mutations on the gene for each patient, or only report the highest priority mutation per patient. 
        show_location (bool): Whether to include the Location column from the mutation dataframe.

        Returns:
        pandas.core.frame.DataFrame: The mutations for the specified gene, appended to all or part of the omics dataframe.
        """
        omics = self.get_omics_from_str_or_list(omics_df, omics_genes)
        mutations = self.get_mutations_from_str_or_list(somatic, mutation_genes, multiple_mutations)

        if (omics is not None) and (mutations is not None): # If either selector returned None, then there were gene(s) that didn't match anything, and an error message was printed. We'll return None.
            merge = omics.join(mutations, how = "left") # Left join omics data and mutation data (left being the omics data)

            merge["Sample_Status"] = np.where(merge.index <= "S104", "Tumor", "Normal") # Add Sample_Status column indicating sample type. Set all samples up to S104 as Tumor, others as normal.

            mutation_regex = r'.*_Mutation' # Construct regex to find all mutation columns
            mutation_cols = [col for col in merge.columns.values if re.match(mutation_regex, col)] # Get a list of all mutation columns
            for mutation_col in mutation_cols:
                merge.loc[(merge['Sample_Status'] == "Normal") & (pd.isnull(merge[mutation_col])), mutation_col] = "Wildtype_Normal" # Change all NaN mutation values (i.e., no mutation data for that sample) for Normal samples to Wildtype_Normal
                merge.loc[(merge['Sample_Status'] == "Tumor") & (pd.isnull(merge[mutation_col])), mutation_col] = "Wildtype_Tumor" # Change all NaN mutation values (i.e., no mutation data for that sample) for Tumor samples to Wildtype_Tumor

            location_regex = r'.*_Location' # Construct regex to find all location columns
            location_cols = [col for col in merge.columns.values if re.match(location_regex, col)] # Get a list of all location columns
            for location_col in location_cols:
                if show_location:
                    merge = merge.fillna(value={location_col:'No_mutation'}) # If there's no location, there wasn't a mutation--make it easier for people to understand what that means.
                else:
                    merge = merge.drop(columns=[location_col]) # Drop the location column, if the caller wanted us to.

            merge.name = "{}, with {}".format(omics.name, mutations.name) # Give it a name identifying the data in it
            return merge
