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

# Next 4 functions are for working with omics data
    def get_col_from_omics(self, omics_df, gene): 
        """Based on a single gene, select a column or columns from an omics dataframe. If dataframe is phosphoproteomics or acetylproteomics, grabs all columns that match the gene.

        Parameters:
        omics_df (pandas DataFrame): omics dataframe to select colum(s) from.
        gene (str): gene to use to select column(s). 

        Returns: 
        pandas DataFrame: The selected column(s) from the dataframe.
        """
        if omics_df.name == 'phosphoproteomics' or omics_df.name == 'acetylproteomics':
            col_regex = gene + "-.*" # Build a regex to get all columns that match the gene
        else:
            col_regex = '^{}$'.format(gene)

        selected = omics_df.filter(regex=col_regex) # Find all columns that match the gene. If only one column matches, DataFrame.filter will still return a dataframe, not a series :)
        if len(selected.columns) == 0: # If none of the columns matched the gene, print a warning message and return a column of NaN
            print('{0} did not match any columns in {1} dataframe. {0}_{1} column inserted, but filled with NaN.'.format(gene, omics_df.name))
            omics_index_copy = omics_df.index.copy()
            empty_omics_df = pd.DataFrame(index=omics_index_copy)
            selected = empty_omics_df.assign(**{gene:np.nan}) # ** unpacks the dictionary, creating a column with the value of the gene variable as the header, and filled with NaNs. Weird pandas syntax.
        selected = selected.rename(columns=lambda x:'{}_{}'.format(x, omics_df.name)) # Append dataframe name to end of each column header, to preserve info when we merge dataframes
        selected.name = "{} for {}".format(omics_df.name, gene) # Give the dataframe a name!
        return selected

    def get_cols_from_omics(self, omics_df, genes):
        """Based on a list or array-like of genes, select multiple columns from an omics dataframe, and return the selected columns as one dataframe.

        Parameters:
        omics_df (pandas DataFrame): omics dataframe to select column(s) from.
        genes (str, or list or array-like of str): list or array-like of genes, as strings, to use to select columns.

        Returns:
        pandas DataFrame: The selected columns from the dataframe.
        """
        df = pd.DataFrame(index=omics_df.index) # Create an empty dataframe, which we'll fill with the columns we select using our genes, and then return.
        for gene in genes:
            selected = self.get_col_from_omics(omics_df, gene) # Get the columns for that gene from the dataframe
            df = df.join(selected, how='left') # Otherwise, append the columns to our dataframe we'll return.
        df.name = "{} for {} genes".format(omics_df.name, len(genes)) # Name the dataframe!
        return df

    def get_omics_from_str_or_list(self, omics_df, genes):
        """Determines whether you passed it a single gene or a list or array-like of genes, selects the corresponding columns from the omics dataframe, and returns them.

        Parameters:
        omics_df (pandas DataFrame): omics dataframe to select columns from
        genes (str, or list or array-like of str): gene(s) to use to select columns from omics_df. str if one gene, list or array-like if multiple. Passing None will select the entire omics dataframe.

        Returns:
        pandas DataFrame: the columns from omics_df corresponding to the gene(s), as a dataframe.
        """
        if isinstance(genes, str): # If it's a single gene, feed it to the proper function
            return self.get_col_from_omics(omics_df, genes) 
        elif isinstance(genes, (list, pd.core.series.Series, pd.core.indexes.base.Index)): # If it's a str, or list or array-like of str of genes, feed it to the proper function
            return self.get_cols_from_omics(omics_df, genes)
        elif genes is None: # If it's the default of None, rename columns and return the entire dataframe
            return_df = omics_df.rename(columns=lambda x:'{}_{}'.format(x, omics_df.name)) # Append dataframe name to end of each column header, to preserve info when we merge dataframes
            return_df.name = omics_df.name # Name the return dataframe
            return return_df
        else: # If it's none of those, they done messed up. Tell 'em.
            print("Genes parameter \n{}\nis of invalid type {}. Valid types: str, list or array-like of str, or NoneType.".format(genes, type(genes)))

    def compare_omics(self, df1, df2, genes1, genes2):
        """Select columns for one gene or a list or array-like of genes from one omics dataframe, and columns for another gene or list or array-like of genes from another omics dataframe, and join them into one dataframe. Intersection (inner join) of indicies is used.

        Parameters:
        df1 (pandas DataFrame): first omics dataframe to select from.
        genes1 (str, or list or array-like of str): gene(s) for column(s) to select from the first omics dataframe. str if one gene, list or array-like of strings if multiple genes. Passing None will select the entire omics dataframe.
        df2 (pandas DataFrame): second omics dataframe to select from.
        genes2 (str, or list or array-like of str): gene(s) for column(s) to select from the second omics dataframe. str if one gene, list or array-like of strings if multiple genes. Passing None will select the entire omics dataframe.

        Returns:
        pandas DataFrame: The data from the selected columns from the two dataframes, joined into one dataframe.
        """
        selected1 = self.get_omics_from_str_or_list(df1, genes1)
        selected2 = self.get_omics_from_str_or_list(df2, genes2)

        if (selected1 is not None) and (selected2 is not None): # If either selector returned None, the gene(s) didn't match any columns, and it printed an informative error message already. We'll return None.
            df = selected1.join(selected2, how='inner') # Join the rows common to both dataframes
            df = df.sort_index() # Sort rows in ascending order
            df.name = "{}, with {}".format(selected1.name, selected2.name) # Give it a nice name identifying the data in it.
            return df

# Next 4 functions are for working with metadata dataframes
    def get_col_from_metadata(self, df, col):
        """Get a single column from a metadata dataframe.

        Parameters:
        df (pandas DataFrame): The dataframe to select the column from.
        col (str): The column to select from the dataframe.

        Returns:
        pandas DataFrame: The specified column from the given dataframe. We return it as a dataframe for easy merging later.
        """
        if col not in df.columns.values: # If they didn't give us one of the actual columns, tell them and return None.
            print('{} column not found in {} dataframe. Please double check that it is included in the dataframe.'.format(col, df.name))
            return
        
        selected = df.loc[:, [col]] # Select the column from the dataframe, keeping it as a dataframe
        selected.name = '{} from {}'.format(col, df.name) # Give it a name, identifying which dataframe it came from
        return selected

    def get_cols_from_metadata(self, df, cols):
        """Select several columns, given as a list, from a metadata dataframe.

        Parameters:
        df (pandas DataFrame): The dataframe to select the column from.
        cols (str, or list or array-like of str): A list or array-like of the columns, as strings, to select from the dataframe.

        Returns:
        pandas DataFrame: The specified columns from the given dataframe.
        """
        return_df = pd.DataFrame(index=df.index) # Create an empty dataframe, which we'll fill with the columns we select, and then return.
        for col in cols:
            selected = self.get_col_from_metadata(df, col) # Get the columns from the given dataframe
            if selected is None: # If it didn't match any columns, get_col_from_omics will have printed an error message. Return None.
                return
            return_df = return_df.join(selected, how='left') # Otherwise, append the columns to our dataframe we'll return.
        return_df.name = "{} columns from {}".format(len(cols), df.name) # Name the dataframe!
        return return_df

    def get_metadata_from_str_or_list(self, df, cols):
        """Determines whether you passed it a single column or a list or array-like of columns, then selects them from the given dataframe, and returns them.

        Parameters:
        df (pandas DataFrame): The dataframe to select column(s) from.
        cols (str, or list or array-like of str): Column(s) to select from the dataframe. str if one, list or array-like if multiple.

        Returns:
        pandas DataFrame: The specificed columns from the given dataframe.
        """
        if isinstance(cols, str): # If it's a single column, feed it to the proper function
            return self.get_col_from_metadata(df, cols)
        elif isinstance(cols, (list, pd.core.series.Series, pd.core.indexes.base.Index)): # If it's a str, or list or array-like of str of columns, feed it to the proper function
            return self.get_cols_from_metadata(df, cols)
        elif cols is None: # If it's the default of None, return the entire dataframe
            return df
        else: # If it's none of those, they done messed up. Tell 'em.
            print("Columns parameter {} is of invalid type {}. Valid types: str, or list or array-like of str.".format(cols, type(cols)))

    def append_metadata_to_omics(self, df, omics_df, df_cols, omics_genes):
        """Append all or part of the given metadata dataframe to all or part of the given omics dataframe. Intersection (inner join) of indicies is used.

        Parameters:
        df (pandas DataFrame): The metadata dataframe from which we'll select our columns to append.
        omics_df (pandas DataFrame): The omics dataframe to append to.
        df_cols (str, or list or array-like of str): Column(s) to select from the metadata dataframe. str if one column, list or array-like of strings if multiple.
        omics_genes (str, or list or array-like of str): Gene(s) to select data for from the omics dataframe. Passing None will select the entire omics dataframe.

        Returns:
        pandas DataFrame: The selected columns from the metadata dataframe, appended to the selected columns from the omics dataframe.
        """  
        df_selected = self.get_metadata_from_str_or_list(df, df_cols)
        omics_selected = self.get_omics_from_str_or_list(omics_df, omics_genes)

        if (df_selected is not None) and (omics_selected is not None): # If either selector returned None, the key(s) didn't match any columns, and it printed an informative error message already. We'll return None.
            df_joined = df_selected.join(omics_selected, how='inner') # Join the rows common to both dataframes
            df_joined = df_joined.sort_index() # Sort rows in ascending order
            df_joined.name = "{}, with {}".format(df_selected.name, omics_selected.name) # Give it a nice name identifying the data in it.
            return df_joined

# Next 4 functions are for working with mutation data
    def get_mutations_for_gene(self, mutations, gene):
        """Gets all the mutations for a specific gene, for all patients.

        Parameters:
        mutations (pandas DataFrame): The somatic_mutation dataframe, which we'll grab the mutation data from.
        gene (str): The gene to grab mutations for.

        Returns:
        pandas DataFrame: The mutations in each patient for the specified gene.
        """
        # Set some column names for use later
        gene_col = "Gene"
        mutation_col = "Mutation"
        location_col = "Location"
        mutation_status_col = "Mutation_Status"

        gene_mutations = mutations[mutations[gene_col] == gene]

        if len(gene_mutations) == 0: # If the gene doesn't match any in the dataframe, tell them, and return None.
            print("{} gene not found in somatic_mutation data.".format(gene))
            return
        
        gene_mutations = gene_mutations.drop(columns=[gene_col]) # Gene column is same for every sample, so we don't need it anymore.
        
        # Create an empty dataframe, which we'll fill with the mutation and location data as lists
        prep_index = gene_mutations.index.copy().drop_duplicates()
        prep_columns = gene_mutations.columns
        prep_cols_with_mut_status = prep_columns.union([mutation_status_col], sort=False) # Add a mutation_status column, which will indicate if there are 1 or multiple mutations
        mutation_lists = pd.DataFrame(index=prep_index, columns=prep_cols_with_mut_status)

        # Get the mutation(s), mutation status, and location information for this gene and sample
        for sample in mutation_lists.index: 
            sample_data = gene_mutations.loc[sample] # Get slice of dataframe for the sample
            sample_mutations = sample_data[mutation_col] # Get mutation(s)
            sample_locations = sample_data[location_col] # Get location(s)

            # Make the mutations a list (even if there's only one)
            if isinstance(sample_mutations, pd.core.series.Series):
                sample_mutations_list = sample_mutations.tolist()
            else:
                sample_mutations_list = [sample_mutations]

            # Make the locations a list (even if there's only one)
            if isinstance(sample_locations, pd.core.series.Series):
                sample_locations_list = sample_locations.tolist()
            else:
                sample_locations_list = [sample_locations]

            # Figure out what our mutation status is (either single_mutation or multiple_mutation)
            if len(sample_mutations_list) > 1:
                sample_mutation_status = "Multiple_mutation"
            else:
                sample_mutation_status = "Single_mutation"

            # Put in our template dataframe
            mutation_lists.at[sample, mutation_col] = sample_mutations_list
            mutation_lists.at[sample, location_col] = sample_locations_list
            mutation_lists.at[sample, mutation_status_col] = sample_mutation_status

        mutation_lists = mutation_lists.rename(columns=lambda x:'{}_{}'.format(gene, x)) # Add the gene name to end beginning of each column header, to preserve info when we merge dataframes.
        mutation_lists.name = 'Somatic mutation data for {} gene'.format(gene)
        return mutation_lists

    def get_mutations_for_genes(self, somatic_mutation, genes):
        """Gets all the mutations for a list or array-like of genes, for all patients.

        Parameters:
        somatic_mutation (pandas DataFrame): The somatic_mutation dataframe that we'll grab the mutation data from.
        genes (str, or list or array-like of str): The genes, as strings, to grab mutations for.

        Returns:
        pandas DataFrame: The mutations in each patient for the specified genes.
        """
        df = pd.DataFrame(index=somatic_mutation.index.drop_duplicates()) # Create an empty dataframe, which we'll fill with the columns we select using our genes, and then return.
        for gene in genes:
            selected = self.get_mutations_for_gene(somatic_mutation, gene) # Get the mutations for our gene
            if selected is None: # If there's no mutation data for that gene, get_mutations_for_gene will have printed an error message. Return None.
                return
            df = df.join(selected, how='left') # Otherwise, append the columns to our dataframe we'll return.
        df.name = "Somatic mutation data for {} genes".format(len(genes)) # Name the dataframe!
        return df

    def get_mutations_from_str_or_list(self, somatic_mutation, genes):
        """Determines whether you passed it a single gene or a list (or pandas.core.series.Series or pandas.core.indexes.base.Index of genes), selects the corresponding mutations from the somatic_mutation dataframe, and returns them.

        Parameters:
        somatic_mutation (pandas DataFrame): The somatic_mutation dataframe we'll grab the mutation data from.
        genes (str, or list or array-like of str): gene(s) to select mutations for. str if one gene, list or array-like if multiple.

        Returns:
        pandas DataFrame: the mutation data corresponding to the gene(s), as a dataframe.
        """
        if isinstance(genes, str): # If it's a single gene, feed it to the proper function
            return self.get_mutations_for_gene(somatic_mutation, genes) 
        elif isinstance(genes, (list, pd.core.series.Series, pd.core.indexes.base.Index)): # If it's a str, or list or array-like of str of genes, feed it to the proper function
            return self.get_mutations_for_genes(somatic_mutation, genes)
        else: # If it's neither of those, they done messed up. Tell 'em.
            print("Genes parameter {} is of invalid type {}. Valid types: str, or list or array-like of str.".format(genes, type(genes)))

    def append_mutations_to_omics(self, somatic_mutation, omics_df, mutation_genes, omics_genes, show_location, sample_status_map):
        """Select all mutations for specified gene(s), and append to all or part of the given omics dataframe. Intersection (inner join) of indicies is used. Each location or mutation cell contains a list, which contains the one or more location or mutation values corresponding to that sample for that gene, or a value indicating that the sample didn't have a mutation in that gene.

        Parameters:
        somatic_mutation (pandas DataFrame): Somatic mutation dataframe we'll get the mutation data from.
        omics_df (pandas DataFrame): Omics dataframe to append the mutation data to.
        mutation_genes (str, or list or array-like of str): The gene(s) to get mutation data for. str if one gene, list or array-like if multiple.
        omics_genes (str, or list or array-like of str): Gene(s) to select from the omics dataframe. str if one gene, list or array-like if multiple. Passing None will select the entire omics dataframe.
        show_location (bool): Whether to include the Location column from the mutation dataframe.
        sample_status_map (pandas Series): A series with the dataset's sample IDs as the indicies, and each sample's status (tumor or normal) as the values. Used to generate the Sample_Status column.

        Returns:
        pandas DataFrame: The mutations for the specified gene, appended to all or part of the omics dataframe. Each location or mutation cell contains a list, which contains the one or more location or mutation values corresponding to that sample for that gene, or a value indicating that the sample didn't have a mutation in that gene.
        """
        omics = self.get_omics_from_str_or_list(omics_df, omics_genes)
        mutations = self.get_mutations_from_str_or_list(somatic_mutation, mutation_genes)

        if (omics is not None) and (mutations is not None): # If either selector returned None, then there were gene(s) that didn't match anything, and an error message was printed. We'll return None.
            merge = omics.join(mutations, how = "left") # Left join omics data and mutation data (left being the omics data)

            # Add Sample_Status column by joining the sample_status_map to the merged mutation dataframe. Do a left join so we drop any indicies not in the mutations dataframe.
            merge = merge.join(sample_status_map, how="left") 

            # Fill in Wildtype_Normal or Wildtype_Tumor for NaN values (i.e., no mutation data for that sample) in merged dataframe mutation columns
            mutation_regex = r'.*_Mutation$' # Construct regex to find all mutation columns
            mutation_cols = [col for col in merge.columns.values if re.match(mutation_regex, col)] # Get a list of all mutation columns
            for mutation_col in mutation_cols:
                merge.loc[(merge['Sample_Status'] == "Normal") & (pd.isnull(merge[mutation_col])), mutation_col] = [[["Wildtype_Normal"]]] # Change all NaN mutation values for Normal samples to Wildtype_Normal. Triple nested list causes .loc to insert the value as ['Wildtype_Normal'], like we want it to, instead of unpacking the list.
                merge.loc[(merge['Sample_Status'] == "Tumor") & (pd.isnull(merge[mutation_col])), mutation_col] = [[["Wildtype_Tumor"]]] # Change all NaN mutation values for Tumor samples to Wildtype_Tumor

            # Depending on show_location, either fill NaN values in the merged dataframe location columns with "No_mutation", or just drop the location columns altogether
            location_regex = r'.*_Location$' # Construct regex to find all location columns
            location_cols = [col for col in merge.columns.values if re.match(location_regex, col)] # Get a list of all location columns
            for location_col in location_cols:
                if show_location:
                    merge.loc[pd.isnull(merge[location_col]), location_col] = [[["No_mutation"]]] # If there's no location, there wasn't a mutation--make it easier for people to understand that.
                else:
                    merge = merge.drop(columns=[location_col]) # Drop the location column, if the caller wanted us to.

            # Fill NaN values in Mutation_Status column with either Wildtype_Tumor or Wildtype_Normal
            mutation_status_regex = r".*_Mutation_Status$" # Construct a regex to find all Mutation_Status columns
            mutation_status_cols = [col for col in merge.columns.values if re.match(mutation_status_regex, col)] # Get a list of all Mutation_Status columns
            for mutation_status_col in mutation_status_cols:
                merge.loc[(merge['Sample_Status'] == "Normal") & (pd.isnull(merge[mutation_status_col])), mutation_status_col] = "Wildtype_Normal" # Change all NaN mutation status values for Normal samples to Wildtype_Normal
                merge.loc[(merge['Sample_Status'] == "Tumor") & (pd.isnull(merge[mutation_status_col])), mutation_status_col] = "Wildtype_Tumor" # Change all NaN mutation status values for Tumor samples to Wildtype_Tumor

            merge.name = "{}, with {}".format(omics.name, mutations.name) # Give it a name identifying the data in it
            return merge
