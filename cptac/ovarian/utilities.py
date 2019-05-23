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

# Next 2 functions are for working with omics data
    def get_omics_cols(self, omics_df, genes):
        """Based on a single gene, or a list or array-like of genes, select multiple columns from an omics dataframe, and return the selected columns as one dataframe.

        Parameters:
        omics_df (pandas DataFrame): Omics dataframe to select column(s) from.
        genes (str, or list or array-like of str): Gene(s) to use to select columns from omics_df. str if one gene, list or array-like if multiple. Passing None will select the entire omics dataframe.

        Returns:
        pandas DataFrame: The selected columns from the dataframe.
        """
        # Process genes parameter
        if isinstance(genes, str): # If it's a single gene, make it a list so we can treat everything the same
            genes = [genes]
        elif isinstance(genes, (list, pd.core.series.Series, pd.core.indexes.base.Index)): # If it's already a list or array-like, we're all good
            pass
        elif genes is None: # If it's the default of None, rename columns and return the entire dataframe
            return_df = omics_df.rename(columns=lambda x:'{}_{}'.format(x, omics_df.name)) # Append dataframe name to end of each column header, to preserve info when we merge dataframes
            return_df.name = omics_df.name # Name the return dataframe
            return return_df
        else: # If it's none of those, they done messed up. Tell 'em.
            print("Genes parameter \n{}\nis of invalid type {}. Valid types: str, list or array-like of str, or NoneType.".format(genes, type(genes)))
            return

        df = pd.DataFrame(index=omics_df.index.copy()) # Create an empty dataframe, which we'll fill with the columns we select using our genes, and then return.
        for gene in genes:
            if omics_df.name == 'phosphoproteomics' or omics_df.name == 'acetylproteomics':
                col_regex = "^{}-.*$".format(gene) # Build a regex to get all columns that match the gene
            else:
                col_regex = '^{}$'.format(gene)

            selected = omics_df.filter(regex=col_regex) # Find all columns that match the gene. 
            if len(selected.columns) == 0: # If none of the columns matched the gene, generate a column of NaN and print a warning message
                empty_omics_df = pd.DataFrame(index=omics_df.index.copy())
                selected = empty_omics_df.assign(**{gene:np.nan}) # Create a column with gene as the name, filled with NaN
                print('{0} did not match any columns in {1} dataframe. {0}_{1} column inserted, but filled with NaN.'.format(gene, omics_df.name))

            selected = selected.rename(columns=lambda x:'{}_{}'.format(x, omics_df.name)) # Append dataframe name to end of each column header, to preserve info when we merge dataframes
            df = df.join(selected, how='left') # Append the columns to our dataframe we'll return.

        # Give the dataframe a name!
        if len(genes) == 1:
            df.name = "{} for {}".format(omics_df.name, genes[0]) 
        else:
            df.name = "{} for {} genes".format(omics_df.name, len(genes)) 
        return df

    def compare_omics(self, df1, df2, genes1, genes2):
        """Select columns for one gene or a list or array-like of genes from one omics dataframe, and columns for another gene or list or array-like of genes from another omics dataframe, and join them into one dataframe. Intersection (inner join) of indicies is used.

        Parameters:
        df1 (pandas DataFrame): First omics dataframe to select from.
        genes1 (str, or list or array-like of str): Gene(s) for column(s) to select from the first omics dataframe. str if one gene, list or array-like of strings if multiple genes. Passing None will select the entire omics dataframe.
        df2 (pandas DataFrame): Second omics dataframe to select from.
        genes2 (str, or list or array-like of str): Gene(s) for column(s) to select from the second omics dataframe. str if one gene, list or array-like of strings if multiple genes. Passing None will select the entire omics dataframe.

        Returns:
        pandas DataFrame: The data from the selected columns from the two dataframes, joined into one dataframe.
        """
        selected1 = self.get_omics_cols(df1, genes1)
        selected2 = self.get_omics_cols(df2, genes2)

        if (selected1 is not None) and (selected2 is not None): # If either selector returned None, the gene(s) didn't match any columns, and it printed an informative error message already. We'll return None.
            df = selected1.join(selected2, how='inner') # Join the rows common to both dataframes
            df = df.sort_index() # Sort rows in ascending order
            df.name = "{}, with {}".format(selected1.name, selected2.name) # Give it a nice name identifying the data in it.
            return df

# Next 2 functions are for working with metadata dataframes
    def get_metadata_cols(self, df, cols):
        """Select a single column or several columns from a metadata dataframe.

        Parameters:
        df (pandas DataFrame): The dataframe to select the column(s) from.
        cols (str, or list or array-like of str): The column(s) to select from the dataframe. str if single, list or array-like of str if multiple. Passing None will select the entire dataframe.

        Returns:
        pandas DataFrame: The specified columns from the given dataframe.
        """
        # Process genes parameter
        if isinstance(cols, str): # If it's a single column, make it a list so we can treat everything the same
            cols = [cols]
        elif isinstance(cols, (list, pd.core.series.Series, pd.core.indexes.base.Index)): # If it's already a list or array-like, we're all good
            pass
        elif cols is None: # If it's the default of None, return the entire dataframe
            return df
        else: # If it's none of those, they done messed up. Tell 'em.
            print("Columns parameter {} is of invalid type {}. Valid types: str, or list or array-like of str.".format(cols, type(cols)))
            return

        return_df = pd.DataFrame(index=df.index.copy()) # Create an empty dataframe, which we'll fill with the columns we select, and then return.
        for col in cols:
            if col not in df.columns.values: # If they didn't give us one of the actual columns, tell them and return None.
                print('{} column not found in the {} dataframe. Please double check that it is included in the dataframe.'.format(col, df.name))
                return
            selected = df.loc[:, [col]] # Select the column from the dataframe, keeping it as a dataframe
            return_df = return_df.join(selected, how='left') # Append the columns to our dataframe we'll return.

        # Give return_df a name, identifying which dataframe it came from
        if len(cols) == 1:
            return_df.name = '{} from {}'.format(cols[0], df.name) 
        else:
            return_df.name = "{} columns from {}".format(len(cols), df.name) 
        return return_df

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
        df_selected = self.get_metadata_cols(df, df_cols)
        omics_selected = self.get_omics_cols(omics_df, omics_genes)

        if (df_selected is not None) and (omics_selected is not None): # If either selector returned None, the key(s) didn't match any columns, and it printed an informative error message already. We'll return None.
            df_joined = df_selected.join(omics_selected, how='inner') # Join the rows common to both dataframes
            df_joined = df_joined.sort_index() # Sort rows in ascending order
            df_joined.name = "{}, with {}".format(df_selected.name, omics_selected.name) # Give it a nice name identifying the data in it.
            return df_joined

# Next 2 functions are for working with mutation data
    def get_genes_mutations(self, somatic_mutation, genes):
        """Gets all the mutations for one or multiple genes, for all patients.

        Parameters:
        somatic_mutation (pandas DataFrame): The somatic_mutation dataframe that we'll grab the mutation data from.
        genes (str, or list or array-like of str): The gene(s) to grab mutations for. str if one, list or array-like of str if multiple.

        Returns:
        pandas DataFrame: The mutations in each patient for the specified gene(s).
        """
        # Process genes parameter
        if isinstance(genes, str): # If it's a single gene, make it a list so we can treat everything the same
            genes = [genes]
        elif isinstance(genes, (list, pd.core.series.Series, pd.core.indexes.base.Index)): # If it's already a list or array-like, we're all good
            pass
        else: # If it's neither of those, they done messed up. Tell 'em.
            print("Genes parameter {} is of invalid type {}. Valid types: str, or list or array-like of str.".format(genes, type(genes)))
            return

        # Set some column names for use later
        gene_col = "Gene"
        mutation_col = "Mutation"
        location_col = "Location"
        mutation_status_col = "Mutation_Status"

        # Create an empty dataframe, which we'll fill with the columns we select using our genes, and then return.
        df = pd.DataFrame(index=somatic_mutation.index.copy().drop_duplicates())
        for gene in genes:
            gene_mutations = somatic_mutation[somatic_mutation[gene_col] == gene] # Get all the mutations for that gene
            if len(gene_mutations) == 0: # If the gene doesn't match any genes in the dataframe, tell them, and return None.
                print("{} gene not found in somatic_mutation data.".format(gene))
                return
            gene_mutations = gene_mutations.drop(columns=[gene_col]) # Gene column is same for every sample, so we don't need it anymore.
            
            # Create another empty dataframe, which we'll fill with the mutation and location data for this gene, as lists
            prep_index = gene_mutations.index.copy().drop_duplicates()
            prep_columns = gene_mutations.columns.copy()
            mutation_status_idx = pd.Index([mutation_status_col]) # Prep mutation_status_col to be appended
            prep_cols_with_mut_status = prep_columns.append(mutation_status_idx) # Add a mutation_status column, which will indicate if there are 1 or multiple mutations
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

                # Fill our template dataframe
                mutation_lists.at[sample, mutation_col] = sample_mutations_list
                mutation_lists.at[sample, location_col] = sample_locations_list
                mutation_lists.at[sample, mutation_status_col] = sample_mutation_status

            mutation_lists = mutation_lists.rename(columns=lambda x:'{}_{}'.format(gene, x)) # Add the gene name to end beginning of each column header, to preserve info when we merge dataframes.
            df = df.join(mutation_lists, how='left') # Append the columns to our dataframe we'll return.

        # Name the dataframe!
        if len(genes) == 1:
            df.name = 'somatic mutation data for {} gene'.format(genes[0])
        else:
            df.name = "somatic mutation data for {} genes".format(len(genes)) 
        return df

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
        omics = self.get_omics_cols(omics_df, omics_genes)
        mutations = self.get_genes_mutations(somatic_mutation, mutation_genes)

        if (omics is not None) and (mutations is not None): # If either selector returned None, then there were gene(s) that didn't match anything, and an error message was printed. We'll return None.
            merge = omics.join(mutations, how = "left") # Left join omics data and mutation data (left being the omics data)

            # Add Sample_Status column by joining the sample_status_map to the merged mutation dataframe. Do a left join so we drop any indicies not in the mutations dataframe.
            merge = merge.join(sample_status_map, how="left") 

            # Fill in Wildtype_Normal or Wildtype_Tumor for NaN values (i.e., no mutation data for that sample) in merged dataframe mutation columns
            mutation_regex = r'^.*_Mutation$' # Construct regex to find all mutation columns
            mutation_cols = [col for col in merge.columns.values if re.match(mutation_regex, col)] # Get a list of all mutation columns
            for mutation_col in mutation_cols:
                merge.loc[(merge['Sample_Status'] == "Normal") & (pd.isnull(merge[mutation_col])), mutation_col] = [[["Wildtype_Normal"]]] # Change all NaN mutation values for Normal samples to Wildtype_Normal. Triple nested list causes .loc to insert the value as ['Wildtype_Normal'], like we want it to, instead of unpacking the list.
                merge.loc[(merge['Sample_Status'] == "Tumor") & (pd.isnull(merge[mutation_col])), mutation_col] = [[["Wildtype_Tumor"]]] # Change all NaN mutation values for Tumor samples to Wildtype_Tumor

            # Depending on show_location, either fill NaN values in the merged dataframe location columns with "No_mutation", or just drop the location columns altogether
            location_regex = r'^.*_Location$' # Construct regex to find all location columns
            location_cols = [col for col in merge.columns.values if re.match(location_regex, col)] # Get a list of all location columns
            for location_col in location_cols:
                if show_location:
                    merge.loc[pd.isnull(merge[location_col]), location_col] = [[["No_mutation"]]] # If there's no location, there wasn't a mutation--make it easier for people to understand that.
                else:
                    merge = merge.drop(columns=[location_col]) # Drop the location column, if the caller wanted us to.

            # Fill NaN values in Mutation_Status column with either Wildtype_Tumor or Wildtype_Normal
            mutation_status_regex = r"^.*_Mutation_Status$" # Construct a regex to find all Mutation_Status columns
            mutation_status_cols = [col for col in merge.columns.values if re.match(mutation_status_regex, col)] # Get a list of all Mutation_Status columns
            for mutation_status_col in mutation_status_cols:
                merge.loc[(merge['Sample_Status'] == "Normal") & (pd.isnull(merge[mutation_status_col])), mutation_status_col] = "Wildtype_Normal" # Change all NaN mutation status values for Normal samples to Wildtype_Normal
                merge.loc[(merge['Sample_Status'] == "Tumor") & (pd.isnull(merge[mutation_status_col])), mutation_status_col] = "Wildtype_Tumor" # Change all NaN mutation status values for Tumor samples to Wildtype_Tumor

            merge.name = "{}, with {}".format(omics.name, mutations.name) # Give it a name identifying the data in it
            return merge
