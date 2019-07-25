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
import webbrowser
import re

class DataSet:
    """
    Note that all cancer datasets are class objects that inherit from cptac.dataset. Therefore
    the same function calls exist for cptac.Endometrial, cptac.Colon, etc.
    """

    def __init__(self, cancer_type):
        """Initialize variables for a DataSet object."""

        # Initialize dataframe and definitions dicts as empty for this parent class
        self._data = {}
        self._definitions = {}

        # Initialize the _cancer_type instance variable
        self._cancer_type = cancer_type.lower()

        # Initialize the _version instance variable
        self._version = None

        # Assign the gene separator for searching columns of phosphoproteomics and acetylproteomics dataframes. Child class can override if needed.
        self._gene_separator = "-"

        # Assign the valid dfs lists, but make them instance variables so they're easy to override if needed
        # These are the omics dataframes that are valid for use in the utilities functions
        self._valid_omics_dfs = [
            'acetylproteomics',
            'proteomics',
            'transcriptomics', # But not circular_RNA or miRNA--they have incompatible column names.
            'CNV',
            'phosphoproteomics',
            'phosphoproteomics_gene']

        # These are the metadata dataframes that are valid for use in the utilities functions
        self._valid_metadata_dfs = [
            'clinical',
            'derived_molecular',
            'experimental_setup'] # We don't allow the treatment df, as in Ovarian, or medical_history df, as in RenalCcrcc, because they both have multiple rows for each sample.

    # Methods to get metadata dataframes
    def get_clinical(self):
        """Get the clinical dataframe."""
        return self._get_dataframe("clinical")

    def get_derived_molecular(self):
        """Get the derived_molecular dataframe."""
        return self._get_dataframe("derived_molecular")

    def get_experimental_setup(self):
        """Get the experimental_setup dataframe."""
        return self._get_dataframe("experimental_setup")

    def get_medical_history(self):
        """Get the medical_history dataframe."""
        return self._get_dataframe("medical_history")

    def get_treatment(self):
        """Get the treatment dataframe."""
        return self._get_dataframe("treatment")

    # Methods to get omics dataframes
    def get_acetylproteomics(self):
        """Get the acetylproteomics dataframe."""
        return self._get_dataframe("acetylproteomics")

    def get_circular_RNA(self):
        """Get the circular_RNA dataframe."""
        return self._get_dataframe("circular_RNA")

    def get_CNV(self):
        """Get the CNV dataframe."""
        return self._get_dataframe("CNV")

    def get_methylation(self):
        """Get the methylation dataframe."""
        return self._get_dataframe("methylation")

    def get_miRNA(self):
        """Get the miRNA dataframe."""
        return self._get_dataframe("miRNA")

    def get_phosphoproteomics(self):
        """Get the phosphoproteomics dataframe."""
        return self._get_dataframe("phosphoproteomics")

    def get_phosphoproteomics_gene(self):
        """Get the phosphoproteomics_gene dataframe. The gene level phosphorylation measurement is an aggregate metric which potentially averages together individual measurements of different sites. Use get_phosphoproteomics() to view the data for individual sites."""
        return self._get_dataframe("phosphoproteomics_gene")

    def get_phosphosites(self, genes):
        """Returns dataframe with all phosphosites of specified gene or list of genes.

        Parameters:
        genes (str, or list or array-like of str): gene or list of genes to use to select phosphosites. str if single, list or array-like of str if multiple.

        Returns:
        pandas DataFrame: The phosphoproteomics for the specified gene(s).
        """
        return self._get_omics_cols("phosphoproteomics", genes)

    def get_proteomics(self):
        """Get the proteomics dataframe."""
        return self._get_dataframe("proteomics")

    def get_transcriptomics(self):
        """Get the transcriptomics dataframe."""
        return self._get_dataframe("transcriptomics")

    # Methods to get mutations dataframes
    def get_mutations(self):
        """Get the somatic_mutation dataframe."""
        return self._get_dataframe("somatic_mutation")

    def get_mutations_binary(self):
        """Get the somatic_mutation_binary dataframe, which has a binary value indicating, for each location on each gene, whether there was a mutation in that gene at that location, for each sample."""
        return self._get_dataframe("somatic_mutation_binary")

    # Help methods
    def define(self, term):
        """Print the definition a term, if it is in the dataset's list of definitions.

        Parameters:
        term (str): term to be defined

        Returns: None
        """
        if term in self._definitions.keys():
            print(self._definitions[term])
        else:
            print("{} not found in definitions. Check capitalization. Alternatively, the dataset's search(term) method can be used to perform a web search of the term provided.".format(term))

    def get_cancer_type(self):
        """Return the cancer type for this dataset, as a string."""
        return self._cancer_type

    def get_version(self):
        """Return the dataset version of this instance, as a string."""
        return self._version

    def how_to_cite(self):
        """Print instructions for citing the data."""
        print('Please include the following statement in publications using data accessed through this module:\n"Data used in this publication were generated by the Clinical Proteomic Tumor Analysis Consortium (NCI/NIH, <https://proteomics.cancer.gov/programs/cptac/>). Data were accessed through the Python module cptac, available at <https://pypi.org/project/cptac/>."')

    def list_data(self):
        """Print list of loaded dataframes and dimensions."""
        print("Below are the dataframes contained in this dataset:")
        for name in sorted(self._data.keys(), key=str.lower):
            df = self._data[name]
            print("\t{}\n\t\tDimensions: {}".format(name, df.shape))

    def list_definitions(self):
        """Print all terms defined in the dataset's list of definitions."""
        if len(self._definitions.keys()) > 0:
            for term in sorted(self._definitions.keys(), key=str.lower):
                print(term)
        else:
            print("No definitions provided for this dataset.")

    def search(self, term):
        """Search for a term in a web browser.

        Parameters:
        term (str): term to be searched

        Returns: None
        """
        url = "https://www.google.com/search?q=" + term
        print("Searching for {} in web browser...".format(term))
        webbrowser.open(url)

    # Utilities methods
    def join_omics_to_omics(self, df1_name, df2_name, genes1=None, genes2=None):
        """Take specified column(s) from one omics dataframe, and join to specified columns(s) from another omics dataframe. Intersection (inner join) of indices is used.

        Parameters:
        df1_name (str): Name of first omics dataframe to select columns from.
        df2_name (str): Name of second omics dataframe to select columns from.
        genes1 (str, or list or array-like of str, optional): Gene(s) for column(s) to select from df1_name. str if one key, list or array-like of str if multiple. Default of None will select entire dataframe.
        genes2 (str, or list or array-like of str, optional): Gene(s) for Column(s) to select from df2_name. str if one key, list or array-like of str if multiple. Default of None will select entire dataframe.

        Returns:
        pandas DataFrame: The selected columns from the two omics dataframes, joined into one dataframe.
        """
        # Make sure they gave us valid dataframe names
        df1_valid = self._is_valid_df(df1_name, "omics")
        df2_valid = self._is_valid_df(df2_name, "omics")

        if (not df1_valid) or (not df2_valid):
            return

        # Select the columns from each dataframe
        selected1 = self._get_omics_cols(df1_name, genes1)
        selected2 = self._get_omics_cols(df2_name, genes2)

        if (selected1 is not None) and (selected2 is not None): # If either selector returned None, the gene(s) didn't match any columns, or a non-included dataframe, and it printed an informative error message already. We'll return None.
            df = selected1.join(selected2, how='inner') # Join the rows common to both dataframes
            df = df.sort_index() # Sort rows in ascending order
            return df

    def join_omics_to_mutations(self, omics_df_name, mutations_genes, omics_genes=None, mutations_filter=None, show_location=True):
        """Select all mutations for specified gene(s), and joins them to all or part of the given omics dataframe. Intersection (inner join) of indices is used. Each location or mutation cell contains a list, which contains the one or more location or mutation values corresponding to that sample for that gene, or a value indicating that the sample didn't have a mutation in that gene.

        Parameters:
        omics_df (str): Name of omics dataframe to join the mutation data to.
        mutations_genes (str, or list or array-like of str): The gene(s) to get mutation data for. str if one gene, list or array-like of str if multiple.
        omics_genes (str, or list or array-like of str, optional): Gene(s) to select from the omics dataframe. str if one gene, list or array-like of str if multiple. Default will select entire dataframe.
        mutations_filter (list, optional): List of mutations to prioritize when filtering out multiple mutations, in order of priority. If none of the multiple mutations in a sample are included in filter_prefer, the function will automatically prioritize truncation over missense mutations, and then mutations earlier in the sequence over later mutations. Passing an empty list will cause this default hierarchy to be applied to all samples. Default parameter of None will cause no filtering to be done, and all mutation data will be included, in a list.
        show_location (bool, optional): Whether to include the Location column from the mutation dataframe. Defaults to True.

        Returns:
        pandas DataFrame: The mutations for the specified gene, joined to all or part of the omics dataframe. Each location or mutation cell contains a list, which contains the one or more location or mutation values corresponding to that sample for that gene, or a value indicating that the sample didn't have a mutation in that gene.
        """
        # Make sure omics_df is valid for this function
        if not (self._is_valid_df(omics_df_name, "omics")):
            return

        # Select the data from each dataframe
        omics = self._get_omics_cols(omics_df_name, omics_genes)
        mutations = self._get_genes_mutations(mutations_genes, mutations_filter)

        if (omics is not None) and (mutations is not None): # If either selector returned None, then there were gene(s) that didn't match anything, or a non-included dataframe, and an error message was printed. We'll return None.
            mutations_were_filtered = mutations_filter is not None
            joined = self._join_other_to_mutations(omics, mutations, mutations_were_filtered, show_location)
            return joined

    def join_metadata_to_metadata(self, df1_name, df2_name, cols1=None, cols2=None):
        """Take specified column(s) from one metadata dataframe, and join to specified columns(s) from another metadata dataframe. Intersection (inner join) of indices is used.

        Parameters:
        df1_name (str): Name of first metadata dataframe to select columns from.
        df2_name (str): Name of second metadata dataframe to select columns from.
        cols1 (str, or list or array-like of str, optional): Column(s) to select from df1_name. str if one key, list or array-like of str if multiple. Default of None will select entire dataframe.
        cols2 (str, or list or array-like of str, optional): Column(s) to select from df2_name. str if one key, list or array-like of str if multiple. Default of None will select entire dataframe.

        Returns:
        pandas DataFrame: The selected columns from the two metadata dataframes, joined into one dataframe.
        """
        # Make sure they gave us valid dataframe names
        df1_valid = self._is_valid_df(df1_name, "metadata")
        df2_valid = self._is_valid_df(df2_name, "metadata")

        if (not df1_valid) or (not df2_valid):
            return

        # Select the columns from each dataframe
        selected1 = self._get_metadata_cols(df1_name, cols1)
        selected2 = self._get_metadata_cols(df2_name, cols2)

        if (selected1 is not None) and (selected2 is not None): # If either selector returned None, the gene(s) didn't match any columns, or a non-included dataframe, and it printed an informative error message already. We'll return None.
            df = selected1.join(selected2, how='inner') # Join the rows common to both dataframes
            df = df.sort_index() # Sort rows in ascending order
            return df

    def join_metadata_to_omics(self, metadata_df_name, omics_df_name, metadata_cols=None, omics_genes=None):
        """Joins columns from a metadata dataframe (clinical, derived_molecular, or experimental_setup) to part or all of an omics dataframe. Intersection (inner join) of indices is used.

        Parameters:
        metadata_df_name (str): Name of metadata dataframe to select columns from.
        omics_df_name (str): Name of omics dataframe to join the metadata columns to.
        metadata_cols (str, or list or array-like of str, optional): Column(s) to select from the metadata dataframe. str if one gene, list or array-like of str if multiple. Default is None, which will select the entire metadata dataframe.
        omics_genes (str, or list or array-like of str, optional): Gene(s) to select data for from the omics dataframe. str if one gene, list or array-like of str if multiple. Default is None, which will select entire dataframe.

        Returns:
        pandas DataFrame: The selected metadata columns, joined with all or part of the omics dataframe.
        """
        # Make sure metadata_df_name and omics_df_name are valid for this function
        metadata_df_valid = self._is_valid_df(metadata_df_name, "metadata")
        omics_df_valid = self._is_valid_df(omics_df_name, "omics")

        if (not metadata_df_valid) or (not omics_df_valid):
            return

        # Select the columns from each dataframe
        metadata_selected = self._get_metadata_cols(metadata_df_name, metadata_cols)
        omics_selected = self._get_omics_cols(omics_df_name, omics_genes)

        if (metadata_selected is not None) and (omics_selected is not None): # If either selector returned None, the key(s) didn't match any columns, or a non-included dataframe, and it printed an informative error message already. We'll return None.
            joined = metadata_selected.join(omics_selected, how='inner') # Join the rows common to both dataframes
            joined = joined.sort_index() # Sort rows in ascending order
            return joined

    def join_metadata_to_mutations(self, metadata_df_name, mutations_genes, metadata_cols=None, mutations_filter=None, show_location=True):
        """Select all mutations for specified gene(s), and joins them to all or part of the given metadata dataframe. Intersection (inner join) of indices is used. Each location or mutation cell contains a list, which contains the one or more location or mutation values corresponding to that sample for that gene, or a value indicating that the sample didn't have a mutation in that gene.

        Parameters:
        metadata_df_name (str): Name of metadata dataframe to join the mutation data to.
        mutations_genes (str, or list or array-like of str): The gene(s) to get mutation data for. str if one gene, list or array-like of str if multiple.
        metadata_cols (str, or list or array-like of str, optional): Gene(s) to select from the metadata dataframe. str if one gene, list or array-like of str if multiple. Default will select entire dataframe.
        mutations_filter (list, optional): List of mutations to prioritize when filtering out multiple mutations, in order of priority. If none of the multiple mutations in a sample are included in filter_prefer, the function will automatically prioritize truncation over missense mutations, and then mutations earlier in the sequence over later mutations. Passing an empty list will cause this default hierarchy to be applied to all samples. Default parameter of None will cause no filtering to be done, and all mutation data will be included, in a list.
        show_location (bool, optional): Whether to include the Location column from the mutation dataframe. Defaults to True.

        Returns:
        pandas DataFrame: The mutations for the specified gene, joined to all or part of the metadata dataframe. Each location or mutation cell contains a list, which contains the one or more location or mutation values corresponding to that sample for that gene, or a value indicating that the sample didn't have a mutation in that gene.
        """
        # Make sure metadata_df is valid for this function
        if not (self._is_valid_df(metadata_df_name, "metadata")):
            return

        # Select the data from each dataframe
        metadata = self._get_metadata_cols(metadata_df_name, metadata_cols)
        mutations = self._get_genes_mutations(mutations_genes, mutations_filter)

        if (metadata is not None) and (mutations is not None): # If either selector returned None, then there were gene(s) that didn't match anything, or a non-included dataframe, and an error message was printed. We'll return None.
            mutations_were_filtered = mutations_filter is not None
            joined = self._join_other_to_mutations(metadata, mutations, mutations_were_filtered, show_location)
            return joined

    # "Private" methods
    def _get_dataframe(self, name):
        """Check if a dataframe with the given name exists, and return a copy of it if it does.

        Parameters:
        name (str): The name of the dataframe to get.

        Returns:
        pandas DataFrame: A copy of the desired dataframe, if it exists in this dataset.
        """
        if name in self._data.keys():
            df = self._data[name]
            return_df = df.copy() # We copy it, with default deep=True, so edits on their copy don't affect the master
            return return_df
            
        else:
            print("{} dataframe not included in this dataset.".format(name))
            return

    def _get_sample_status_map(self):
        """Get a pandas Series from the clinical dataframe, with sample ids as the index, and each sample's status (tumor or normal) as the values."""
        clinical = self.get_clinical()
        status_map = clinical["Sample_Tumor_Normal"] 
        status_map.name = "Sample_Status"
        return status_map

    def _is_valid_df(self, df_name, df_type):
        """Tells you whether a dataframe with this name is valid for use as an omics or metadata dataframe in one of the utilties functions. Also prints message informing user.

        Parameters:
        df_name (str): The dataframe name to check.
        df_type (str): Which type of dataframe we're validating--either "omics" or "metadata"

        Returns:
        bool: Indicates whether the dataframe of that name would be valid for use as an omics or metadata dataframe in a utilities function.
        """
        if not isinstance(df_name, str): # Check that they passed a str, since utilities functions used to directly accept dataframes
            print(f"Please pass a str for {df_type} dataframe name parameter. You passed a {type(df_name)}")
            return False

        if df_type == "omics":
            valid_dfs = self._valid_omics_dfs
        elif df_type == "metadata":
            valid_dfs = self._valid_metadata_dfs
        else:
            return

        if (df_name not in valid_dfs) or (df_name not in self._data.keys()):
            print(f"{df_name} is not a valid {df_type} dataframe for this function in this dataset. Valid options:")
            for valid_name in valid_dfs:
                if valid_name in self._data.keys(): # Only print it if it's included in this dataset
                    print('\t' + valid_name)
            return False
        else:
            return True

    def _get_omics_cols(self, omics_df_name, genes):
        """Based on a single gene, or a list or array-like of genes, select multiple columns from an omics dataframe, and return the selected columns as one dataframe.

        Parameters:
        omics_df_name (str): Name of omics dataframe to select column(s) from.
        genes (str, or list or array-like of str): Gene(s) to use to select columns from omics_df. str if one gene, list or array-like if multiple. Passing None will select the entire omics dataframe.

        Returns:
        pandas DataFrame: The selected columns from the dataframe.
        """
        # Get our omics df
        omics_df = self._get_dataframe(omics_df_name)
        if omics_df is None: # Dataframe not included in this dataset. _get_dataframe already printed an error message.
            return

        # Process genes parameter
        if isinstance(genes, str): # If it's a single gene, make it a list so we can treat everything the same
            genes = [genes]
        elif isinstance(genes, (list, pd.core.series.Series, pd.core.indexes.base.Index)): # If it's already a list or array-like, we're all good
            pass
        elif genes is None: # If it's the default of None, rename columns and return the entire dataframe
            return_df = omics_df.rename(columns=lambda x:'{}_{}'.format(x, omics_df_name)) # Append dataframe name to end of each column header, to preserve info when we merge dataframes
            return return_df
        else: # If it's none of those, they done messed up. Tell 'em.
            print("Genes parameter \n{}\nis of invalid type {}. Valid types: str, list or array-like of str, or NoneType.".format(genes, type(genes)))
            return

        df = pd.DataFrame(index=omics_df.index.copy()) # Create an empty dataframe, which we'll fill with the columns we select using our genes, and then return.
        for gene in genes:
            if omics_df_name == 'phosphoproteomics' or omics_df_name == 'acetylproteomics':
                col_regex = r"^{0}{1}[^{1}]*$".format(gene, self._gene_separator) # Build a regex to get all columns that match the gene. Ending with "[^{1}]*$", where {1} corresponds to  self._gene_separator, makes sure that we're matching all the way up to the last occurrence of the gene separator, so that if there's a gene name with the gene separator in it, it's not matched by another gene name that's the part of the name before the separator--e.g., if the gene separator is a dash, "ANKHD1-EIF4EBP3-S2539" is matched by ANKHD1-EIF4EBP3 but not by ANKHD1.
            else:
                col_regex = r'^{}$'.format(gene)

            selected = omics_df.filter(regex=col_regex) # Find all columns that match the gene.
            if len(selected.columns) == 0: # If none of the columns matched the gene, generate a column of NaN and print a warning message
                empty_omics_df = pd.DataFrame(index=omics_df.index.copy())
                selected = empty_omics_df.assign(**{gene:np.nan}) # Create a column with gene as the name, filled with NaN
                print('{0} did not match any columns in {1} dataframe. {0}_{1} column inserted, but filled with NaN.'.format(gene, omics_df_name))

            selected = selected.rename(columns=lambda x:'{}_{}'.format(x, omics_df_name)) # Append dataframe name to end of each column header, to preserve info when we join dataframes
            df = df.join(selected, how='left') # Append the columns to our dataframe we'll return.

        return df

    def _get_metadata_cols(self, df_name, cols):
        """Select a single column or several columns from a metadata dataframe.

        Parameters:
        df_name (str): The name of the metadata dataframe to select the column(s) from.
        cols (str, or list or array-like of str): The column(s) to select from the dataframe. str if single, list or array-like of str if multiple. Passing None will select the entire dataframe.

        Returns:
        pandas DataFrame: The specified columns from the given dataframe.
        """
        # Get our dataframe
        df = self._get_dataframe(df_name)
        if df is None: # Dataframe not included in this dataset. _get_dataframe already printed error message.
            return

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
                print('{} column not found in the {} dataframe. Please double check that it is included in the dataframe.'.format(col, df_name))
                return
            selected = df.loc[:, [col]] # Select the column from the dataframe, keeping it as a dataframe
            return_df = return_df.join(selected, how='left') # Append the columns to our dataframe we'll return.

        return return_df

    def _get_genes_mutations(self, genes, mutations_filter):
        """Gets all the mutations for one or multiple genes, for all patients.

        Parameters:
        genes (str, or list or array-like of str): The gene(s) to grab mutations for. str if one, list or array-like of str if multiple.
        mutations_filter (list, optional): List of mutations to prioritize when filtering out multiple mutations, in order of priority. If none of the multiple mutations in a sample are included in filter_prefer, the function will automatically prioritize truncation over missense mutations, and then mutations earlier in the sequence over later mutations. Passing an empty list will cause this default hierarchy to be applied to all samples. Passing None will cause no filtering to be done, and all mutation data will be included, in a list.

        Returns:
        pandas DataFrame: The mutations in each patient for the specified gene(s).
        """
        somatic_mutation = self.get_mutations()

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

        # Check that they didn't make any typos in specifying filter values
        invalid_filter = False
        if mutations_filter is not None:
            for filter_val in mutations_filter:
                if (filter_val not in somatic_mutation[mutation_col].values) and (filter_val not in somatic_mutation[location_col].values):
                    print(f"Error: Filter value {filter_val} does not exist in the mutations dataframe for this dataset. Check for typos and existence.")
                    invalid_filter = True
        if invalid_filter:
            return

        # Create an empty dataframe, which we'll fill with the columns we select using our genes, and then return.
        df = pd.DataFrame(index=somatic_mutation.index.copy().drop_duplicates())
        for gene in genes:
            gene_mutations = somatic_mutation[somatic_mutation[gene_col] == gene] # Get all the mutations for that gene
            if len(gene_mutations) == 0: # If the gene doesn't match any genes in the dataframe, tell them, and return None.
                print("{} gene not found in somatic_mutation data.".format(gene))
                return
            gene_mutations = gene_mutations.drop(columns=[gene_col]) # Gene column is same for every sample, so we don't need it anymore.

            # Check whether all filter values exist for this particular gene. If not, that's fine, we just want to warn the user.
            if mutations_filter is not None:
                for filter_val in mutations_filter:
                    if (filter_val not in gene_mutations[mutation_col].values) and (filter_val not in gene_mutations[location_col].values):
                        print(f"Warning: Filter value {filter_val} does not exist in the mutations data for the {gene} gene, though it exists for other genes.")

            # Create another empty dataframe, which we'll fill with the mutation and location data for this gene, as lists
            prep_index = gene_mutations.index.copy().drop_duplicates()
            prep_columns = gene_mutations.columns.copy()
            mutation_status_idx = pd.Index([mutation_status_col]) # Prep mutation_status_col to be joined
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

                if mutations_filter is not None: # Filter multiple mutations down to just one
                    chosen_mutation, chosen_location = self._filter_multiple_mutations(mutations_filter, sample_mutations_list, sample_locations_list)
                    mutation_lists.at[sample, mutation_col] = chosen_mutation
                    mutation_lists.at[sample, location_col] = chosen_location
                else: # Include all the mutations!
                    mutation_lists.at[sample, mutation_col] = sample_mutations_list
                    mutation_lists.at[sample, location_col] = sample_locations_list

                # Also add the mutations status column
                mutation_lists.at[sample, mutation_status_col] = sample_mutation_status

            mutation_lists = mutation_lists.rename(columns=lambda x:'{}_{}'.format(gene, x)) # Add the gene name to end beginning of each column header, to preserve info when we join dataframes.
            df = df.join(mutation_lists, how='left') # Append the columns to our dataframe we'll return.

        return df

    def _join_other_to_mutations(self, other, mutations, mutations_were_filtered, show_location):
        """Join selected mutations data to selected other omics or metadata, add a Sample_Status column, fill in NaNs with Wildtype_Normal or Wildtype_Tumor, and name the dataframe.

        Parameters:
        other (pandas DataFrame): The selected data from the other type of dataframe (omics or metadata) to join with the selected mutations.
        mutations (pandas DataFrame): The selected mutations data to join with.
        mutations_were_filtered (bool): Whether multiple mutations in the mutations data were filtered down to just one, or not. Determines whether fill values are in lists or not.
        show_location (bool): Whether to include the Location column from the mutation dataframe.

        Returns:
        pandas DataFrame: The joined dataframe, with a Sample_Status column added and NaNs filled.
        """
        joined = other.join(mutations, how="left") # Left join other data and mutation data (left being the other data)

        # Add Sample_Status column by joining the sample_status_map to the joined mutation dataframe. Do a left join so we drop any indices not in the mutations dataframe.
        sample_status_map = self._get_sample_status_map()
        joined = joined.join(sample_status_map, how="left")

        # If there's no mutations filter, then based on the dtypes in the dataframe, set our fill values so that .loc will insert the value as a single item in a list, instead of unpacking the list.
        if mutations_were_filtered:
            wildtype_normal_fill = "Wildtype_Normal"
            wildtype_tumor_fill = "Wildtype_Tumor"
            no_mutation_fill = "No_mutation"
        elif (joined.dtypes == "object").all(): # If all columns in "joined" have a dtype of "object"
            wildtype_normal_fill = [["Wildtype_Normal"]]
            wildtype_tumor_fill = [["Wildtype_Tumor"]]
            no_mutation_fill = [["No_mutation"]]
        else:
            wildtype_normal_fill = [[["Wildtype_Normal"]]]
            wildtype_tumor_fill = [[["Wildtype_Tumor"]]]
            no_mutation_fill = [[["No_mutation"]]]

        # Fill in Wildtype_Normal or Wildtype_Tumor for NaN values (i.e., no mutation data for that sample) in joined dataframe mutation columns
        mutation_regex = r'^.*_Mutation$' # Construct regex to find all mutation columns
        mutation_cols = [col for col in joined.columns.values if re.match(mutation_regex, col)] # Get a list of all mutation columns
        for mutation_col in mutation_cols:
            joined.loc[(joined['Sample_Status'] == "Normal") & (pd.isnull(joined[mutation_col])), mutation_col] = wildtype_normal_fill # Change all NaN mutation values for Normal samples to Wildtype_Normal. 
            joined.loc[(joined['Sample_Status'] == "Tumor") & (pd.isnull(joined[mutation_col])), mutation_col] = wildtype_tumor_fill # Change all NaN mutation values for Tumor samples to Wildtype_Tumor

        # Depending on show_location, either fill NaN values in the joined dataframe location columns with "No_mutation", or just drop the location columns altogether
        location_regex = r'^.*_Location$' # Construct regex to find all location columns
        location_cols = [col for col in joined.columns.values if re.match(location_regex, col)] # Get a list of all location columns
        for location_col in location_cols:
            if show_location: # If we're including the location column, fill NaN with "No_mutation", since that's what it means, so things are clearer to the user.
                joined.loc[(pd.isnull(joined[location_col])) & (pd.notnull(joined['Sample_Status'])), location_col] = no_mutation_fill # Make sure Sample Status is not NaN, though--if it is, we have no mutation data at all for that sample, so we can't say "No_mutation". It must have been a sample that was in the other dataframe, but not the mutations.
            else:
                joined = joined.drop(columns=[location_col]) # Drop the location column, if the caller wanted us to.

        # Fill NaN values in Mutation_Status column with either Wildtype_Tumor or Wildtype_Normal
        mutation_status_regex = r"^.*_Mutation_Status$" # Construct a regex to find all Mutation_Status columns
        mutation_status_cols = [col for col in joined.columns.values if re.match(mutation_status_regex, col)] # Get a list of all Mutation_Status columns
        for mutation_status_col in mutation_status_cols:
            joined.loc[(joined['Sample_Status'] == "Normal") & (pd.isnull(joined[mutation_status_col])), mutation_status_col] = "Wildtype_Normal" # Change all NaN mutation status values for Normal samples to Wildtype_Normal
            joined.loc[(joined['Sample_Status'] == "Tumor") & (pd.isnull(joined[mutation_status_col])), mutation_status_col] = "Wildtype_Tumor" # Change all NaN mutation status values for Tumor samples to Wildtype_Tumor

        return joined

    def _filter_multiple_mutations(self, mutations_filter, sample_mutations_list, sample_locations_list):
        """Based on a mutations filter, choose one mutation and its location from two lists of mutations and locations.

        Parameters:
        mutations_filter (list of str): A list of mutations to prioritize, in order of priority. Passing an empty list will cause truncations to be chosen over missense, and mutations earlier in the sequence over later ones.
        sample_mutations_list (list of str): The mutations to filter.
        sample_locations_list (list of str): The locations to filter, in the same order as the mutations.

        Returns:
        str: The chosen mutation
        str: The chosen location
        """
        # Based on the cancer type, define which mutation types are truncations, for sorting later
        if self._cancer_type == 'colon':
            truncations = ['frameshift deletion', 'frameshift insertion', 'frameshift substitution', 'stopgain', 'stoploss']
            missenses = ['nonframeshift deletion', 'nonframeshift insertion', 'nonframeshift substitution', 'nonsynonymous SNV']
        else:
            truncations = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site']
            missenses = ['In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation']

        # Filter the mutations!!
        chosen_indices = []
        for filter_val in mutations_filter: # This will start at the beginning of the filter list, thus filters earlier in the list are prioritized, like we want
            if filter_val in sample_mutations_list:
                chosen_indices = [index for index, value in enumerate(sample_mutations_list) if value == filter_val]
            elif filter_val in sample_locations_list: # The mutations filter can also filter by location
                chosen_indices = [index for index, value in enumerate(sample_locations_list) if value == filter_val]
            if len(chosen_indices) > 0: # We found at least one mutation from the filter to prioritize, so we don't need to worry about later values in the filter priority list
                break

        if len(chosen_indices) == 0: # None of the mutations for the sample were in the filter, so we're going to have to use our default hierarchy
            for mutation in sample_mutations_list:
                if mutation in truncations:
                    chosen_indices += [index for index, value in enumerate(sample_mutations_list) if value == mutation]

        if len(chosen_indices) == 0: # None of them were in the filter, nor were truncations, so we'll grab all the missenses
            for mutation in sample_mutations_list:
                if mutation in missenses:
                    chosen_indices += [index for index, value in enumerate(sample_mutations_list) if value == mutation]

        if len(chosen_indices) == 0: # There were no truncations or missenses, so they should all be Silent mutations
            for mutation in sample_mutations_list:
                if mutation != "Silent":
                    print(f"Warning: Unknown mutation type {mutation}. Assigned lowest priority in filtering.")
            chosen_indices = range(len(sample_mutations_list)) # We'll sort them all by location

        # If there are multiple mutations in chosen_indices, the following code will pick the one soonest in the peptide sequence.
        soonest_mutation = sample_mutations_list[chosen_indices[0]]
        soonest_location = sample_locations_list[chosen_indices[0]]
        for index in chosen_indices:
            mutation = sample_mutations_list[index]
            location = sample_locations_list[index]                            

            if pd.isnull(location): # Some of the mutations have no location. We'll de-prioritize those.
                continue
            if pd.isnull(soonest_location): # This would happen if our initial value for soonest_location was NaN. If we got here, then the one we're testing isn't null, and we'll automatically prefer it
                soonest_location = location
                soonest_mutation = mutation
                continue

            num_location = self._parse_mutation_location(location) # Here, we're parsing the numerical position out of the location strings, for comparisons
            num_soonest_location = self._parse_mutation_location(soonest_location)
            if num_location < num_soonest_location:
                soonest_location = location
                soonest_mutation = mutation
        return soonest_mutation, soonest_location

    def _parse_mutation_location(self, location):
        """Parse the number out of the location for a mutation.

        Parameters:
        location (str): The location to parse.

        Returns:
        int: The numerical part of the location.
        """
        if pd.isnull(location):
            return location
        num = ""
        found_digits = False
        for char in location:
            if char.isdigit():
                num = num + char
                found_digits = True
            else:
                if found_digits: # We only want the first block of numbers
                    return int(num)
        return int(num) # We get here if the location ended with a digit
