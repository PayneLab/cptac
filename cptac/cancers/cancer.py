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

import logging
import pandas as pd
import numpy as np
import re
import warnings

from functools import reduce

import cptac.utils as ut

from cptac.exceptions import *
from cptac.tools.dataframe_tools import add_index_levels, join_col_to_dataframe, sort_df_by_sample_status
from cptac.tools.dataframe_tools import unionize_indices, generate_sample_status_col

class Cancer:
    """Note that all cancer datasets are class objects that inherit from cptac.dataset. Therefore
    the same function calls exist for cptac.Brca, cptac.Gbm, etc.
    """

    def __init__(self, cancer_type, attempt_update_index=True, skip_init=False):
        """Initialize variables for a Cancer object.

        Parameters:
        cancer_type (str): The cancer type requested for this dataset

        """

        self._cancer_type = cancer_type
        self._sources = {} # Child class __init__ needs to fill this

        self._valid_omics_dfs = [
            'acetylproteomics',
            'circular_RNA',
            'CNV',
            'lincRNA',
            'lipidomics',
            'metabolomics',
            'miRNA',
            'phosphoproteomics',
            'phosphoproteomics_gene',
            'proteomics',
            'somatic_mutation_binary',
            'transcriptomics', 
            'CNV_log2ratio',
            'CNV_gistic'
             ]

        self._valid_metadata_dfs = [
            "clinical",
            "derived_molecular",
            "experimental_design",
            #"followup", # Right now there are duplicate rows, so don't include follow up tables for joins.
            ] # We don't allow the treatment df, as in Ovarian, or medical_history df, as in Ccrcc, because they both have multiple rows for each sample.


        #ignore logging messages
        logger = logging.getLogger()
        logger.setLevel(logging.CRITICAL)

    def set_source_version(self, source, version):
        """Set the data version you wish to use for a single source

        Parameters:
        source (string): the desired source (e.g. 'broad' or 'washu')
        version (string): the desired version (e.g. '3.1' or 'latest')
        """
        self._sources[source].set_version(version)
        print(f"{self._cancer_type} cancer data from source {source} now using data freeze {self._sources[source].get_version()}")

    def delete_df(self, df_type, source='all'):
        '''This function enables users to delete dataframes they no longer need to free up RAM

        Parameters:
        df_type (string): The type of dataframe to delete. For example, if the proteomics dataframe was loaded 
            with get_proteomics() then they would pass 'proteomics' to this function to delete it.
        source (string): The source from which to delete the dataframe of type df_type
            Default: 'all' loops through all sources, checks for that df_type and deletes it if present
        '''
        # list of sources from which the df_type was deleted
        deleted_from = None

        if source == 'all':
            deleted_from = list()
            for key in self._sources.keys():
                if df_type in self._sources[key]._data:
                    del self._sources[key]._data[df_type]
                    deleted_from.append(key)
        else:
            if df_type in self._sources[source]._data:
                del self._sources[source]._data[df_type]
                deleted_from = source

        if deleted_from:
            print(f"{df_type} deleted from {deleted_from}.")
        else:
            warnings.warn(f"{df_type} not found for deletion. Perhaps you misspelled the df_type ({df_type}) or meant to delete a different dataframe?")


    # Clinical table getters
    def get_clinical(self, source=None, tissue_type="both", imputed=False):
        """Get the clinical dataframe from the specified data source."""
        return self.get_dataframe("clinical", source, tissue_type, imputed=imputed)

    def get_demographic(self, source=None, tissue_type="both", imputed=False):
        """Get the demographic dataframe from the specified data source."""
        return self.get_dataframe("demographic", source, tissue_type, imputed=imputed)

    def get_derived_molecular(self, source=None, tissue_type="both"):
        """Get the derived_molecular dataframe."""
        return self.get_dataframe("derived_molecular", source, tissue_type)

    def get_experimental_design(self, source=None, tissue_type="both"):
        """Get the experimental_design dataframe."""
        return self.get_dataframe("experimental_design", source, tissue_type)

    def get_medical_conditions(self, source=None, tissue_type="both", imputed=False):
        """Get the medical_conditions dataframe from the specified data source."""
        return self.get_dataframe("medical_conditions", source, tissue_type, imputed=imputed)

    def get_medical_history(self, source=None, tissue_type="both"):
        """Get the medical_history dataframe."""
        return self.get_dataframe("medical_history", source, tissue_type)

    def get_previous_cancer(self, source=None, tissue_type="both", imputed=False):
        """Get the previous_cancer dataframe from the specified data source."""
        return self.get_dataframe("previous_cancer", source, tissue_type, imputed=imputed)

    def get_cancer_diagnosis(self, source=None, tissue_type="both", imputed=False):
        """Get the cancer_diagnosis dataframe from the specified data source."""
        return self.get_dataframe("cancer_diagnosis", source, tissue_type, imputed=imputed)

    def get_followup(self, source=None, tissue_type="both", imputed=False):
        """Get the followup dataframe from the specified data source."""
        return self.get_dataframe("followup", source, tissue_type, imputed=imputed)



    # Quantitative table getters
    def get_acetylproteomics(self, source, tissue_type="both", imputed=False):
        """Get the acetylproteomics dataframe from the specified data source."""
        return self.get_dataframe("acetylproteomics", source, tissue_type, imputed=imputed)

    def get_acetylproteomics_gene(self, source=None, tissue_type="both"):
        """Get the acetylproteomics_gene dataframe. The gene level acetylation measurement is an aggregate metric which potentially averages together individual measurements of different sites. Use get_acetylproteomics() to view the data for individual sites."""
        return self.get_dataframe("acetylproteomics_gene", source, tissue_type)

    def get_circular_RNA(self, source=None, tissue_type="both", imputed=False):
        """Get a circular RNA dataframe from the specified data source."""
        return self.get_dataframe("circular_RNA", source, tissue_type, imputed=imputed)

    def get_CNV(self, source=None, tissue_type="both", imputed=False):
        """Get a CNV dataframe from the specified data source."""
        return self.get_dataframe("CNV", source, tissue_type, imputed=imputed)

    def get_deconvolution_cibersort(self, source=None, tissue_type="both", imputed=False):
        """Get a deconvolution dataframe from the specified data source using the cibersort algorithm."""
        return self.get_dataframe('xcell', source, tissue_type, imputed=imputed)

    def get_deconvolution_xcell(self, source=None, tissue_type="both", imputed=False):
        """Get a deconvolution dataframe from the specified data source using the cibersort algorithm."""
        return self.get_dataframe('cibersort', source, tissue_type, imputed=imputed)

    def get_deconvolution(self, deconv_algorithm=None, source=None, tissue_type="both", imputed=False):
        """Get a deconvolution dataframe from the specified data source.

        Parameters:
        deconv_algorithm (str):  Choose an alorithm. Acceptable values are ['cibersort', 'xcell']. 
        source (str): Select data generated by a certain institution. Available sources are ['washu']. Defaults to 'washu'. 
        """
        valid_algs = ['cibersort', 'xcell']
        if deconv_algorithm is None or deconv_algorithm not in valid_algs:
            raise InvalidParameterError(f"Please pass a valid value to the 'deconv_algorithm' parameter to specify which algorithm you want deconvolution data from. Valid options are {valid_algs}.")

        return self.get_dataframe(deconv_algorithm, source, tissue_type, imputed=imputed)

    def get_lincRNA(self, source=None, tissue_type="both"):
        """Get the lincRNA dataframe."""
        return self.get_dataframe("lincRNA",source, tissue_type)

    def get_lipidomics(self, source=None, tissue_type="both"):
        """Get the lipidomics dataframe."""
        return self.get_dataframe("lipidomics", source, tissue_type)

    def get_gene_fusion(self, source=None):
        """Get the gene_fusion dataframe."""
        return self.get_dataframe("gene_fusion", source)

    def get_metabolomics(self, source=None, tissue_type="both"):
        """Get the metabolomics dataframe."""
        return self.get_dataframe("metabolomics", source, tissue_type)

    def get_methylation(self, source=None, tissue_type="both"):
        """Get the methylation dataframe."""
        return self.get_dataframe("methylation", source, tissue_type)

    def get_miRNA(self, source=None, miRNA_type='total', tissue_type="both", imputed=False):
        """Get miRNA dataframe from the specified data source.

        Parameters:
        source (str): Select data generated by a certain institution. Available sources are ['washu']. Defaults to 'washu'.
        miRNA_type (str): Choose the type of miRNA molecules measured. Acceptable values are ['total', 'precursor', 'mature'].
        """
        return self.get_dataframe(miRNA_type+'_miRNA', source, tissue_type, imputed=imputed)

    def get_phosphoproteomics(self, source=None, tissue_type="both", imputed=False):
        """Get the phosphoproteomics dataframe from the specified data source."""
        return self.get_dataframe("phosphoproteomics", source, tissue_type, imputed=imputed)

    def get_phosphoproteomics_gene(self, source=None, tissue_type="both"):
        """Get the phosphoproteomics_gene dataframe. The gene level phosphorylation measurement is an aggregate metric which potentially averages together individual measurements of different sites. Use get_phosphoproteomics() to view the data for individual sites."""
        return self.get_dataframe("phosphoproteomics_gene", source, tissue_type)

    def get_phosphosites(self, genes, source=None):
        """Returns dataframe with all phosphosites of specified gene or list of genes.

        Parameters:
        genes (str, or list or array-like of str): gene or list of genes to use to select phosphosites. str if single, list or array-like of str if multiple.

        Returns:
        pandas.DataFrame: The phosphoproteomics for the specified gene(s).
        """
        return self._get_omics_cols("phosphoproteomics", genes, source)

    def get_proteomics(self, source=None, tissue_type="both", imputed=False):
        """Get the proteomics dataframe from the specified data source."""
        return self.get_dataframe("proteomics", source, tissue_type, imputed=imputed)

    def get_somatic_mutation(self, source=None, tissue_type="both", imputed=False):
        """Get the somatic mutation dataframe from the specified data source."""
        return self.get_dataframe("somatic_mutation", source, tissue_type, imputed=imputed)

    def get_somatic_mutation_binary(self, source=None):
        """Get the somatic_mutation_binary dataframe, which has a binary value indicating, for each location on each gene, whether there was a mutation in that gene at that location, for each sample."""
        return self.get_dataframe("somatic_mutation_binary", source)

    def get_targeted_phosphoproteomics(self, source=None, tissue_type="both"):
        """Get the targeted_phosphoproteomics dataframe."""
        return self.get_dataframe("targeted_phosphoproteomics", source, tissue_type)

    def get_targeted_proteomics(self, source=None, tissue_type="both"):
        """Get the targeted_proteomics dataframe."""
        return self.get_dataframe("targeted_proteomics", source, tissue_type)

    def get_transcriptomics(self, source=None, tissue_type="both", imputed=False):
        """Get the transcriptomics dataframe from the specified data source."""
        """source (str): Select data generated by a certain institution. Available sources are ['washu','bcm','broad']."""
        return self.get_dataframe("transcriptomics", source, tissue_type, imputed=imputed)

    def get_treatment(self, source=None, tissue_type="both"):
        """Get the treatment dataframe."""
        return self.get_dataframe("treatment", source, tissue_type)

    def get_tumor_purity(self, source=None, tissue_type="both", imputed=False):
        """Get the tumor purity dataframe from the specified data source."""
        return self.get_dataframe("tumor_purity", source, tissue_type, imputed=imputed)

    def get_ubiquitylomics(self, source=None, tissue_type="both", imputed=False):
        """Get the ubiquitylomics dataframe from the specified data source."""
        return self.get_dataframe("ubiquitylomics", source, tissue_type, imputed=imputed)

    def get_docs(self, data, source=None, tissue_type="both", imputed=False):
        """Get the readme docs for the specified data type and source."""
        if source in self._datasets.keys():
            obj = self._datasets[source] 
            print('README file for ' + data, '\n')
            print(obj._readme_files['readme_'+data])
            return 0
        else:
            raise DataSourceNotFoundError(f"Data source {source} not found for the {self._cancer_type} dataset.")

    def define(self, term):
        """Print the definition a term, if it is in the dataset's list of definitions.

        Parameters:
        term (str): term to be defined

        Returns: None
        """
        if len(self._definitions.keys()) == 0:
            raise NoDefinitionsError("No definitions provided for this dataset.")
        elif term in self._definitions.keys():
            print(self._definitions[term])
        else:
            raise InvalidParameterError("{} not found in definitions. Check capitalization. Alternatively, the 'search(<your term>)' function, available for import from the cptac.utils sub-module, can be used to perform a web search of the term provided.".format(term))

    def list_definitions(self):
        """Print all terms defined in the dataset's list of definitions."""
        if len(self._definitions.keys()) > 0:
            for term in sorted(self._definitions.keys(), key=str.lower):
                print(term)
        else:
            raise NoDefinitionsError("No definitions provided for this dataset.")

    # Join functions
    # Note: These are now helper functions that call multi_join on awg data by default
    def join_omics_to_omics(self, df1_name, df2_name, df1_source=None, df2_source=None, genes1=None, genes2=None, how="outer", quiet=False, tissue_type="both"):
        """Take specified column(s) from one omics dataframe, and join to specified columns(s) from another omics dataframe. Intersection (inner join) of indices is used.

        Parameters:
        df1_name (str): Name of first omics dataframe to select columns from.
        df2_name (str): Name of second omics dataframe to select columns from.
        genes1 (str, or list or array-like of str, optional): Gene(s) for column(s) to select from df1_name. str if one key, list or array-like of str if multiple. Default of None will select entire dataframe.
        genes2 (str, or list or array-like of str, optional): Gene(s) for Column(s) to select from df2_name. str if one key, list or array-like of str if multiple. Default of None will select entire dataframe.
        how (str, optional): How to perform the join, acceptable values are from ['outer', 'inner', 'left', 'right']. Defaults to 'outer'.
        quiet (bool, optional): Whether to warn when inserting NaNs. Defaults to False.
        tissue_type (str): Acceptable values in ["tumor","normal","both"]. Specifies the desired tissue type desired in the dataframe. Defaults to "both".

        Returns:
        pandas.DataFrame: The selected columns from the two omics dataframes, joined into one dataframe.
        """
        # Check to make sure that the "how" parameter is valid
        self._check_how_parameter(how)

        if df1_source is None:
            df1_source = "awg"
            warnings.warn(f"No source specified for {df1_name} data. Source awg used, pass a source to the df1_source parameter to prevent this warning", stacklevel=3)
        if df2_source is None:
            df2_source = "awg"
            warnings.warn(f"No source specified for {df1_name} data. Source awg used, pass a source to the df2_source parameter to prevent this warning", stacklevel=3)

        # Set up parameters to work with multi_join
        df1_name = f"{df1_source} {df1_name}"
        df2_name = f"{df2_source} {df2_name}"

        if genes1 is None:
            genes1 = []
        if genes2 is None:
            genes2 = []

        return self.multi_join({df1_name:genes1, df2_name:genes2}, how=how, tissue_type=tissue_type)

    def join_omics_to_mutations(self, omics_name, mutations_genes, omics_source=None, mutations_source=None, omics_genes=None, mutations_filter=None, show_location=True, how="outer", quiet=False, tissue_type="both",mutation_cols=["Mutation","Location"]):
        """Select all mutations for specified gene(s), and joins them to all or part of the given omics dataframe. Intersection (inner join) of indices is used. Each location or mutation cell contains a list, which contains the one or more location or mutation values corresponding to that sample for that gene, or a value indicating that the sample didn't have a mutation in that gene.

        Parameters:
        omics_df (str): Name of omics dataframe to join the mutation data to.
        mutations_genes (str, or list or array-like of str): The gene(s) to get mutation data for. str if one gene, list or array-like of str if multiple.
        omics_genes (str, or list or array-like of str, optional): Gene(s) to select from the omics dataframe. str if one gene, list or array-like of str if multiple. Default will select entire dataframe.
        mutations_filter (list, optional): List of mutations to prioritize when filtering out multiple mutations, in order of priority. If none of the multiple mutations in a sample are included in mutations_filter, the function will automatically prioritize truncation over missense mutations, and then mutations earlier in the sequence over later mutations. Passing an empty list will cause this default hierarchy to be applied to all samples. Default parameter of None will cause no filtering to be done, and all mutation data will be included, in a list.
        show_location (bool, optional): Whether to include the Location column from the mutation dataframe. Defaults to True.
        how (str, optional): How to perform the join, acceptable values are from ['outer', 'inner', 'left', 'right']. Defaults to 'outer'.
        quiet (bool, optional): Whether to warn when inserting NaNs. Defaults to False.
        tissue_type (str): Acceptable values in ["tumor","normal","both"]. Specifies the desired tissue type desired in the dataframe. Defaults to "both".

        Returns:
        pandas.DataFrame: The mutations for the specified gene, joined to all or part of the omics dataframe. Each location or mutation cell contains a list, which contains the one or more location or mutation values corresponding to that sample for that gene, or a value indicating that the sample didn't have a mutation in that gene.
        """

        # Check to make sure that the "how" parameter is valid
        self._check_how_parameter(how)

        # Figure out sources, default to awg
        if omics_source is None:
            omics_source = "awg"
            warnings.warn(f"No source specified for {omics_name} data. Source awg used, pass a source to the omics_source parameter to prevent this warning", stacklevel=3)
        if mutations_source is None:
            mutations_source = "awg"
            warnings.warn(f"No source specified for mutations data. Source awg used, pass a source to the mutations_source parameter to prevent this warning", stacklevel=3)

        # Set up parameters to work with multi_join
        df1_name = f"{omics_source} {omics_name}"
        df2_name = f"{mutations_source} somatic_mutation"

        if omics_genes is None:
            genes1 = []
        else:
            genes1 = omics_genes

        return self.multi_join({df1_name:genes1, df2_name:mutations_genes}, mutations_filter=mutations_filter, how=how, tissue_type=tissue_type)

    def join_metadata_to_metadata(self, df1_name, df2_name, df1_source=None, df2_source=None, cols1=None, cols2=None, how="outer", quiet=False, tissue_type="both"):
        """Take specified column(s) from one metadata dataframe, and join to specified columns(s) from another metadata dataframe. Intersection (inner join) of indices is used.

        Parameters:
        df1_name (str): Name of first metadata dataframe to select columns from.
        df2_name (str): Name of second metadata dataframe to select columns from.
        cols1 (str, or list or array-like of str, optional): Column(s) to select from df1_name. str if one key, list or array-like of str if multiple. Default of None will select entire dataframe.
        cols2 (str, or list or array-like of str, optional): Column(s) to select from df2_name. str if one key, list or array-like of str if multiple. Default of None will select entire dataframe.
        how (str, optional): How to perform the join, acceptable values are from ['outer', 'inner', 'left', 'right']. Defaults to 'outer'.
        quiet (bool, optional): Whether to warn when inserting NaNs. Defaults to False.
        tissue_type (str): Acceptable values in ["tumor","normal","both"]. Specifies the desired tissue type desired in the dataframe. Defaults to "both".

        Returns:
        pandas.DataFrame: The selected columns from the two metadata dataframes, joined into one dataframe.
        """

        # Check to make sure that the "how" parameter is valid
        self._check_how_parameter(how)

        # Figure out sources, default to awg
        if df1_source is None:
            df1_source = "awg"
            warnings.warn(f"No source specified for {df1_name} data. Source awg used, pass a source to the df1_source parameter to prevent this warning", stacklevel=3)
        if df2_source is None:
            df2_source = "awg"
            warnings.warn(f"No source specified for {df1_name} data. Source awg used, pass a source to the df2_source parameter to prevent this warning", stacklevel=3)

        # Set up parameters to work with multi_join
        df1_name = f"{df1_source} {df1_name}"
        df2_name = f"{df2_source} {df2_name}"

        if cols1 is None:
            cols1 = []
        if cols2 is None:
            cols2 = []

        return self.multi_join({df1_name:cols1, df2_name:cols2}, how=how, tissue_type=tissue_type)

    def join_metadata_to_omics(self, metadata_name, omics_name, omics_source=None, metadata_source=None, metadata_cols=None, omics_genes=None, how="outer", quiet=False, tissue_type="both"):
        """Joins columns from a metadata dataframe (clinical, derived_molecular, or experimental_design) to part or all of an omics dataframe. Intersection (inner join) of indices is used.

        Parameters:
        metadata_name (str): Name of metadata dataframe to select columns from.
        omics_name (str): Name of omics dataframe to join the metadata columns to.
        metadata_cols (str, or list or array-like of str, optional): Column(s) to select from the metadata dataframe. str if one gene, list or array-like of str if multiple. Default is None, which will select the entire metadata dataframe.
        omics_genes (str, or list or array-like of str, optional): Gene(s) to select data for from the omics dataframe. str if one gene, list or array-like of str if multiple. Default is None, which will select entire dataframe.
        how (str, optional): How to perform the join, acceptable values are from ['outer', 'inner', 'left', 'right']. Defaults to 'outer'.
        quiet (bool, optional): Whether to warn when inserting NaNs. Defaults to False.
        tissue_type (str): Acceptable values in ["tumor","normal","both"]. Specifies the desired tissue type desired in the dataframe. Defaults to "both".

        Returns:
        pandas.DataFrame: The selected metadata columns, joined with all or part of the omics dataframe.
        """

        # Check to make sure that the "how" parameter is valid
        self._check_how_parameter(how)

        # Figure out sources, default to awg
        if omics_source is None:
            omics_source = "awg"
            warnings.warn(f"No source specified for {omics_name} data. Source awg used, pass a source to the omics_source parameter to prevent this warning", stacklevel=3)
        if metadata_name is None:
            metadata_name = "awg"
            warnings.warn(f"No source specified for {metadata_name} data. Source awg used, pass a source to the mutations_source parameter to prevent this warning", stacklevel=3)

        # Set up parameters to work with multi_join
        df1_name = f"{metadata_source} {metadata_name}"
        df2_name = f"{omics_source} {omics_name}"

        if metadata_cols is None:
            genes1 = []
        else:
            genes1 = metadata_cols
        if omics_genes is None:
            genes2 = []
        else:
            genes2 = omics_genes

        return self.multi_join({df1_name:genes1, df2_name:genes2}, how=how, tissue_type=tissue_type)

    def join_metadata_to_mutations(self, metadata_name, mutations_genes, metadata_source=None, mutations_source=None,  metadata_cols=None, mutations_filter=None, show_location=True, how="outer", quiet=False, tissue_type="both"):
        """Select all mutations for specified gene(s), and joins them to all or part of the given metadata dataframe. Intersection (inner join) of indices is used. Each location or mutation cell contains a list, which contains the one or more location or mutation values corresponding to that sample for that gene, or a value indicating that the sample didn't have a mutation in that gene.

        Parameters:
        metadata_name (str): Name of metadata dataframe to join the mutation data to.
        mutations_genes (str, or list or array-like of str): The gene(s) to get mutation data for. str if one gene, list or array-like of str if multiple.
        metadata_cols (str, or list or array-like of str, optional): Gene(s) to select from the metadata dataframe. str if one gene, list or array-like of str if multiple. Default will select entire dataframe.
        mutations_filter (list, optional): List of mutations to prioritize when filtering out multiple mutations, in order of priority. If none of the multiple mutations in a sample are included in mutations_filter, the function will automatically prioritize truncation over missense mutations, and then mutations earlier in the sequence over later mutations. Passing an empty list will cause this default hierarchy to be applied to all samples. Default parameter of None will cause no filtering to be done, and all mutation data will be included, in a list.
        show_location (bool, optional): Whether to include the Location column from the mutation dataframe. Defaults to True.
        how (str, optional): How to perform the join, acceptable values are from ['outer', 'inner', 'left', 'right']. Defaults to 'outer'.
        quiet (bool, optional): Whether to warn when inserting NaNs. Defaults to False.
        tissue_type (str): Acceptable values in ["tumor","normal","both"]. Specifies the desired tissue type desired in the dataframe. Defaults to "both".

        Returns:
        pandas.DataFrame: The mutations for the specified gene, joined to all or part of the metadata dataframe. Each location or mutation cell contains a list, which contains the one or more location or mutation values corresponding to that sample for that gene, or a value indicating that the sample didn't have a mutation in that gene.
        """

        # Check to make sure that the "how" parameter is valid
        self._check_how_parameter(how)

        # Figure out sources, default to awg
        if metadata_source is None:
            metadata_source = "awg"
            warnings.warn(f"No source specified for {metadata_name} data. Source awg used, pass a source to the metadata_source parameter to prevent this warning", stacklevel=3)
        if mutations_source is None:
            mutations_source = "awg"
            warnings.warn(f"No source specified for mutations data. Source awg used, pass a source to the mutations_source parameter to prevent this warning", stacklevel=3)

        # Set up parameters to work with multi_join
        df1_name = f"{metadata_source} {metadata_name}"
        df2_name = f"{mutations_source} somatic_mutation"

        if metadata_cols is None:
            metadata_cols = []

        return self.multi_join({df1_name:metadata_cols, df2_name:mutations_genes}, how=how, tissue_type=tissue_type, mutations_filter=mutations_filter)

    def multi_join(self, join_dict, mutations_filter=None, flatten=False, levels_to_drop=[], how="outer", tissue_type="both"):    
        """Takes a dictionary where keys include a source and datatype (either in a ('source', 'datatype') tuple or as a space separated string "source datatype"), and values are columns from those dataframes. Joins all the columns into one dataframe. If the value is an empty list it will join the entire dataframe

        Parameters:
        join_dict (dict): A dictionary with the dataframe and columns to join. Keys are the names of the dataframes and the value is a list of string with the name of the columns corresponding to each dataframe. Example: {('awg', 'phosphoproteomics'):['A2M', 'AAAS'],'awg proteomics':['AAAS', 'ZZZ3'], 'awg somatic_mutation':['AHCTF1', 'ZFHX3']}.

            Valid dataframes are: acetylproteomics, CNV, phosphoproteomics, phosphoproteomics_gene, proteomics,
            somatic_mutation_binary, somatic_mutation, transcriptomics, clinical, derived_molecular and experimental_design.

            For somatic_mutation_binary it joins all columns that match a gene. Example {'awg somatic_mutation_binary' : ['A1CF', 'ZYG11B']} It returns a dataframe with all columns that contain those genes.

        mutations_filter (list, optional): List of mutations to prioritize when filtering out multiple mutations, in order of priority. If none of the multiple mutations in a sample are included in mutations_filter, the function will automatically prioritize truncation over missense mutations, and then mutations earlier in the sequence over later mutations. Passing an empty list will cause this default hierarchy to be applied to all samples. Default parameter of None will cause no filtering to be done, and all mutation data will be included, in a list.

        flatten (bool, optional): Defaults to False and will flatten the multiindexes if set to True.

        levels_to_drop (list, optional): Defaults to empty list. Takes a list of strings. Strings are levels to be dropped. If empty it will not drop any.

        tissue_type (str): Acceptable values in ["tumor","normal","both"]. Specifies the desired tissue type desired in the dataframe. Defaults to "both"

        Returns: pandas.DataFrame
        """
        column_names = []
        to_join=[]

        # a flag for if mutations data is included, for formatting later
        format_mutations = False

        for source_data_key in join_dict.keys():
            if isinstance(source_data_key, tuple):
                source, datatype = source_data_key
            else:
                # TODO: remind user to input "source datatype" for dictionary keys in error
                source, datatype = source_data_key.split()

            # Make sure all requested data exists and is valid
            if datatype not in self._sources[source].load_functions:
                raise DataFrameNotIncludedError(f"{source} {datatype} is not a valid dataframe in the {self.get_cancer_type()} dataset.")

            ## If key belongs to omics
            if datatype in self._valid_omics_dfs:
                # If key is somatic_mutation_binary it will join all columns that match a gene
                if datatype == "somatic_mutation_binary":
                    binary_data = self.get_somatic_mutation_binary(source)
                    binary_data = "\n".join(list(binary_data.columns))
                    found_genes = []
                    if len(join_dict[source_data_key]) != 0:
                        for gene in join_dict[source_data_key]:
                            found_genes += (re.findall((gene+".+\n"), binary_data))
                        found_genes = list(map(lambda x: x.strip("\n"), found_genes))#returns a list of columns that match the given gene
                        columns = self._get_omics_cols(datatype, source, found_genes, tissue_type= tissue_type)
                    else:
                        columns = self._get_omics_cols(datatype, source, None, tissue_type= tissue_type)
                else:
                    if len(join_dict[source_data_key]) != 0:# If there are values to join it will get the columns
                        columns = self._get_omics_cols(datatype, source, join_dict[source_data_key], tissue_type = tissue_type)
                    else:# Else join all the dataframe
                        columns = self._get_omics_cols(datatype, source, None, tissue_type = tissue_type)
            ## If key belongs to metadata 
            elif datatype in self._valid_metadata_dfs:
                if len(join_dict[source_data_key]) != 0:# If there are values to join it will get the columns
                    columns = self._get_metadata_cols(datatype, source, join_dict[source_data_key], tissue_type = tissue_type)
                else:# Else join all the dataframe
                    columns = self._get_metadata_cols(datatype, source, None, tissue_type = tissue_type)
            ## If key is somatic_mutation 
            elif datatype == "somatic_mutation":
                columns = self._get_genes_mutations(source, join_dict[source_data_key], mutations_filter = mutations_filter)
                # Set flag that mutations data in join_dict
                format_mutations = True

            ### Checks if there are columns with the same name and adds the name of the
            for i in columns.columns:
                if type(i) == tuple:
                    i = reduce(lambda  x, y: str(x)+str(y), i)#returns a flattened column name
                if i in column_names:
                    columns = columns.rename(columns={i: str(i)+'_'+source+'_'+datatype})
                column_names.append(i)
            ###

            to_join.append(columns)

        joined, how = reduce(self._join_dataframe, to_join, how)

        # Format any included mutations data
        if format_mutations:
            mutations_were_filtered = mutations_filter is not None
            joined = self._format_mutations_data(joined, mutations_were_filtered, how=how, tissue_type=tissue_type)

        if len(levels_to_drop) != 0:
            joined = ut.reduce_multiindex(joined, levels_to_drop=levels_to_drop)

        if flatten == True:
            joined = ut.reduce_multiindex(joined, flatten=flatten)
        return joined

    # Help functions
    def get_cancer_type(self):
        """Return the cancer type for this dataset, as a string."""
        return self._cancer_type

    def get_data_list(self):
        """Return a list of all data currently loaded into memory"""
        complete_list = {}
        for source in self._sources.keys():
            data_list = {}
            for name in sorted(self._sources[source]._data.keys(), key=str.lower):
                df = self._sources[source]._data[name]
                data_list[name] = {'rows': df.shape[0], 'columns': df.shape[1]}
            complete_list[source] = data_list
        return complete_list

    def how_to_cite(self, cancer_type='', pmid='', unpublished=False):
        """Print instructions for citing the data."""

        # current main message
        main_message = ('Please include the following statement(s) in publications using data '
         'accessed through this module:\n"Data used in this publication were generated '
         'by the Clinical Proteomic Tumor Analysis Consortium (NCI/NIH). '
         'Data were accessed through the Python module cptac, PMID: 33560848."')
        print(main_message)

        # extra message specific to cancer type
        if (cancer_type and pmid):
            cancer_specific_message = (f"Data from {cancer_type} were originally published in PMID: {pmid}")
            print(cancer_specific_message)
        elif (cancer_type and unpublished):
            print(f"Data for {cancer_type} has not been published yet.")
        else:
            # no additional message will be printed if we have not passed in parameters
            pass

    def list_data_sources(self):
        """Print which sources provide each data type.

        Parameters:
        print_list (bool, optional): Whether to print the list. Default is True. Otherwise, it's returned as a string.
        """

        # This dict will be keyed by data type, and the values will be each source that provides that data type
        data_sources = {}

        for source in sorted(self._sources.keys()):
            for df_name in sorted(self._sources[source].load_functions.keys()):
                if df_name in ["cibersort", "xcell"]:
                    df_name = f"deconvolution_{df_name}" # For clarity

                if df_name in data_sources.keys():
                    data_sources[df_name][0] += f", {source}"
                else:
                    data_sources[df_name] = [source]

        data_sources = pd.\
        DataFrame(data_sources).\
        transpose().\
        sort_index().\
        reset_index()

        data_sources.columns=["Data type", "Available sources"]

        return data_sources

    def version(self, source):
        """Return the dataset version of this instance, as a string."""
        return self._sources[source]

    def get_dataframe(self, name, source=None, tissue_type="both", imputed=False):
        """Check that a given dataframe from a given source exists, and return a copy if it does.

        Parameters:
        name (str): The datatype for which you want the dataframe.
        source (str): The source of the dataframe.
        tissue_type (str, optional): Acceptable values in ["tumor","normal","both"]. Specifies the tissue type desired in the dataframe. Defaults to "both".
        imputed (bool, optional): whether the data is imputed. Defaults to False.

        Returns: pandas.DataFrame
        """

        if imputed:
            name = name + "_imputed"

        # If no source specified, tell user what sources are available for that datatype
        if source is None:
            sources_for_data = []
            for src in self._sources.keys():
                if name in self._sources[src].load_functions.keys():
                    sources_for_data.append(src)
            if len(sources_for_data) == 0:
                # Desired datatype does not exist
                raise DataFrameNotIncludedError(f"{name} datatype not included in the {self._cancer_type} dataset. Use <cancer object>.list_data_sources() to see which data are available.")
            elif len(sources_for_data) == 1:
                # Warn the user that a default value is being used
                source = sources_for_data[0]
                warnings.warn(f"Using source {source} for {name} data as no other sources provide this data. To remove this warning, pass {source} as the source parameter.", ParameterWarning, stacklevel=3)
            else:
                # Raise error and let user know what sources are available
                raise DataSourceNotFoundError(f"No source selected. Available sources for {self._cancer_type} {name} data are: {sources_for_data}.")

        if source in self._sources.keys():
            df = self._sources[source].get_df(name)

            # Handle tissue type and filter df as specified
            if tissue_type == "normal":
                df = self._normal_only(df)
            elif tissue_type == "tumor":
                df = self._tumor_only(df)

            return df
        else:
            raise DataSourceNotFoundError(f"Data source {source} not found for the {self._cancer_type} dataset.")

    # "Private" methods
    def _check_df_valid(self, df_name, source, df_type):
        """Checks whether a dataframe with this name is valid for use as an omics or metadata dataframe in one of the utilties functions. Throws an InvalidParameterError if it isn't.

        Parameters:
        df_name (str): The dataframe name to check.
        source (str): The source of the dataframe.
        df_type (str): Which type of dataframe we're validating--either "omics" or "metadata"

        Returns: None
        """
        if not isinstance(df_name, str): # Check that they passed a str, since utilities functions used to directly accept dataframes
            raise InvalidParameterError(f"Please pass a str for dataframe name parameter. You passed {df_name}, which is a {type(df_name)}")

        if df_type == "omics":
            valid_dfs = self._valid_omics_dfs
        elif df_type == "metadata":
            valid_dfs = self._valid_metadata_dfs
        else:
            raise CptacDevError(f"Invalid df_type of {df_type} passed to cptac.Dataset._check_df_valid.")

        if df_name not in self._sources[source].load_functions:
            raise DataFrameNotIncludedError(f"{source} {df_name} dataframe not included in the {self.get_cancer_type()} dataset.")
        elif df_name not in valid_dfs:
            error_msg = f"{df_name} is not a valid {df_type} dataframe for this function in this dataset. Valid options:"
            for valid_name in valid_dfs:
                if valid_name in self._sources[source].load_functions: # Only print it if it's included in this dataset
                    error_msg = error_msg + '\n\t' + valid_name
            raise InvalidParameterError(error_msg)

    def _tumor_only(self, df):
        """For a given dataframe, keep only the tumor samples."""
        normal_endings = ('.N', '.C') # HNSCC data has cored normal samples marked .C
        tumor_df = df[~df.index.str.endswith(normal_endings)]
        return tumor_df

    def _normal_only(self, df):
        """For a given dataframe, keep only the normal samples."""
        normal_endings = ('.N', '.C') # HNSCC data has cored normal samples marked .C
        normal_df = df[df.index.str.endswith(normal_endings)]
        return normal_df

    def _get_omics_cols(self, omics_df_name, source, genes, tissue_type="both"):
        """Based on a single gene, or a list or array-like of genes, select multiple columns from an omics dataframe, and return the selected columns as one dataframe.

        Parameters:
        omics_df_name (str): Name of omics dataframe to select column(s) from.
        source (str): Source for data to select column(s) from.
        genes (str, or list or array-like of str): Gene(s) to use to select columns from omics_df. str if one gene, list or array-like if multiple. Passing None will select the entire omics dataframe.
        tissue_type (str): Acceptable values in ["tumor","normal","both"]. Specifies the tissue type desired in the dataframe. Defaults to "both".

        Returns:
        pandas.DataFrame: The selected columns from the dataframe.
        """
        # Check that they passed a valid omics df
        self._check_df_valid(omics_df_name, source, "omics")

        # Get our omics df, using get_dataframe to catch invalid requests
        omics_df = self.get_dataframe(omics_df_name, source, tissue_type).copy()

        # Process genes parameter
        if isinstance(genes, str): # If it's a single gene, make it a list so we can treat everything the same
            genes = [genes]
        elif isinstance(genes, (list, pd.Series, pd.Index)): # If it's already a list or array-like, we're all good
            pass
        elif genes is None: # If it's the default of None, rename columns and return the entire dataframe
            # Add the gene name to end beginning of each column header, to preserve info when we join dataframes.
            if isinstance(omics_df.columns, pd.MultiIndex):
                omics_df.columns = omics_df.columns.set_levels(omics_df.columns.levels[0] + '_' + source + '_' + omics_df_name, level=0)
            else:
                omics_df = omics_df.add_suffix('_' + source + '_' + omics_df_name)
            return omics_df
        else: # If it's none of those, they done messed up. Tell 'em.
            raise InvalidParameterError("Genes parameter \n{}\nis of invalid type {}. Valid types: str, list or array-like of str, or NoneType.".format(genes, type(genes)))

        genes = pd.Index(genes, name="Name")

        if isinstance(omics_df.columns, pd.MultiIndex):
            contained = genes.intersection(omics_df.columns.get_level_values("Name")).drop_duplicates() # Get the genes that actually exist in the dataframe's columns
            mi_contained = omics_df.columns[omics_df.columns.get_level_values("Name").isin(genes)]

            not_contained = genes.difference(contained).drop_duplicates() # So we can warn the user later
            arrays = [not_contained] + [[np.nan] for i in range(omics_df.columns.nlevels - 1)]
            mi_not_contained = pd.MultiIndex.from_product(arrays, names=omics_df.columns.names)

            genes = mi_contained.union(mi_not_contained) # To use for reindexing the dataframe
        else:
            contained = genes.intersection(omics_df.columns).drop_duplicates() # Get the genes that actually exist in the dataframe's columns
            not_contained = genes.difference(contained).drop_duplicates() # So we can warn the user later

        selected = omics_df[contained]
        selected = selected.reindex(columns=genes) # This will add the columns not included in the dataframe, and fill them with NaN.

        # Warn the user about columns filled with NaN
        if len(not_contained) > 0:
            warnings.warn(f"The following columns were not found in the {source} {omics_df_name} dataframe, so they were inserted into joined table, but filled with NaN: {', '.join(not_contained)}", ParameterWarning, stacklevel=3)

        # Append dataframe name to end of each column header, to preserve info when we merge dataframes
        if isinstance(omics_df.columns, pd.MultiIndex):
            selected.columns = selected.columns.set_levels(selected.columns.levels[0] + '_' + source + '_' + omics_df_name, level=0)
        else:
            selected = selected.add_suffix('_' + source + '_' + omics_df_name)

        selected.columns.name = "Name"
        return selected

    def _get_metadata_cols(self, df_name, source, cols, tissue_type="both"):
        """Select a single column or several columns from a metadata dataframe.

        Parameters:
        df_name (str): The name of the metadata dataframe to select the column(s) from.
        source (str): The source for the metadata.
        cols (str, or list or array-like of str): The column(s) to select from the dataframe. str if single, list or array-like of str if multiple. Passing None will select the entire dataframe.
        tissue_type (str): Acceptable values in ["tumor","normal","both"]. Specifies the tissue type desired in the dataframe. Defaults to "both".

        Returns:
        pandas.DataFrame: The specified columns from the given dataframe.
        """
        # Check that they passed a valid metadata df
        self._check_df_valid(df_name, source, "metadata")

        # Get our dataframe, using get_dataframe to catch invalid requests
        df = self.get_dataframe(df_name, source, tissue_type).copy()

        # Process genes parameter
        if isinstance(cols, str): # If it's a single column, make it a list so we can treat everything the same
            cols = [cols]
        elif isinstance(cols, (list, pd.Series, pd.Index)): # If it's already a list or array-like, we're all good
            pass
        elif cols is None: # If it's the default of None, return the entire dataframe
            return df
        else: # If it's none of those, they done messed up. Tell 'em.
            raise InvalidParameterError("Columns parameter {} is of invalid type {}. Valid types: str, or list or array-like of str.".format(cols, type(cols)))

        cols = pd.Index(cols).drop_duplicates()

        # Check that they didn't pass any invalid columns
        not_contained = cols.difference(df.columns)
        if len(not_contained) > 0:
            raise InvalidParameterError(f'The following columns were not found in the {source} {df_name} dataframe: {", ".join(not_contained)}')

        selected = df[cols]
        return selected

    def _get_genes_mutations(self, source, genes, mutations_filter, mutation_cols = ["Mutation","Location"]):
            """Gets all the mutations for one or multiple genes, for all patients.

            Parameters:
            genes (str, or list or array-like of str): The gene(s) to grab mutations for. str if one, list or array-like of str if multiple.
            mutations_filter (list[str] or str, optional): List of mutations to prioritize when filtering out multiple mutations, in order of priority. If none of the multiple mutations in a sample are included in mutations_filter, the function will automatically prioritize truncation over missense mutations, and then mutations earlier in the sequence over later mutations. Passing an empty list will cause this default hierarchy to be applied to all samples. Passing None will cause no filtering to be done, and all mutation data will be included, in a list.
            mutation_cols (list or str): List of columns to include in joined df. Default is the Mutation and Location column. The str "All" returns all mutation columns. 

            Returns:
            pandas.DataFrame: The mutations in each patient for the specified gene(s).
            """

            somatic_mutation = self.get_somatic_mutation(source)

            # Process genes parameter
            if isinstance(genes, str): # If it's a single gene, make it a list so we can treat everything the same
                genes = [genes]
            elif isinstance(genes, (list, pd.Series, pd.Index)): # If it's already a list or array-like, we're all good
                pass
            else: # If it's neither of those, they done messed up. Tell 'em.
                raise InvalidParameterError("Genes parameter {} is of invalid type {}. Valid types: str, or list or array-like of str.".format(genes, type(genes)))

            # Process mutations_filter parameter
            if mutations_filter is not None:
                if isinstance(mutations_filter, str): # If it's a single gene, make it a list so we can treat everything the same
                    mutations_filter = [mutations_filter]
                elif isinstance(mutations_filter, (list, pd.Series, pd.Index)): # If it's already a list or array-like, we're all good
                    pass
                else: # Again, let them know if their parameter is not going to cut it
                    raise InvalidParameterError("mutations_filter parameter {} is of invalid type {}. Valid types: str, or list or array-like of str.".format(mutations_filter, type(mutations_filter)))

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
                        raise InvalidParameterError(f"Filter value {filter_val} does not exist in the mutations dataframe for this dataset. Check for typos and existence. Merge aborted.")

            # Create an empty dataframe, which we'll fill with the columns we select using our genes, and then return.
            df = pd.DataFrame(index=somatic_mutation.index.copy().drop_duplicates())
            genes = pd.Series(genes).drop_duplicates()
            for gene in genes:
                gene_mutations = somatic_mutation[somatic_mutation[gene_col] == gene] # Get all the mutations for that gene
                if len(gene_mutations) == 0: # If the gene doesn't match any genes in the dataframe, tell them
                    raise InvalidParameterError("{} gene not found in somatic_mutation data.".format(gene))
                gene_mutations = gene_mutations.drop(columns=[gene_col]) # Gene column is same for every sample, so we don't need it anymore.

                # Check whether all filter values exist for this particular gene. If not, that's fine, we just want to warn the user.
                if mutations_filter is not None:
                    for filter_val in mutations_filter:
                        if (filter_val not in gene_mutations[mutation_col].values) and (filter_val not in gene_mutations[location_col].values):
                            warnings.warn(f"Filter value {filter_val} does not exist in the mutations data for the {gene} gene, though it exists for other genes.", ParameterWarning, stacklevel=3)

                # Create another empty dataframe, which we'll fill with the mutation and location data for this gene, as lists
                prep_index = gene_mutations.index.copy().drop_duplicates()
                prep_columns = gene_mutations.columns.copy()
                mutation_status_idx = pd.Index([mutation_status_col]) # Prep mutation_status_col to be joined
                prep_cols_with_mut_status = prep_columns.append(mutation_status_idx) # Add a mutation_status column, which will indicate if there are 1 or multiple mutations
                mutation_lists = pd.DataFrame(index=prep_index, columns=prep_cols_with_mut_status)

                if mutation_cols == "All":
                    mutation_cols = somatic_mutation.columns.to_list()
                    mutation_cols = mutation_cols[1:] #drop column name "Gene"  

                # Get the mutation(s), mutation status, and location information for this gene and sample
                # Yes, I know I'm doing that horrible thing, using nested for loops to work with dataframes. However, I tried refactoring it to use DataFrame.groupby and DataFrame.apply, and both actually made it slower. Go figure.
                for sample in mutation_lists.index: # samples ids with mutation

                    sample_data = gene_mutations.loc[sample] # Get slice of dataframe for the sample

                    for col in mutation_cols:

                        sliced_col = sample_data[col] # Get single col 

                        #Make Mutations List 
                        if col == mutation_col:
                            if isinstance(sliced_col, pd.Series):
                                sample_mutations_list = sliced_col.tolist()
                            else:
                                sample_mutations_list = [sliced_col] 

                            # Figure out what our mutation status is (either single_mutation or multiple_mutation)
                            if len(sample_mutations_list) > 1:
                                sample_mutation_status = "Multiple_mutation"
                            else:
                                sample_mutation_status = "Single_mutation"
                        #Make Location List        
                        if col == location_col:
                            if isinstance(sliced_col, pd.Series):
                                sample_locations_list = sliced_col.tolist()
                            else:
                                sample_locations_list = [sliced_col]
                        else:
                            if isinstance(sliced_col, pd.Series):
                                sample_mut_col_list = sliced_col.tolist()
                            else:
                                sample_mut_col_list = [sliced_col]

                            mutation_lists.at[sample, col] = sample_mut_col_list  

                    necessary_cols = ["Mutation","Location"] #must be included in order to filter mutations                   
                    if mutations_filter is not None: # Filter multiple mutations down to just one
                        #Check that Mutation and Location is mutation_cols so mutations_filter can work 
                        if (all(x in mutation_cols for x in necessary_cols)):
                            chosen_mutation, chosen_location = self._filter_multiple_mutations(mutations_filter, sample_mutations_list, sample_locations_list)
                            mutation_lists.at[sample, mutation_col] = chosen_mutation
                            mutation_lists.at[sample, location_col] = chosen_location

                    else: # Include all the mutations!
                        if (all(x in mutation_cols for x in necessary_cols)):
                            mutation_lists.at[sample, mutation_col] = sample_mutations_list
                            mutation_lists.at[sample, location_col] = sample_locations_list 


                    # Also add the mutations status column
                    mutation_lists.at[sample, mutation_status_col] = sample_mutation_status

                # Print warning if either Location and Mutation not provided and mutation_filter provided                       
                if (mutations_filter is not None) & ((all(x in mutation_cols for x in necessary_cols)) == False):
                    warnings.warn(f"mutations_filter was not applied because columns 'Location' and 'Mutation' were not included in mutation_cols.", ParameterWarning, stacklevel=3)                        

                mutation_lists = mutation_lists[mutation_cols + [mutation_status_col]] #only include user specified columns and Mutation_Status 
                mutation_lists = mutation_lists.add_prefix(gene + '_') # Add the gene name to end beginning of each column header, to preserve info when we join dataframes.
                df = df.join(mutation_lists, how='outer') # Append the columns to our dataframe we'll return.
                df.columns.name = "Name"

            return df

    def _format_mutations_data(self, mutations, mutations_were_filtered, show_location=True, how="outer", quiet=False, tissue_type="both", mutation_cols=["Mutation","Location"]):
        """Format mutations data. Add a Sample_Status column, fill in NaNs with Wildtype_Normal or Wildtype_Tumor.
        Note: This is mostly so that the multi_join function can behave the same way as the old join functions.
        join_other_to_mutations does this same formatting, and should probably be refactored to use this.
        Or, since none of those old functions work currently anyway, they could be removed. Or become helper functions that call multi_join

        Parameters:
        mutations (pandas.DataFrame): The selected mutations data to format.
        mutations_were_filtered (bool): Whether multiple mutations in the mutations data were filtered down to just one, or not.
            Determines whether fill values are in lists or not.
        show_location (bool): Whether to include the Location column from the mutation dataframe.
        quiet (bool): Whether to show warning when filling in rows with no mutation data with "Wildtype_Tumor" or "Wildtype_Normal"

        Returns:
        pandas.DataFrame: The joined dataframe, with a Sample_Status column added and NaNs filled.
        """

        # Add Sample_Status column
        mutations["Sample_Status"] = "Tumor"
        mutations.loc[mutations.index.str.endswith('.N'), "Sample_Status"] = "Normal"
        sample_status_map = mutations["Sample_Status"]

        if tissue_type == "normal":
            mutations = mutations.iloc[0:0] #If tissue type is normal, we drop all of the mutations rows and join only with the columns.

        # Set our fill values
        wildtype_normal_fill = "Wildtype_Normal"
        wildtype_tumor_fill = "Wildtype_Tumor"
        no_mutation_fill = "No_mutation"

        # Fill in Wildtype_Normal or Wildtype_Tumor for NaN values (i.e., no mutation data for that sample) in mutations dataframe mutation columns
        mutation_regex = r'^.*_Mutation$' # Construct regex to find all mutation columns
        mutation_cols = mutations.columns[mutations.columns.get_level_values("Name").str.match(mutation_regex)] # Get a list of all mutation columns

        fill_log = [] # We're going to keep track of value filling, and let the user know we did it.
        for mutation_col in mutation_cols:

            # See how many values we'll fill by using sum to get number of "True" in array
            num_filled = (((sample_status_map == "Normal") | (sample_status_map == "Tumor")) & (pd.isnull(mutations[mutation_col]))).sum()
            if num_filled > 0:
                if isinstance(mutation_col, tuple):
                    gene = mutation_col[0].rsplit("_", maxsplit=1)[0]
                else:
                    gene = mutation_col.rsplit("_", maxsplit=1)[0]
                # Log how many values we're going to fill for this gene
                fill_log.append(f"{num_filled} samples for the {gene} gene")

            # Impute values
            # Change all NaN mutation values for Normal samples to Wildtype_Normal.
            mutations.loc[(sample_status_map == "Normal") & (pd.isnull(mutations[mutation_col])), mutation_col] = wildtype_normal_fill
            # Change all NaN mutation values for Tumor samples to Wildtype_Tumor
            mutations.loc[(sample_status_map == "Tumor") & (pd.isnull(mutations[mutation_col])), mutation_col] = wildtype_tumor_fill

            # If we didn't filter mutations, encapsulate the fill values in lists, to match the other values in the column
            if not mutations_were_filtered:
                mutations[mutation_col] = mutations[mutation_col].apply(lambda x: x if isinstance(x, list) else [x])

        if len(fill_log) > 0 and not quiet:
            warnings.warn(f"In joining the somatic_mutation table, no mutations were found for the following samples, so they were filled with Wildtype_Tumor or Wildtype_Normal: {', '.join(fill_log)}", FilledMutationDataWarning, stacklevel=3)

        # Depending on show_location, either fill NaN values in the mutations dataframe location columns with "No_mutation", or just drop the location columns altogether
        location_regex = r'^.*_Location$' # Construct regex to find all location columns
        location_cols = mutations.columns[mutations.columns.get_level_values("Name").str.match(location_regex)] # Get a list of all location columns
        for location_col in location_cols:
            if show_location: # If we're including the location column, fill NaN with "No_mutation", since that's what it means, so things are clearer to the user.

                # Make sure Sample Status is not NaN, though--if it is, we have no mutation data at all for that sample, so we can't say "No_mutation".
                # It must have been a sample that was in the other dataframe, but not the mutations.
                mutations.loc[(pd.isnull(mutations[location_col])) & (pd.notnull(sample_status_map)), location_col] = no_mutation_fill 

                # If we didn't filter mutations, encapsulate the fill values in lists, to match the other values in the column
                if not mutations_were_filtered:
                    mutations[location_col] = mutations[location_col].apply(lambda x: x if isinstance(x, list) else [x])
            else:
                mutations = mutations.drop(columns=[location_col]) # Drop the location column, if the caller wanted us to.

        # Fill NaN values in Mutation_Status column with either Wildtype_Tumor or Wildtype_Normal
        mutation_status_regex = r"^.*_Mutation_Status$" # Construct a regex to find all Mutation_Status columns
        mutation_status_cols = mutations.columns[mutations.columns.get_level_values("Name").str.match(mutation_status_regex)] # Get a list of all Mutation_Status columns
        for mutation_status_col in mutation_status_cols:
            # Change all NaN mutation status values for Normal samples to Wildtype_Normal
            mutations.loc[(sample_status_map == "Normal") & (pd.isnull(mutations[mutation_status_col])), mutation_status_col] = wildtype_normal_fill
            # Change all NaN mutation status values for Tumor samples to Wildtype_Tumor
            mutations.loc[(sample_status_map == "Tumor") & (pd.isnull(mutations[mutation_status_col])), mutation_status_col] = wildtype_tumor_fill

        return mutations


    def _join_other_to_mutations(self, other, mutations, mutations_were_filtered, show_location, how, quiet):
        """Join selected mutations data to selected other omics or metadata, add a Sample_Status column, fill in NaNs with Wildtype_Normal or Wildtype_Tumor, and name the dataframe.

        Parameters:
        other (pandas.DataFrame): The selected data from the other type of dataframe (omics or metadata) to join with the selected mutations.
        mutations (pandas.DataFrame): The selected mutations data to join with.
        mutations_were_filtered (bool): Whether multiple mutations in the mutations data were filtered down to just one, or not. Determines whether fill values are in lists or not.
        show_location (bool): Whether to include the Location column from the mutation dataframe.
        quiet (bool): Whether to show warning when filling in rows with no mutation data with "Wildtype_Tumor" or "Wildtype_Normal"

        Returns:
        pandas.DataFrame: The joined dataframe, with a Sample_Status column added and NaNs filled.
        """
        # Make the indices the same
        if mutations.columns.nlevels != other.columns.nlevels:
            mutations.columns = add_index_levels(to=mutations.columns, source=other.columns)
        joined = other.join(mutations, how=how)

        # Add Sample_Status column by joining the sample_status_map to the joined mutation dataframe.
        sample_status_map = self._get_sample_status_map()
        joined = join_col_to_dataframe(joined, sample_status_map)
        joined.columns.name = "Name" # This attribute gets lost in the join above

        # Set our fill values
        wildtype_normal_fill = "Wildtype_Normal"
        wildtype_tumor_fill = "Wildtype_Tumor"
        no_mutation_fill = "No_mutation"

        # Fill in Wildtype_Normal or Wildtype_Tumor for NaN values (i.e., no mutation data for that sample) in joined dataframe mutation columns
        mutation_regex = r'^.*_Mutation$' # Construct regex to find all mutation columns
        mutation_cols = joined.columns[joined.columns.get_level_values("Name").str.match(mutation_regex)] # Get a list of all mutation columns

        fill_log = [] # We're going to keep track of value filling, and let the user know we did it.
        for mutation_col in mutation_cols:

            # Log how many values we're going to fill for this gene
            num_filled = (((sample_status_map == "Normal") | (sample_status_map == "Tumor")) & (pd.isnull(joined[mutation_col]))).sum() # See how many values we'll fill by using sum to get number of "True" in array
            if num_filled > 0:
                if isinstance(mutation_col, tuple):
                    gene = mutation_col[0].rsplit("_", maxsplit=1)[0]
                else:
                    gene = mutation_col.rsplit("_", maxsplit=1)[0]
                fill_log.append(f"{num_filled} samples for the {gene} gene")

            # Impute values
            joined.loc[(sample_status_map == "Normal") & (pd.isnull(joined[mutation_col])), mutation_col] = wildtype_normal_fill # Change all NaN mutation values for Normal samples to Wildtype_Normal.
            joined.loc[(sample_status_map == "Tumor") & (pd.isnull(joined[mutation_col])), mutation_col] = wildtype_tumor_fill # Change all NaN mutation values for Tumor samples to Wildtype_Tumor

            # If we didn't filter mutations, encapsulate the fill values in lists, to match the other values in the column
            if not mutations_were_filtered:
                joined[mutation_col] = joined[mutation_col].apply(lambda x: x if isinstance(x, list) else [x])

        if len(fill_log) > 0 and not quiet:
            warnings.warn(f"In joining the somatic_mutation table, no mutations were found for the following samples, so they were filled with Wildtype_Tumor or Wildtype_Normal: {', '.join(fill_log)}", FilledMutationDataWarning, stacklevel=3)

        # Depending on show_location, either fill NaN values in the joined dataframe location columns with "No_mutation", or just drop the location columns altogether
        location_regex = r'^.*_Location$' # Construct regex to find all location columns
        location_cols = joined.columns[joined.columns.get_level_values("Name").str.match(location_regex)] # Get a list of all location columns
        for location_col in location_cols:
            if show_location: # If we're including the location column, fill NaN with "No_mutation", since that's what it means, so things are clearer to the user.

                # Make sure Sample Status is not NaN, though--if it is, we have no mutation data at all for that sample, so we can't say "No_mutation". It must have been a sample that was in the other dataframe, but not the mutations.
                joined.loc[(pd.isnull(joined[location_col])) & (pd.notnull(sample_status_map)), location_col] = no_mutation_fill 

                # If we didn't filter mutations, encapsulate the fill values in lists, to match the other values in the column
                if not mutations_were_filtered:
                    joined[location_col] = joined[location_col].apply(lambda x: x if isinstance(x, list) else [x])
            else:
                joined = joined.drop(columns=[location_col]) # Drop the location column, if the caller wanted us to.

        # Fill NaN values in Mutation_Status column with either Wildtype_Tumor or Wildtype_Normal
        mutation_status_regex = r"^.*_Mutation_Status$" # Construct a regex to find all Mutation_Status columns
        mutation_status_cols = joined.columns[joined.columns.get_level_values("Name").str.match(mutation_status_regex)] # Get a list of all Mutation_Status columns
        for mutation_status_col in mutation_status_cols:
            joined.loc[(sample_status_map == "Normal") & (pd.isnull(joined[mutation_status_col])), mutation_status_col] = "Wildtype_Normal" # Change all NaN mutation status values for Normal samples to Wildtype_Normal
            joined.loc[(sample_status_map == "Tumor") & (pd.isnull(joined[mutation_status_col])), mutation_status_col] = "Wildtype_Tumor" # Change all NaN mutation status values for Tumor samples to Wildtype_Tumor

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
        elif self._cancer_type == 'hnscc' and self.version() == "0.1":
            truncations =["stopgain", "stoploss"]
            missenses = ["nonframeshift insertion", "nonframeshift deletion"]
        else:
            truncations = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site']
            missenses = ['In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation']

        if self._cancer_type == "gbm":
            noncodings = ["Intron", "RNA", "3'Flank", "Splice_Region", "5'UTR", "5'Flank", "3'UTR"]

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

        if self._cancer_type == "gbm" and len(chosen_indices) == 0: # None of them were in the filter, nor were truncations, nor missenses, so we'll grab all the noncodings
            for mutation in sample_mutations_list:
                if mutation in noncodings:
                    chosen_indices += [index for index, value in enumerate(sample_mutations_list) if value == mutation]

        if len(chosen_indices) == 0: # There were no truncations or missenses, so they should all be Silent mutations
            for mutation in sample_mutations_list:
                if mutation not in ["Silent", "synonymous SNV"]:
                    warnings.warn(f"Unknown mutation type {mutation}. Assigned lowest priority in filtering.", ParameterWarning, stacklevel=4)
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

    def _check_how_parameter(self, given_how):
        possible_values = ['outer', 'inner', 'left', 'right']
        if given_how not in possible_values:
            raise InvalidParameterError("'{}' is not a valid value for 'how'. Possible values are 'outer', 'inner', 'left', 'right'.".format(given_how))

    def _join_dataframe(self, df1, df2):
        """Joins a dataframe to another dataframe.

        Returns: pandas.DataFrame
        """
        if type(df1) == str:
            return [df2, df1]
        elif type(df1) == list:
            how = df1[1]
            if df1[0].columns.names != df2.columns.names:
                    df1[0].columns = add_index_levels(to=df1[0].columns, source=df2.columns)
                    df2.columns = add_index_levels(to=df2.columns, source=df1[0].columns)

            joined = df1[0].join(df2, how= how)
        return [joined,how]

    def _pancan_unionize_indices(self):
        '''Gets a master index of all IDs, adds IDs to clinical (clinical will have all IDs in other data), 
        adds Sample_Tumor_Normal column, and sorts clinical rows'''
        master_index = pd.Index([])

        for name, ds in self._datasets.items():
            ds_index = unionize_indices(ds._data) # unionize_indices takes 

            master_index = master_index.\
            union(ds_index).\
            drop_duplicates()

        # For some cases, there needs to be a multiindex with patient and aliquot IDs. In this case, we drop
        # tuples so they are not included in the master index
        master_index = [x for x in master_index if not isinstance(x, tuple)] # drop tuples from index (Patient_ID, Aliquot)

        new_clinical = self._datasets["mssm"]._data["clinical"] # Get the clinical dataframe
        new_clinical = new_clinical.reindex(master_index) # Take the big, comprehensive index generated by unionize_indices and use it to reindex the clinical dataframe. This will insert empty rows for any samples in the master index that weren't in the clinical dataframe before.
        sample_status_col = generate_sample_status_col(new_clinical, normal_test=lambda sample: sample[-2:] == '.N') # Parse the patient IDs and based on that generate a column that says which rows are tumor or normal. Note that if the normal IDs are marked by some format other than a ".N" at the end, you'll need to edit the "normal_test" parameter. But otherwise you're fine.
        new_clinical.insert(0, "Sample_Tumor_Normal", sample_status_col) # Add this new column into the clinical dataframe
        self._datasets["mssm"]._data['clinical'] = new_clinical # Save the edited clinical dataframe

        #TODO don't forget this sorting in whatever save_df does. We won't be unionizing indices at all though.
        #This will no longer work since we deprecated the sort_all_rows_pancan function
        #self._datasets["mssm"]._data = sort_all_rows_pancan(self._datasets["mssm"]._data) # sort clinical after adding master

    def _get_sample_status_map(self):
        """Get a pandas Series from the clinical dataframe, with sample ids as the index, and each sample's status (tumor or normal) as the values."""
        clinical = self.get_clinical()
        status_map = clinical["Sample_Tumor_Normal"]
        status_map.name = "Sample_Status"
        return status_map

    def get_genotype_all_vars(self, mutations_genes, omics_source, mutations_filter=None, show_location=True, mutation_hotspot=None):
        """Return a dataframe that has the mutation type and wheather or not it is a multiple mutation
        Parameters:
        mutation_genes (str, or list or array-like of str): The gene(s) to get mutation data for.
        mutations_filter (list, optional):  List of mutations to prioritize when filtering out multiple mutations, in order of priority.
        omics_source(str): Source of omics data ex "bcm", "washu", "broad", "umich"
        show_location (bool, optional): Whether to include the Location column from the mutation dataframe. Defaults to True.
        mutation_hotspot (optional): a list of hotspots
        """

        #If they don't give us a filter, this is the default.

        mutations_filter = ["Deletion",
                                    'Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Nonstop_Mutation', #tuncation
                                    'Missense_Mutation_hotspot',
                                    'Missense_Mutation',
                                    'Amplification',
                                    'In_Frame_Del', 'In_Frame_Ins', 'Splice_Site' ,
                                    'De_Novo_Start_Out_Frame' ,'De_Novo_Start_In_Frame', 
                                    'Start_Codon_Ins', 'Start_Codon_SNP', 
                                    'Silent',
                                    'Wildtype']

        truncations = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site']
        missenses = ['In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation']
        noncodings = ["Intron", "RNA", "3'Flank", "Splice_Region", "5'UTR", "5'Flank", "3'UTR"]



        #check that gene is in the somatic_mutation DataFrame
        somatic_mutation = self.get_somatic_mutation()
        if mutations_genes not in somatic_mutation["Gene"].unique(): #if the gene isn't in the somacic mutations df it will still have CNV data that we want
            def add_del_and_amp_no_somatic(row):
                if row[mutations_genes] <= -.2:
                    mutations = 'Deletion'

                elif row[mutations_genes] >= .2:
                    mutations = 'Amplification'
                else:
                    mutations = "No_Mutation"

                return mutations

            cnv = self.get_CNV(source = omics_source)
            #drop the database index from ccrcc and brca
            if isinstance(cnv.keys(), pd.core.indexes.multi.MultiIndex):
                drop = ['Database_ID']
                cnv = ut.reduce_multiindex(df=cnv, levels_to_drop=drop)       
            gene_cnv = cnv[[mutations_genes]]
            mutation_col = gene_cnv.apply(add_del_and_amp_no_somatic, axis=1)
            df = gene_cnv.assign(Mutation = mutation_col)
            return df


        #combine the cnv and mutations dataframe
        combined = self.join_omics_to_mutations(omics_df_name="CNV", mutations_genes=mutations_genes, omics_genes=mutations_genes, omics_source = omics_source)


        #drop the database index 
        drop = ['Database_ID']
        combined = ut.reduce_multiindex(df=combined, levels_to_drop=drop)


        #If there are hotspot mutations, append 'hotspot' to the mutation type so that it's prioritized correctly
        def mark_hotspot_locations(row):
            #iterate through each location in the current row
            mutations = []
            for location in row[mutations_genes+'_Location']:
                if location in mutation_hotspot: #if it's a hotspot mutation
                    #get the position of the location
                    position = row[mutations_genes+'_Location'].index(location)
                    #use that to change the correct mutation
                    mutations.append(row[mutations_genes+"_Mutation"][position] + "_hotspot")
                else:
                    # get the position of the location
                    position = row[mutations_genes+'_Location'].index(location)
                    mutations.append(row[mutations_genes+"_Mutation"][position])
            return mutations

        if mutation_hotspot is not None:
            combined['hotspot'] = combined.apply(mark_hotspot_locations, axis=1)
            combined[mutations_genes+"_Mutation"] = combined['hotspot']
            combined = combined.drop(columns='hotspot')


        # Based on cnv make a new column with mutation type that includes deletions and amplifications
        def add_del_and_amp(row):
            if row[mutations_genes+ "_" + omics_source + "_CNV"] <= -.2:
                mutations = row[mutations_genes+"_Mutation"] + ['Deletion']
                locations = row[mutations_genes+'_Location']+['Deletion']

            elif row[mutations_genes + "_" + omics_source+"_CNV"] >= .2:
                mutations = row[mutations_genes+"_Mutation"] + ['Amplification']
                locations = row[mutations_genes+'_Location']+['Amplification']
            else:
                mutations = row[mutations_genes+"_Mutation"]
                locations = row[mutations_genes+"_Location"]

            return mutations, locations


        combined['mutations'], combined['locations'] = zip(*combined.apply(add_del_and_amp, axis=1))

        #now that we have the deletion and amplifications, we need to prioritize the correct mutations.
        def sort(row):
            sortedcol = []
            location = []
            chosen_indices = []
            sample_mutations_list = row['mutations']
            sample_locations_list = row['locations']
            if len(sample_mutations_list) == 1: #if there's only one mutation in the list
                sortedcol.append(sample_mutations_list[0])
                location.append(sample_locations_list[0])

            else:
                for filter_val in mutations_filter: # This will start at the beginning of the filter list, thus filters earlier in the list are prioritized, like we want
                    if filter_val in sample_mutations_list:
                        chosen_indices = [index for index, value in enumerate(sample_mutations_list) if value == filter_val]
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

                if len(chosen_indices) == 0: # None of them were in the filter, nor were truncations, nor missenses, so we'll grab all the noncodings
                    for mutation in sample_mutations_list:
                        if mutation in noncodings:
                            chosen_indices += [index for index, value in enumerate(sample_mutations_list) if value == mutation]

                soonest_mutation = sample_mutations_list[chosen_indices[0]]
                soonest_location = sample_locations_list[chosen_indices[0]]
                chosen_indices.clear()
                sortedcol.append(soonest_mutation)
                location.append(soonest_location)

            return pd.Series([sortedcol, location],index=['mutations', 'locations'])

        df = combined.apply(sort, axis=1)
        combined['Mutation'] = df['mutations']
        combined['Location'] = df['locations']

        #get a sample_status column that says if the gene has multiple mutations (including dletion and amplification)
        def sample_status(row):
            if len(row['mutations']) > 1: #if there's more than one mutation
                if len(row['mutations']) == 2 and "Wildtype_Tumor" in row['mutations']: #one of the mutations might be a "wildtype tumor"
                    status ="Single_mutation"

                elif len(row['mutations']) == 2 and "Wildtype_Normal" in row['mutations']:
                    status ="Single_mutation"

                else:
                    status = "Multiple_mutation"
            else:
                if row["mutations"] == ["Wildtype_Normal"]:
                    status = "Wildtype_Normal"
                elif row['mutations'] == ['Wildtype_Tumor']:
                    status = "Wildtype_Tumor"
                else:
                    status = "Single_mutation"

            return status
        combined['Mutation_Status'] = combined.apply(sample_status, axis=1)

        #drop all the unnecessary Columns
        df = combined.drop(columns=[mutations_genes+ "_" + omics_source +"_CNV", mutations_genes+"_Mutation", mutations_genes+"_Location", mutations_genes+"_Mutation_Status", 'Sample_Status', 'mutations','locations'])
        df['Mutation'] = [','.join(map(str, l)) for l in df['Mutation']]
        df['Location'] = [','.join(map(str, l)) for l in df['Location']]
        if show_location == False: df = df.drop(columns="Location") #if they don't want us to show the location, drop it

        return df

    def _warn_inserted_nans(self, name1, name2, index1, index2):
        """Compare two indices from two dataframes, and warn the user that any rows with index values not in both indices were filled with NaNs in a join function.

        Parameters:
        name1 (str): Name of the dataframe the first index came from
        name2 (str): Name of the dataframe the second index came from
        index1 (pandas.Index): First index to compare
        index2 (pandas.Index): Second index to compare

        Returns: None
        """
        unique1 = index1.difference(index2)
        unique2 = index2.difference(index1)

        self._issue_inserted_nans_warning(unique1, name2)
        self._issue_inserted_nans_warning(unique2, name1)

    def _issue_inserted_nans_warning(self, unique, other_name):
        """Issue a warning that the samples in unique were not found in the other_name dataframe, and those column(s) were filled with NaN.

        Parameters:
        unique (list or array-like of str): The samples that weren't in the other_name dataframe.
        other_name (str): The name of the dataframe the samples weren't found in.

        Returns: None
        """
        if other_name == "somatic_mutation":
            return # This will have separate fill warnings printed, because we use different fill values.
        elif len(unique) > 0:
            warnings.warn(f"{other_name} data was not found for the following samples, so {other_name} data columns were filled with NaN for these samples: {', '.join(unique)}", InsertedNanWarning, stacklevel=4)
