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
import numpy as mp

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
        if (type(df1) != pd.DataFrame) or (type(df2) != pd.DataFrame):
            print("Provided data not a dataframe, please check that both data inputs are dataframes")
            return

    def compare_genes(self, df1, df2, genes):
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


    def compare_clinical(self, clinical, data, clinical_col):
        """
        Parameters
        clinical: clinical dataframe for omics data to be appended with
        data: omics data for clinical data to be appended with
        clinical_col: column in clinical dataframe to be inserted into provided omics data

        Returns
        Dataframe with specified column from clinical dataframe added to specified dataframe (i.e., proteomics) for comparison and easy plotting
        """
        pass

    def compare_phosphosites(self, proteomics, phosphoproteomics, gene):
         """
        Parameters
        proteomics: the proteomics dataframe
        phosphoproteomics: the phosphoproteomics dataframe
        gene: proteomics gene to query phosphoproteomics dataframe

        Searches for any phosphosites on the gene provided

        Returns
        Dataframe with a column from proteomics for the gene specified, as well as columns for all phosphoproteomics columns beginning with the specified gene
        """
        pass

    def add_mutation_hierarchy(self, somatic): # private
        """
        Parameters
        somatic: somatic data to add mutation hierarchy to

        Retunrs
        Somatic mutation dataframe with added mutation hierarchy
        """
        pass

    def merge_somatic(self, somatic, gene, df_gene, multiple_mutations = False): # private
        """
        Parameters
        somatic: somatic mutations dataframe that will be parsed for specified gene data
        gene: string of gene to be selected for in somatic mutation data
        df_gene: selection of omics data for particular gene to be merged with somatic data
        multiple_mutations: boolean indicating whether to include multiple mutations for specified gene in an individual

        Returns
        Dataframe of merged somatic and omics dataframes based on gene provided
        """
        pass

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
        pass

    def merge_mutations_trans(self, omics, omics_gene, somatic, somatic_gene, duplicates = False):
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
        pass
