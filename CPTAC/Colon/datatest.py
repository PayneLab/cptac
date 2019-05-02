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
import CPTAC.Colon as co

def check_df_name(df, expected_name):
    """Checks that a dataframe has a "name" attribute, and that it has the proper value."""

    PASS = True

    # Check that the dataframe has a name
    has_name = True
    if not hasattr(df, 'name'):
        print('Dataframe did not have a "name" attribute.')
        has_name = False
        PASS = False

    # Check that the dataframe has the correct name
    if has_name:
        if df.name != expected_name:
            print("Dataframe had incorrect name.\n\tExpected: {}\n\tActual: {}".format(expected_name, df.name))
            PASS = False

    return PASS

def check_merged_column(original_df, merged_df, original_header, merged_header): # private
    """Checks that when a column was taken from one dataframe and added to another, no data was lost or changed.

    Parameters
    original_df: the dataframe the column was taken from
    merged_df: the merged dataframe with the column
    original_header: the column's header in the original dataframe
    merged_header: the column's header in the merged dataframe

    Returns
    bool indicating if the data in the merged column was accurate.
    """
    PASS = True

    # For each sample, check that the value in the column in the merged dataframe specified by merged_header matches the value for that sample in the column in the source dataframe specified by original_header.
    for sample in merged_df.index.values:
        original_value = original_df.loc[sample, original_header]
        merged_value = merged_df.loc[sample, merged_header]
        if (merged_value != original_value) and (pd.notna(merged_value) or pd.notna(original_value)):
            print("Merged dataframe had incorrect values.\n\tSample: {}\tColumn: {}\n\tExpected: {}\tActual: {}\n".format(sample, merged_header, original_value, merged_value))
            PASS = False

    return PASS

def check_merged_column_from_somatic(source_df, merged_df, ID_column, filter_column, filter_value, source_column, dest_column):
    """
    Parameters
    source_df: dataframe the data came from
    merged_df: dataframe the data was put in
    ID_column: string indiating the column in source_df that has the ID values that are the indices in merged_df
    filter_column: string indicating the column whose value was looked at to decide whether to take the value in the source column for a particular row, and put it in merged_df
    filter_value: the value in filter_column that indicates the data from source_column for that sample should go in merged_df
    source_column: string indicating the column in source_df from which data was taken
    dest_column: string indicating the column in merged_df where the data from source_df was put

    Returns
    bool indicating whether the data in the merged column was preserved.
    """
    PASS = True

    for sample in merged_df.index.values:
        sample_source_df = source_df.loc[source_df[ID_column] == sample] # Load a dataframe with all just the values from source_df for this sample
        source_filtered_df = sample_source_df.loc[sample_source_df[filter_column] == filter_value]
        original_values = source_filtered_df[source_column].values

        if len(original_values) == 0:
            if sample[-1] == 'N':
                original_value = 'Wildtype_Normal'
            else:
                original_value = 'Wildtype_Tumor'
        elif len(original_values) == 1:
            original_value = original_values[0]
        else:
            source_filtered_with_hierarchy = Utilities().add_mutation_hierarchy(source_filtered_df)
            source_filtered_with_hierarchy = source_filtered_with_hierarchy.sort_values(by = [ID_column, 'Mutation_Hierarchy'], ascending = [True,False]) #sorts by patient key, then by hierarchy so the duplicates will come with the lower number first
            original_value = source_filtered_with_hierarchy[source_column].iloc[0]

        merged_value = merged_df.loc[sample, dest_column]
        if (merged_value != original_value) and (pd.notna(merged_value) or pd.notna(original_value)):
            print("Merged dataframe had incorrect value.\n\tSample: {}\tColumn: {}\n\tExpected: {}\tActual: {}\n".format(sample, dest_column, original_value, merged_value))
            PASS = False

    return PASS

def test_compare_gene_single():
    """Test compare_gene for one gene (we use TP53)."""

    print("Testing compare_gene for one gene (we use TP53)...")
    PASS = True

    # Load our dataframes
    proteomics = co.get_proteomics()
    transcriptomics = co.get_transcriptomics()

    # Set gene name
    gene = 'TP53'

    # Get our compared dataframe
    compared = co.compare_gene(proteomics, transcriptomics, gene)

    # Check that our merged dataframe has a name, and that it's the correct name
    expected_name = gene
    if not check_df_name(compared, expected_name):
        PASS = False

    # Check that data was preserved in the two columns in the dataframe (proteomics and transcriptomics)
    if not check_merged_column(proteomics, compared, gene, compared.columns.values[0]):
        PASS = False
    if not check_merged_column(transcriptomics, compared, gene, compared.columns.values[1]):
        PASS = False

    # Print whether the test passed
    if PASS:
        print('PASS')
    else:
        print('FAIL')
    
def test_compare_gene_list():
    """Test compare_gene for a list of genes."""

    print("Testing compare_gene for a list of genes...")
    PASS = True

    # Load our dataframes
    proteomics = co.get_proteomics()
    transcriptomics = co.get_transcriptomics()

    # Create our gene list
    gene_list = ['TP53', 'PIK3CA', 'AURKA']

    # Get our compared dataframe
    compared = co.compare_gene(proteomics, transcriptomics, gene_list)

    # Check that our compared dataframe has a name, and that it's the correct name
    expected_name = str(len(gene_list)) + ' Genes Combined'
    if not check_df_name(compared, expected_name):
        PASS = False

    # Check that each column's data is accurate
    sorted_gene_list = sorted(gene_list) # compare_gene should sort the list

    for i in range(3):
        if not check_merged_column(proteomics, compared, sorted_gene_list[i], compared.columns.values[i]): # The first 3 columns are the proteomics data for each gene, in order
            PASS = False
        if not check_merged_column(transcriptomics, compared, sorted_gene_list[i], compared.columns.values[i + 3]): # The second 3 columns are the transcriptomics for each gene, in order
            PASS = False

    # Print whether the test passed
    if PASS: 
        print('PASS')
    else:
        print('FAIL\n')

def test_compare_clinical():
    """Test compare_clinical."""

def test_compare_phosphosites():
    """Test compare_phosphosites."""

def test_compare_mutations_cis():
    """Test compare_mutations with same gene for omics and mutation data."""

def test_compare_mutations_trans():
    """Test compare_mutations with one gene for omics, and a different gene for mutation data."""

def test_compare_mutations_full_cis():
    """Test compare_mutations_full with the same gene for omics and mutation data."""

def test_compare_mutations_full_trans():
    """Test compare_mutations_full with one gene for omics, and a different gene for mutation data."""

test_compare_gene_single()
test_compare_gene_list()
