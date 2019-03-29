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

    # Check that the merged dataframe has a name
    if not hasattr(merged_df, 'name'):
        print('Merged dataframe did not have a "name" attribute.')
        PASS = False

    # For each sample, check that the value in the column in the merged dataframe specified by merged_header matches the value for that sample in the column in the source dataframe specified by original_header.
    for sample in merged_df.index.values:
        original_value = original_df.loc[sample, original_header]
        merged_value = merged_df.loc[sample, merged_header]
        if (merged_value != original_value) and (pd.notna(merged_value) or pd.notna(original_value)):
            print("Merged dataframe had incorrect values.\n\tDataframe: {}\n\tSample: {}\tColumn: {}\n\tExpected: {}\tActual: {}\n".format(merged_df.name, sample, merged_header, original_value, merged_value))
            PASS = False

    return PASS

def check_merged_column_from_row(source_df, merged_df, ID_column, filter_column, filter_value, source_column, dest_column):
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
            if source_df.name.startswith('somatic'):
                if sample[-1] == 'N':
                    original_value = 'Wildtype_Normal'
                else:
                    original_value = 'Wildtype_Tumor'
            else:
                original_value == float('NaN')
        elif len(original_values) == 1:
            original_value = original_values[0]
        elif len(original_values) > 1 and source_df.name.startswith('somatic'):
            source_filtered_with_hierarchy = Utilities().add_mutation_hierarchy(source_filtered_df)
            source_filtered_with_hierarchy = source_filtered_with_hierarchy.sort_values(by = [ID_column, 'Mutation_Hierarchy'], ascending = [True,False]) #sorts by patient key, then by hierarchy so the duplicates will come with the lower number first
            original_value = source_filtered_with_hierarchy[source_column].iloc[0]
        else:
            raise ValueError('Unexpected duplicate entries in source dataframe for merged dataframe.\n\tSource dataframe: {}\n\tMerged dataframe: {}\n\tSample: {}\n\tColumn: {}\n\tValues found: {}\n'.format(source_df.name, merged_df.name, sample, source_column, original_values))

        merged_value = merged_df.loc[sample, dest_column]
        if (merged_value != original_value) and (pd.notna(merged_value) or pd.notna(original_value)):
            print("Merged dataframe had incorrect value.\n\tDataframe: {}\n\tSample: {}\tColumn: {}\n\tExpected: {}\tActual: {}\n".format(merged_df.name, sample, dest_column, original_value, merged_value))
            PASS = False

    return PASS

def test_compare_gene_single():
    """Test compare_gene for one gene (we use TP53)."""

    print("Testing compare_gene for one gene (we use TP53)...")

    # Load our dataframes
    proteomics = co.get_proteomics()
    transcriptomics = co.get_transcriptomics()

    # Set gene name
    gene = 'TP53'

    # Get our compared dataframe
    compared = co.compare_gene(proteomics, transcriptomics, gene)

    # Check that data was preserved in the two columns in the dataframe (proteomics and transcriptomics)
    PASS = True

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
