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
import cptac

def print_test_result(PASS):
    """Prints the result of a test, based on a bool.

    Parameters:
    PASS (bool): Whether or not the test passed.
    """
    if PASS:
        print('\tPASS')
    else:
        print('\tFAIL\n')

def check_returned_is_df(returned):
    """Checks that an object is a dataframe. Prints a specific message if it's actually None, or a general message if it's something else.

    Parameters:
    returned: The object to test

    Returns:
    bool: Indicates whether the object was a dataframe.
    """
    if returned is None:
        print("Function under test returned None.")
        return False
    
    if not isinstance(returned, pd.core.frame.DataFrame):
        print("Returned object was not a dataframe. Type of object: {}".format(type(returned)))
        return False
    return True

def check_df_name(df, expected_name):
    """Checks that a dataframe has a "name" attribute, and that it has the proper value.

    Parameters:
    df (pandas.core.frame.DataFrame): The dataframe to test.
    expected_name (str): The expected name of the dataframe.

    Returns:
    bool: Whether the dataframe had a name, and it was the correct name.
    """
    # Check that the dataframe has a name
    if not hasattr(df, 'name'):
        print('Dataframe did not have a "name" attribute.')
        return False

    # Check that the dataframe has the correct name
    if df.name != expected_name:
        print("Dataframe had incorrect name.\n\tExpected: {}\n\tActual: {}".format(expected_name, df.name))
        return False
    return True

def check_df_shape(df, exp_shape):
    """Checks that a dataframe has the proper shape.

    Parameters:
    df (pandas.core.frame.DataFrame): The dataframe to test.
    exp_shape (tuple): A tuple with two elements. First element is expected number of rows, second is expected number of columns.

    Returns:
    bool: Indicates whether the dataframe had the proper shape.
    """
    act_shape = df.shape
    if exp_shape != act_shape:
        print("Dataframe dimensions did not match expected values.\n\tExpected: {}\n\tActual: {}\n".format(exp_shape, act_shape))
        return False
    return True

def check_getter(df, exp_name, exp_dim, exp_headers, coordinates, values): 
    """Test a dataframe's name, dimensions and headers, and three test values, then print whether it passed the test.

    Parameters
    df: the dataframe gotten by the getter we are testing
    exp_name: string containing the expected name of the dataframe gotten by the getter we're testing
    exp_dim: a tuple containing the expected dimensions of the dataframe, in the format (rows, columns)
    exp_headers: if the dataframe has up to 20 columns, all of the headers for the dataframe, in order. If it has more than 20 columns, then a list containing the first ten and last ten headers, in order.
    coordinates: a tuple with three elements, each element being a tuple with two elements, the first element being the int index of the row of a test value, and the second element being the int index of the column of a test value
    values: a tuple with three elements, each element being the expected value of the test value corresponding to the coordinates at the same index in the coordinates parameter 

    Returns
    bool indicating if the dataframe had the correct data.
    """
    PASS = True

    # Check that df is a dataframe, not None or something else.
    if not check_returned_is_df(df):
        return False # End test, because other tests will be useless.

    # Check name
    if not check_df_name(df, exp_name):
        PASS = False

    # Check dimensions
    if not check_df_shape(df, exp_dim):
        PASS = False

    # Check headers
    act_headers_all = list(df.columns.values)
    if len(df.columns.values) <= 20:
        act_headers = act_headers_all
    else:
        act_headers = act_headers_all[:10] + act_headers_all[-10:]

    if len(exp_headers) != len(act_headers):
        print("Unexpected number of test headers in dataframe. Expected number of headers: {}. You passed {} headers.\n".format(len(act_headers), len(exp_headers)))
        PASS = False
    else:
        for i, header in enumerate(exp_headers):
            if header != act_headers[i]:
                print("Dataframe header did not match expected value.\n\tExpected: {}\n\tActual: {}\n".format(header, act_headers[i]))
                PASS = False

    # Check test values
    act_values = [
        df.iloc[coordinates[0][0], coordinates[0][1]],
        df.iloc[coordinates[1][0], coordinates[1][1]],
        df.iloc[coordinates[2][0], coordinates[2][1]]]

    for i, value in enumerate(values):
        if act_values[i] != value:
            print("Dataframe value did not match expected value.\n\tColumn: {}\n\tIndex: {}\n\tExpected: {}\n\tActual: {}\n".format(df.columns.values[coordinates[i][1]], df.index.values[coordinates[i][0]], value, act_values[i]))
            PASS = False

    # Return whether the dataframe passed the test
    return PASS

def build_omics_regex(genes, suffix=""):
    """Builds a regex from a list of genes to grab all columns corresponding to those genes from an omics dataframe. Optional suffix allows you to search acetylproteomics and phosphoproteomics dataframes.

    Parameters:
    genes (str or list): one gene as a string, or a list of the genes, as strings
    suffix (str, optional): a suffix to put on the end of the regex, if needed. Default is empty string.

    Returns:
    str: regex to get the columns for those genes.
    """
    regex = '^('
    if isinstance(genes, str):
        genes = [genes]
    for gene in genes:
        regex = regex + gene + '|'
    regex = regex[:-1] + ')' + suffix + '$'
    return regex

def check_appended_columns(source_df, dest_df, source_headers, dest_headers=None):
    """Checks whether a column or list of columns appended to a dataframe have the same values for each index in that dataframe as they had in the dataframe they were taken from.

    Parameters:
    source_df (pandas.core.frame.DataFrame): The dataframe the columns were taken from.
    dest_df (pandas.core.frame.DataFrame): The dataframe the columns were appended to (with them appended to it).
    source_headers (str, or list or array-like of str): The header(s) of the columns to test in source_df. str if one, list or array-like of str if multiple. 
    dest_headers (str, or list or array-like of str, optional): The header(s) of the columns to test in dest_df. str if one, list or array-like of str if multiple. If provided, must be in the same order as their corresponding headers in source_headers. If not provided, header(s) in dest_df will be constructed by appending an underscore and source_df.name to each of the source_headers.
    
    Returns:
    bool: Indicates whether the specified column(s) in dest_df had the same values for each index as they did in source_df.
    """
    PASS = True
    if isinstance(source_headers, str):
        source_headers = [source_headers]

    original_dest_headers = dest_headers # Copy for later, so we can continue running conditionals based on whether dest_headers was provided
    if dest_headers is None:
        dest_headers = [header + '_' + source_df.name for header in source_headers] # Construct what the source_headers will be in the merged dataframe
    elif isinstance(dest_headers, str):
        dest_headers = [dest_headers]

    source_selection = source_df[source_headers] # Select all of the columns to test from the source dataframe
    dest_selection = dest_df[dest_headers]

    if original_dest_headers is None:
        dest_selection = dest_selection.rename(columns=lambda x: x.split('_')[0])
    else:
        dest_rename_dict = dict(zip(dest_headers, source_headers)) # Create a dict mapping a dest header to its corresponding source header. This is why if dest_headers is provided, the headers must be in the same order as their corresponding headers in source_headers.
        dest_selection = dest_selection.rename(columns=dest_rename_dict)

    if not source_selection.equals(dest_selection):
        PASS = False
        diff_bools = source_selection.eq(dest_selection)
        diff_cols = diff_bools.columns
        for col in diff_cols:
            if not diff_bools[col].all():
                source_col = source_selection[col]
                dest_col = dest_selection[col]
                col_bools = source_col.eq(dest_col)
                source_col_diffs = source_col[~col_bools]
                dest_col_diffs = dest_col[~col_bools]
                diffs_df = pd.DataFrame({"Source":source_col_diffs, "Merged":dest_col_diffs})
                print("Merged dataframe had incorrect values. Column: {}\nValues:\n{}\n\n".format(col, diffs_df))

    return PASS

def check_invalid_columns(invalid_cols, exp_num_invalid):
    """Check that we have the expected number of columns returned for an invalid key from an omics dataframe, and that they are all NaN.
    
    Parameters:
    invalid_cols(pandas DataFrame): The invalid columns found in the merged dataframe.
    exp_num_invalid (int): The expected number of invalid columns.

    Returns:
    bool: Indicates whether the test passed.
    """
    PASS = True

    # Check that we got the expected number of invalid columns
    if len(invalid_cols.columns) != exp_num_invalid:
        print("Unexpected number of columns from invalid key. Expected number: {} Actual number: {}".format(exp_num_invalid, len(invalid_cols.columns)))
        PASS = False
    
    # Check that they were all NaN
    null_bools = invalid_cols.isnull()
    if not null_bools.all().all():
        print("Unexpected values in dataframe of columns for invalid keys. Expected all NaN. Actual:\n{}".format(invalid_cols))
        PASS = False

    if PASS:
        print("(NOTE: The invalid key message above was expected.)")

    return PASS

def check_mutation_columns(mutations, merged_df, genes, show_location=True):
    """
    Parameters
    mutations (pandas.core.frame.DataFrame): The somatic_mutation dataframe.
    merged_df (pandas.core.frame.DataFrame): The merged dataframe.
    gene (str or list): The gene(s) the mutation data was collected for. str if one, list of str if multiple.
    show_location (bool, optional): Whether the location column was included in merged_df. Default is True.

    Returns
    bool: Indicates whether the mutation data for that gene and each sample in the merged dataframe matched the data in the somatic_mutation dataframe.
    """
    PASS = True

    # Set our column names for use later
    gene_col = 'Gene'
    sample_col = 'Sample_ID'
    location_col = 'Location'
    mutation_col = 'Mutation'
    sample_status_col = 'Sample_Status'
    mutation_status_col = "Mutation_Status"

    # Get a map of sample IDs to sample statuses, so we can check the Sample_Status column
    sample_status_map = en._get_sample_status_map()

    if isinstance(genes, str):
        genes = [genes]

    # Loop through all the genes
    for gene in genes:
        # Set two more column names that vary by gene
        merged_location_col = gene + '_' + location_col
        merged_mutation_col = gene + '_' + mutation_col
        merged_mutation_status_col = gene + '_' + mutation_status_col

        # Load all the mutations for the gene
        gene_df = mutations.loc[mutations[gene_col] == gene] 

        for sample in merged_df.index.values:
            # Get the rows for just this sample from our two dataframes
            sample_df = gene_df.loc[gene_df.index == sample]
            merged_sample_df = merged_df.loc[merged_df.index == sample]

            if len(sample_df.index) == 0: # There were no mutations for that gene in this sample
                original_location = [['No_mutation']]
                if sample_status_map[sample] == "Tumor":
                    original_mutation = [['Wildtype_Tumor']]
                    original_mutation_status = "Wildtype_Tumor"
                else:
                    original_mutation = [['Wildtype_Normal']]
                    original_mutation_status = "Wildtype_Normal"
                sample_dict = { # Create a prep dictionary with what the values for the different columns should be
                    sample_col:[sample], # This will be our index
                    mutation_col:original_mutation, 
                    location_col:original_location}
                sample_df = pd.DataFrame(data=sample_dict, index=sample_dict[sample_col]) # Make that dict a dataframe, and set it as our sample_df

            elif len(sample_df.index) == 1: # There was one mutation for that gene in this sample
                original_mutation = sample_df.loc[sample, mutation_col] # Get the values
                original_location = sample_df.loc[sample, location_col]

                sample_df.at[sample, mutation_col] = [original_mutation] # Put them inside of lists, to match the merged dataframe.
                sample_df.at[sample, location_col] = [original_location] 

                original_mutation_status = "Single_mutation"

            else: # There were multiple mutations for that gene in this sample
                original_mutations = sample_df[mutation_col].tolist() 
                original_locations = sample_df[location_col].tolist()

                original_mutations = [original_mutations] # In the merged dataframe, all of the multiple mutations/location values will be in one cell, in a list
                original_locations = [original_locations] # So, we nest each mutation/location list inside another list, so that the inner list isn't unpacked when we create our df.

                sample_dict = { # Create a prep dictionary with the column values
                    sample_col:[sample], # This will be our index
                    mutation_col:original_mutations,
                    location_col:original_locations}
                sample_df = pd.DataFrame(data=sample_dict, index=sample_dict[sample_col]) # Make it a dataframe, and set it as our sample_df

                original_mutation_status = "Multiple_mutation"

            # Add a Sample_Status column to our original data df
            if sample_status_map[sample] == "Tumor": # Figure out what our sample status should be
                original_sample_status = 'Tumor'
            else:
                original_sample_status = 'Normal'
            sample_df = sample_df.assign(**{sample_status_col:original_sample_status}) # Append a Sample_Status column to our original values dataframe, with what the value should be, for comparison to the merged dataframe.

            # Append a Sample_Status column to our original values dataframe, with what the value should be, for comparison to the merged dataframe.
            sample_df = sample_df.assign(**{mutation_status_col:original_mutation_status}) 

            # TODO: Add mutation status column, check, change dest to merge

            # Make a list of our columns to test
            source_cols_to_test = [mutation_col, sample_status_col, mutation_status_col]
            dest_cols_to_test = [merged_mutation_col, sample_status_col, merged_mutation_status_col]
            if show_location: # Only test the location column if it was included
                source_cols_to_test.append(location_col)
                dest_cols_to_test.append(merged_location_col)

            # Test the columns
            if not check_appended_columns(sample_df, merged_sample_df, source_cols_to_test, dest_cols_to_test):
                PASS = False

    return PASS

# Test functions that get dataframes
def test_get_clinical():
    """Test get_clinical."""

    print('Running test_get_clinical...')

    df = en.get_clinical()
    name = "clinical"
    dimensions = (144, 26)
    headers = ['Patient_ID', 'Proteomics_Tumor_Normal', 'Country', 'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity', 'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site', 'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm', 'Num_full_term_pregnancies']
    test_coord = ((79, 16), (15, 25), (88, 2))
    test_vals = (77.0, '3', 'Poland')

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_derived_molecular():
    """Test get_derived_molecular."""

    print('Running test_get_derived_molecular...')

    df = en.get_derived_molecular()
    name = 'derived_molecular'
    dimensions = (144, 125) 
    headers = ['Estrogen_Receptor', 'Estrogen_Receptor_%', 'Progesterone_Receptor', 'Progesterone_Receptor_%', 'MLH1', 'MLH2', 'MSH6', 'PMS2', 'p53', 'Other_IHC_specify', 'Log2_variant_total', 'Log2_SNP_total', 'Log2_INDEL_total', 'Genomics_subtype', 'Mutation_signature_C>A', 'Mutation_signature_C>G', 'Mutation_signature_C>T', 'Mutation_signature_T>C', 'Mutation_signature_T>A', 'Mutation_signature_T>G']
    test_coord = ((3, 4), (30, 117), (80, 52))
    test_vals = ('Intact nuclear expression', 5.459431619, -0.34)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_experimental_setup():
    """Test get_experimental_setup."""

    print('Running test_get_experimental_setup...')

    df = en.get_experimental_setup()
    name = 'experimental_setup'
    dimensions = (144, 26) 
    headers = ['Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs', 'Proteomics_Aliquot_ID', 'Proteomics_OCT', 'WXS_normal_sample_type', 'WXS_normal_filename', 'WXS_normal_UUID', 'WXS_tumor_sample_type', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID', 'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality']
    test_coord = ((2, 13), (143, 2), (67, 25))
    test_vals = ('a16b07d8-46c1-4fd9-8204-4f866aacfbec', '130N', 'PASS')

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_acetylproteomics():
    """Test get_acetylproteomics."""

    print('Running test_get_acetylproteomics...')

    df = en.get_acetylproteomics()
    name = 'acetylproteomics'
    dimensions = (144, 10862)
    headers = ['A2M-K1168', 'A2M-K1176', 'A2M-K135', 'A2M-K145', 'A2M-K516', 'A2M-K664', 'A2M-K682', 'AACS-K391', 'AAGAB-K290', 'AAK1-K201', 'ZSCAN31-K215', 'ZSCAN32-K659', 'ZW10-K634', 'ZYX-K24', 'ZYX-K25', 'ZYX-K265', 'ZYX-K272', 'ZYX-K279', 'ZYX-K533', 'ZZZ3-K117']
    test_coord = ((1, 1), (12, 10861), (90, 5849))
    test_vals = (0.47700000000000004, 0.16, 0.4098)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_proteomics():
    """Test get_proteomics."""

    print('Running test_get_proteomics...')

    df = en.get_proteomics()
    name = "proteomics"
    dimensions = (144, 10999)
    headers = ['A1BG', 'A2M', 'A2ML1', 'A4GALT', 'AAAS', 'AACS', 'AADAT', 'AAED1', 'AAGAB', 'AAK1', 'ZSWIM8', 'ZSWIM9', 'ZW10', 'ZWILCH', 'ZWINT', 'ZXDC', 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3']
    test_coord = ((34, 6003), (99, 9544), (143, 32))
    test_vals = (0.0461, 1.68, 0.904)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_transcriptomics():
    """Test get_transcriptomics."""

    print('Running test_get_transcriptomics...')

    df = en.get_transcriptomics()
    name = "transcriptomics"
    dimensions = (109, 28057)
    headers = ['A1BG', 'A1BG-AS1', 'A1CF', 'A2M', 'A2M-AS1', 'A2ML1', 'A2MP1', 'A3GALT2', 'A4GALT', 'A4GNT', 'ZWILCH', 'ZWINT', 'ZXDA', 'ZXDB', 'ZXDC', 'ZYG11A', 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3']
    test_coord = ((22, 25483), (108, 23), (101, 17748))
    test_vals = (0.82, 12.0, 6.19)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_circular_RNA():
    """Test get_circular_RNA."""

    print('Running test_get_circular_RNA...')

    df = en.get_circular_RNA()
    name = "circular_RNA"
    dimensions = (109, 4945)
    headers = ['circ_chr10_100260218_100262063_CWF19L1', 'circ_chr10_100923975_100926019_SLF2', 'circ_chr10_100923978_100926019_SLF2', 'circ_chr10_100937402_100944128_SLF2', 'circ_chr10_100937402_100950753_SLF2', 'circ_chr10_101584602_101586156_POLL', 'circ_chr10_101667886_101676436_FBXW4', 'circ_chr10_101672915_101676436_FBXW4', 'circ_chr10_101792839_101807901_OGA', 'circ_chr10_101792839_101810314_OGA', 'circ_chrX_80288906_80310233_CHMP1B2P', 'circ_chrX_80289664_80310233_CHMP1B2P', 'circ_chrX_80707427_80719656_BRWD3', 'circ_chrX_80791854_80793772_BRWD3', 'circ_chrX_84096194_84164387_RPS6KA6', 'circ_chrX_84134782_84164387_RPS6KA6', 'circ_chrX_85067127_85074391_APOOL', 'circ_chrX_85978767_85981809_CHM', 'circ_chrX_91414904_91418871_PABPC5-AS1', 'circ_chrX_9691579_9693419_TBL1X']
    test_coord = ((108, 1), (30, 4935), (73, 2003))
    test_vals = (9.08, 6.56, 0.0)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_miRNA():
    """Test get_miRNA."""

    print('Running test_get_miRNA...')

    df = en.get_miRNA()
    name = "miRNA"
    dimensions = (99, 2337)
    headers = ['hsa-let-7a-2-3p', 'hsa-let-7a-3p', 'hsa-let-7a-5p', 'hsa-let-7b-3p', 'hsa-let-7b-5p', 'hsa-let-7c-3p', 'hsa-let-7c-5p', 'hsa-let-7d-3p', 'hsa-let-7d-5p', 'hsa-let-7e-3p', 'hsa-miR-9901', 'hsa-miR-9902', 'hsa-miR-9903', 'hsa-miR-9983-3p', 'hsa-miR-9985', 'hsa-miR-9986', 'hsa-miR-99a-3p', 'hsa-miR-99a-5p', 'hsa-miR-99b-3p', 'hsa-miR-99b-5p']
    test_coord = ((5, 0), (98, 1597), (54, 2231))
    test_vals = (1.79, 1.36, 0.26)
    
    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_CNA():
    """Test get_CNA."""

    print('Running test_get_CNA...')

    df = en.get_CNA()
    name = "CNA"
    dimensions = (95, 28057)
    headers = ['A1BG', 'A1BG-AS1', 'A1CF', 'A2M', 'A2M-AS1', 'A2ML1', 'A2MP1', 'A3GALT2', 'A4GALT', 'A4GNT', 'ZWILCH', 'ZWINT', 'ZXDA', 'ZXDB', 'ZXDC', 'ZYG11A', 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3']
    test_coord = ((12, 27865), (60, 8), (94, 15439))
    test_vals = (0.11, 0.01, -0.01)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_phosphoproteomics():
    """Test get_phosphoproteomics."""

    print('Running test_get_phosphoproteomics...')

    df =  en.get_phosphoproteomics()
    name = "phosphoproteomics"
    dimensions = (144, 73212)
    headers = ['AAAS-S495', 'AAAS-S541', 'AAAS-Y485', 'AACS-S618', 'AAED1-S12', 'AAGAB-S310', 'AAGAB-S311', 'AAK1-S14', 'AAK1-S18', 'AAK1-S20', 'ZZZ3-S397', 'ZZZ3-S411', 'ZZZ3-S420', 'ZZZ3-S424', 'ZZZ3-S426', 'ZZZ3-S468', 'ZZZ3-S89', 'ZZZ3-T415', 'ZZZ3-T418', 'ZZZ3-Y399']
    test_coord = ((36, 46), (12, 72436), (96, 45361))
    test_vals = (0.579, 0.669, 0.156)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_phosphoproteomics_gene():
    """Test get_phosphoproteomics_gene."""

    print('Running test_get_phosphoproteomics_gene...')

    df = en.get_phosphoproteomics_gene()
    name = "phosphoproteomics_gene"
    dimensions = (144, 8466)
    headers = ['AAAS', 'AACS', 'AAED1', 'AAGAB', 'AAK1', 'AAMDC', 'AARS', 'AASDH', 'AATF', 'ABCA1', 'ZSCAN5C', 'ZSWIM3', 'ZSWIM8', 'ZUP1', 'ZW10', 'ZXDA', 'ZXDC', 'ZYX', 'ZZEF1', 'ZZZ3']
    test_coord =  ((2, 7999), (143, 1045), (71, 6543))
    test_vals = (-0.0879, 0.929, 0.153)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_phosphosites():
    """Test get_phosphosites."""

    print('Running test_get_phosphosites...')

    gene = 'AAK1'
    df = en.get_phosphosites(gene)
    name = 'phosphoproteomics for ' + gene
    dimensions = (144, 37)
    headers = ['AAK1-S14_phosphoproteomics', 'AAK1-S18_phosphoproteomics', 'AAK1-S20_phosphoproteomics', 'AAK1-S21_phosphoproteomics', 'AAK1-S26_phosphoproteomics', 'AAK1-S618_phosphoproteomics', 'AAK1-S623_phosphoproteomics', 'AAK1-S624_phosphoproteomics', 'AAK1-S637_phosphoproteomics', 'AAK1-S642_phosphoproteomics', 'AAK1-T448_phosphoproteomics', 'AAK1-T606_phosphoproteomics', 'AAK1-T620_phosphoproteomics', 'AAK1-T640_phosphoproteomics', 'AAK1-T653_phosphoproteomics', 'AAK1-T674_phosphoproteomics', 'AAK1-T681_phosphoproteomics', 'AAK1-T694_phosphoproteomics', 'AAK1-T848_phosphoproteomics', 'AAK1-Y687_phosphoproteomics']
    test_coord = ((5, 33), (64, 14), (128, 0))
    test_vals = (0.547, -0.5379999999999999, 0.1395)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_mutations():
    """Test get_mutations."""

    print('Running test_get_mutations...')

    df = en.get_mutations()
    name = "somatic_mutation"
    dimensions = (52560, 3)
    headers = ['Gene', 'Mutation', 'Location']
    test_coord = ((52000, 2), (12, 0), (34567, 1))
    test_vals = ('p.V167L', 'ARID1A', 'Missense_Mutation')

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_mutations_binary():
    """Test get_mutations_binary."""

    print('Running test_get_mutations_binary...')

    df = en.get_mutations_binary()
    name = "somatic_mutation_binary"
    dimensions = (95, 51559)
    headers = ['A1BG_p.E298K', 'A1BG_p.S181N', 'A1CF_p.F487L', 'A1CF_p.S236Y', 'A2ML1_p.A8V', 'A2ML1_p.G1306D', 'A2ML1_p.L1347F', 'A2ML1_p.L82I', 'A2ML1_p.P712S', 'A2ML1_p.R443Q', 'ZYG11A_p.Q442H', 'ZYG11B_p.H315R', 'ZYG11B_p.R495M', 'ZYG11B_p.R728C', 'ZYX_p.C447Y', 'ZZEF1_p.A2723V', 'ZZEF1_p.D845Y', 'ZZEF1_p.K1251E', 'ZZEF1_p.K2387Sfs*40', 'ZZZ3_p.Y891C']
    test_coord = ((94, 51558), (0, 0), (45, 25436))
    test_vals = (0, 0, 0)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

# Test merging and appending functions
def test_compare_omics_source_preservation():
    """Test that compare_omics does not alter the dataframes it pulls data from."""
    print("Running test_compare_omics_source_preservation...")
    PASS = True

    # Load the source dataframes
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics() # Acetylproteomics and phosphoproteomics have multiple columns for one gene. We use acetylproteomics to make sure compare_omics can grab all those values.

    # Set the names
    prot_name = "proteomics"
    acet_name = "acetylproteomics"

    # Copy the source dataframes so we can make sure later on that compare_omics doesn't alter them.
    prot_copy = prot.copy()
    acet_copy = acet.copy()

    # Call compare_omics on the dataframes, and make sure it doesn't return None.
    compared = en.compare_omics(prot_name, acet_name)
    if compared is None:
        print('compare_omics returned None.')
        PASS = False

    # Use the copies we made at the beginning to make sure that compare_omics didn't alter the source dataframes
    if not prot.equals(prot_copy):
        print("Proteomics dataframe was altered by compare_omics.")
        PASS = False

    if not acet.equals(acet_copy):
        print("Acetylproteomics dataframe was altered by compare_omics.")
        PASS = False

    # Indicate whether the test passed.
    print_test_result(PASS)

def test_compare_omics_default_parameters():
    """Test compare_omics with default parameters cols1=None and cols2=None."""
    print("Running test_compare_omics_default_parameters...")
    PASS = True

    # Load the source dataframes
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics() # Acetylproteomics and phosphoproteomics have multiple columns for one gene. We use acetylproteomics to make sure compare_omics can grab all those values.

    # Set the names
    prot_name = "proteomics"
    acet_name = "acetylproteomics"

    # Run the function, make sure it returned properly
    compared = en.compare_omics(prot_name, acet_name) 
    if not check_returned_is_df(compared):
        PASS = False
        print_test_result(PASS)
        return # Skip other tests, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = prot.name + ', with ' + acet.name
    if not check_df_name(compared, exp_name):
        PASS = False

    # Check dataframe shape
    exp_num_rows = len(prot.index.intersection(acet.index))
    exp_num_cols = len(prot.columns) + len(acet.columns)
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(compared, exp_shape):
        PASS = False

    # Test that all columns from proteomics were appended, and that values were preserved
    if not check_appended_columns(prot, compared, prot.columns):
        PASS = False

    # Test that all columns from acetylproteomics were appended, and that values were preserved
    if not check_appended_columns(acet, compared, acet.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_compare_omics_one_gene():
    """Test compare_omics with one genes for cols1 and cols2."""
    print("Running test_compare_omics_one_gene...")
    PASS = True

    # Load the source dataframes
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics() # Acetylproteomics and phosphoproteomics have multiple columns for one gene. We use acetylproteomics to make sure compare_omics can grab all those values.

    # Set the names
    prot_name = "proteomics"
    acet_name = "acetylproteomics"

    # Run the function, make sure it returned properly
    prot_gene = 'TP53'
    acet_gene = 'A2M'
    compared = en.compare_omics(prot_name, acet_name, prot_gene, acet_gene) 
    if not check_returned_is_df(compared):
        PASS = False
        print_test_result(PASS)
        return # Skip other tests, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{} for {}, with {} for {}".format(prot.name, prot_gene, acet.name, acet_gene)
    if not check_df_name(compared, exp_name):
        PASS = False

    # Figure out which columns from proteomics correspond to prot_gene (should be just one)
    prot_regex = build_omics_regex(prot_gene)
    prot_cols = prot.filter(regex=prot_regex)
    if len(prot_cols.columns) != 1:
        print("Unexpected number of matching proteomics columns in test.\n\tExpected: 1\n\tActual: {}".format(len(prot_cols)))
        PASS = False

    # Figure out which columns from acetylproteomics correspond to acet_gene
    acet_suffix = "-.*"
    acet_regex = build_omics_regex(acet_gene, suffix=acet_suffix)
    acet_cols = acet.filter(regex=acet_regex)

    # Check dataframe shape
    exp_num_rows = len(prot.index.intersection(acet.index))
    exp_num_cols = len(prot_cols.columns) + len(acet_cols.columns)
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(compared, exp_shape):
        PASS = False

    # Check column values
    if not check_appended_columns(prot, compared, prot_cols.columns):
        PASS = False

    if not check_appended_columns(acet, compared, acet_cols.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_compare_omics_multiple_genes():
    """Test compare_omics with lists of genes for cols1 and cols2."""
    print("Running test_compare_omics_multiple_genes...")
    PASS = True

    # Load the source dataframes
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics() # Acetylproteomics and phosphoproteomics have multiple columns for one gene. We use acetylproteomics to make sure compare_omics can grab all those values.

    # Set the names
    prot_name = "proteomics"
    acet_name = "acetylproteomics"

    # Run the function, make sure it returned properly
    prot_genes = ['A4GALT', 'TP53', 'ZSCAN30']
    acet_genes = ['AAGAB', 'AACS', 'ZW10', 'ZYX']
    compared = en.compare_omics(prot_name, acet_name, prot_genes, acet_genes) 
    if not check_returned_is_df(compared):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{} for {} genes, with {} for {} genes".format(prot.name, len(prot_genes), acet.name, len(acet_genes))
    if not check_df_name(compared, exp_name):
        PASS = False

    # Figure out which columns from proteomics correspond to prot_gene (should be the same number as number of genes in prot_genes)
    prot_regex = build_omics_regex(prot_genes)
    prot_cols = prot.filter(regex=prot_regex) # Use the regex to get all matching columns
    if len(prot_cols.columns) != len(prot_genes):
        print("Unexpected number of matching proteomics columns in test.\n\tExpected: {}\n\tActual: {}".format(len(prot_genes), len(prot_cols.columns)))
        PASS = False

    # Figure out which columns from acetylproteomics correspond to acet_gene
    acet_suffix = '-.*'
    acet_regex = build_omics_regex(acet_genes, suffix=acet_suffix)
    acet_cols = acet.filter(regex=acet_regex) # Use the regex to get all matching columns

    # Check dataframe shape
    exp_num_rows = len(prot.index.intersection(acet.index))
    exp_num_cols = len(prot_cols.columns) + len(acet_cols.columns)
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(compared, exp_shape):
        PASS = False

    # Check column values
    if not check_appended_columns(prot, compared, prot_cols.columns):
        PASS = False

    if not check_appended_columns(acet, compared, acet_cols.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_compare_omics_all_dfs():
    """Test that compare_omics works will all dataframes that are valid for the function."""
    print("Running test_compare_omics_all_dfs...")
    PASS = True

    # Load our dataframes to test, and set our genes. We call individual parameters, to make sure the columns are formatted properly.
    acet = "acetylproteomics"
    CNA = "CNA"
    phosg = "phosphoproteomics_gene"
    phoss = "phosphoproteomics"
    prot = "proteomics"
    tran = "transcriptomics"
    gene1 = 'TP53'
    gene2 = 'AAGAB'

    # Call compare_omics on the dataframes
    acet_CNA = en.compare_omics(acet, CNA, gene1, gene2)
    phosg_phoss = en.compare_omics(phosg, phoss, gene1, gene2)
    prot_tran = en.compare_omics(prot, tran, gene1, gene2)

    # Check the return values
    if not check_returned_is_df(acet_CNA):
        print("Dataframes compared: acetylproetomics and CNA.")
        PASS = False
    if not check_returned_is_df(phosg_phoss):
        print("Dataframes compared: phosphoproteomics_gene and phosphoproteomics.")
        PASS = False
    if not check_returned_is_df(prot_tran):
        print("Dataframes compared: proteomics and transcriptomics.")
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_compare_omics_invalid_dfs():
    """Test that compare_omics will not accept non-omics dataframes, and not accept omics dataframes of the wrong format."""
    print("Running test_compare_omics_invalid_dfs...")
    PASS = True

    # Load our dataframes to test with
    prot = en.get_proteomics() # We want to try a mix of valid and invalid dataframes, so we need to load this valid dataframe
    clin = en.get_clinical()
    tran_cir = en.get_circular_RNA() # Although circular_RNA is an omics dataframe, it's of the wrong format to work with compare_omics

    # Set the names
    prot_name = "proteomics"
    clin_name = "clinical"
    tran_cir_name = "circular_RNA"

    # Test with one valid dataframe and one invalid one
    comp = en.compare_omics(prot_name, tran_cir_name)
    if comp is not None:
        print("compare_omics should have returned None when passed the {} dataframe, but instead returned a {}".format(tran_cir.name, type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid dataframe message above was expected.)")

    # Test with two invalid dataframes
    comp = en.compare_omics(clin_name, tran_cir_name)
    if comp is not None:
        print("compare_omics should have returned None when passed the {} and {} dataframes, but instead returned a {}".format(clin.name, tran_cir.name, type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid dataframe message above was expected.)")

    # Print whether the test passed
    print_test_result(PASS)

def test_compare_omics_one_invalid_key():
    """Test that compare_omics will gracefully handle one invalid key."""
    print("Running test_compare_omics_one_invalid_key...")
    PASS = True

    # Load the source dataframes
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics()

    # Set the names
    prot_name = "proteomics"
    acet_name = "acetylproteomics"

    # Run the function, make sure it returned properly
    invalid = 'gobbledegook'
    acet_valid = 'AACS' 
    compared = en.compare_omics(prot_name, acet_name, invalid, acet_valid)
    if not check_returned_is_df(compared):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{} for {}, with {} for {}".format(prot.name, invalid, acet.name, acet_valid)
    if not check_df_name(compared, exp_name):
        PASS = False

    # Get our columns corresponding to the invalid key
    invalid_suffix = "_.*"
    invalid_regex = build_omics_regex(invalid, suffix=invalid_suffix)
    invalid_cols = compared.filter(regex=invalid_regex)
    exp_num_invalid = 1

    # Figure out which columns from acetylproteomics should have been grabbed for acet_valid
    acet_suffix = '-.*'
    acet_regex = build_omics_regex(acet_valid, suffix=acet_suffix)
    acet_cols = acet.filter(regex=acet_regex) # Use the regex to get all matching columns

    # Check dataframe shape
    exp_num_rows = len(prot.index.intersection(acet.index))
    exp_num_cols = exp_num_invalid + len(acet_cols.columns)
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(compared, exp_shape):
        PASS = False

    # Check invalid columns
    if not check_invalid_columns(invalid_cols, exp_num_invalid):
        PASS = False

    # Check columns for acet_gene
    if not check_appended_columns(acet, compared, acet_cols.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_compare_omics_both_invalid_keys():
    """Test that compare_omics will gracefully handle two invalid keys."""
    print("Running test_compare_omics_both_invalid_keys...")
    PASS = True

    # Load dataframes to test with, and set keys to use
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics()

    # Set the names
    prot_name = "proteomics"
    acet_name = "acetylproteomics"

    # Run the function, make sure it returned properly
    invalid = 'gobbledegook'
    compared = en.compare_omics(prot_name, acet_name, invalid, invalid)
    if not check_returned_is_df(compared):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{} for {}, with {} for {}".format(prot.name, invalid, acet.name, invalid)
    if not check_df_name(compared, exp_name):
        PASS = False

    # Get our columns corresponding to the invalid key
    invalid_suffix = "_.*"
    invalid_regex = build_omics_regex(invalid, suffix=invalid_suffix)
    invalid_cols = compared.filter(regex=invalid_regex)
    exp_num_invalid = 2

    # Check dataframe shape
    exp_num_rows = len(prot.index.intersection(acet.index))
    exp_num_cols = exp_num_invalid
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(compared, exp_shape):
        PASS = False

    # Check invalid columns
    if not check_invalid_columns(invalid_cols, exp_num_invalid):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_compare_omics_one_list_with_invalid_key():
    """Test that compare_omics will gracefully handle a list of valid keys containing one invalid key."""
    print("Running test_compare_omics_one_list_with_invalid_key...")
    PASS = True

    # Load dataframes to test with
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics()

    # Set the names
    prot_name = "proteomics"
    acet_name = "acetylproteomics"

    # Set the keys
    invalid = "gobbledegook"
    prot_valid_list = ['A4GALT', 'TP53', 'ZSCAN30']
    prot_invalid_list = prot_valid_list.copy()
    prot_invalid_list.append(invalid)
    acet_valid_list = ['AAGAB', 'AACS', 'ZW10', 'ZYX']
    
    # Run the function, make sure it returned properly
    compared = en.compare_omics(prot_name, acet_name, prot_invalid_list, acet_valid_list)
    if not check_returned_is_df(compared):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{} for {} genes, with {} for {} genes".format(prot.name, len(prot_invalid_list), acet.name, len(acet_valid_list))
    if not check_df_name(compared, exp_name):
        PASS = False

    # Get our columns corresponding to the invalid key
    invalid_suffix = "_.*"
    invalid_regex = build_omics_regex(invalid, suffix=invalid_suffix)
    invalid_cols = compared.filter(regex=invalid_regex)
    exp_num_invalid = 1

    # Figure out which columns from proteomics should have been grabbed for the valid genes from prot_invalid_list (should be the same number as number of valid genes
    prot_regex = build_omics_regex(prot_valid_list)
    prot_cols = prot.filter(regex=prot_regex) # Use the regex to get all matching columns
    if len(prot_cols.columns) != len(prot_valid_list):
        print("Unexpected number of matching proteomics columns in test.\n\tExpected: {}\n\tActual: {}".format(len(prot_valid_list), len(prot_cols.columns)))
        PASS = False

    # Figure out which columns from acetylproteomics should have been grabbed for acet_valid_list
    acet_suffix = '-.*'
    acet_regex = build_omics_regex(acet_valid_list, suffix=acet_suffix)
    acet_cols = acet.filter(regex=acet_regex) # Use the regex to get all matching columns

    # Check dataframe shape
    exp_num_rows = len(prot.index.intersection(acet.index))
    exp_num_cols = len(prot_cols.columns) + len(acet_cols.columns) + exp_num_invalid 
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(compared, exp_shape):
        PASS = False

    # Check invalid columns
    if not check_invalid_columns(invalid_cols, exp_num_invalid):
        PASS = False

    # Check columns for valid genes in prot_invalid_list
    if not check_appended_columns(prot, compared, prot_cols.columns):
        PASS = False

    # Check columns for acet_valid_list
    if not check_appended_columns(acet, compared, acet_cols.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_compare_omics_both_list_with_invalid_key():
    """Test that compare_omics will gracefully handle two lists of valid keys, each containing one invalid key."""
    print("Running test_compare_omics_both_list_with_invalid_key...")
    PASS = True

    # Load dataframes to test with
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics()

    # Set the names
    prot_name = "proteomics"
    acet_name = "acetylproteomics"

    # Set the keys
    invalid = "gobbledegook"
    prot_valid_list = ['A4GALT', 'TP53', 'ZSCAN30']
    prot_invalid_list = prot_valid_list.copy()
    prot_invalid_list.append(invalid)
    acet_valid_list = ['AAGAB', 'AACS', 'ZW10', 'ZYX']
    acet_invalid_list = acet_valid_list.copy()
    acet_invalid_list.append(invalid)

    # Run the function, make sure it returned properly
    compared = en.compare_omics(prot_name, acet_name, prot_invalid_list, acet_invalid_list)
    if not check_returned_is_df(compared):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{} for {} genes, with {} for {} genes".format(prot.name, len(prot_invalid_list), acet.name, len(acet_invalid_list))
    if not check_df_name(compared, exp_name):
        PASS = False

    # Get our columns corresponding to the invalid key
    invalid_suffix = "_.*"
    invalid_regex = build_omics_regex(invalid, suffix=invalid_suffix)
    invalid_cols = compared.filter(regex=invalid_regex)
    exp_num_invalid = 2

    # Figure out which columns from proteomics should have been grabbed for the valid genes from prot_invalid_list (should be the same number as number of valid genes
    prot_regex = build_omics_regex(prot_valid_list)
    prot_cols = prot.filter(regex=prot_regex) # Use the regex to get all matching columns
    if len(prot_cols.columns) != len(prot_valid_list):
        print("Unexpected number of matching proteomics columns in test.\n\tExpected: {}\n\tActual: {}".format(len(prot_valid_list), len(prot_cols.columns)))
        PASS = False

    # Figure out which columns from acetylproteomics should have been grabbed for acet_valid_list
    acet_suffix = '-.*'
    acet_regex = build_omics_regex(acet_valid_list, suffix=acet_suffix)
    acet_cols = acet.filter(regex=acet_regex) # Use the regex to get all matching columns

    # Check dataframe shape
    exp_num_rows = len(prot.index.intersection(acet.index))
    exp_num_cols = len(prot_cols.columns) + len(acet_cols.columns) + exp_num_invalid 
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(compared, exp_shape):
        PASS = False

    # Check invalid columns
    if not check_invalid_columns(invalid_cols, exp_num_invalid):
        PASS = False

    # Check columns for valid genes in prot_invalid_list
    if not check_appended_columns(prot, compared, prot_cols.columns):
        PASS = False

    # Check columns for valid genes in acet_invalid_list
    if not check_appended_columns(acet, compared, acet_cols.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_compare_omics_invalid_key_types():
    """Test that compare_omics will gracefully handle a key of an invalid type."""
    print("Running test_compare_omics_invalid_key_types...")
    PASS = True

    # Load our dataframes
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics()

    # Set the names
    prot_name = "proteomics"
    acet_name = "acetylproteomics"

    # Set our keys to use
    prot_valid = 'TP53'
    acet_valid = 'AACS'
    int_key = 100
    prot_valid_list = ['TP53', 'AURKA', 'PIK3R1']

    # Test a key of type int
    comp = en.compare_omics(prot_name, acet_name, prot_valid, int_key)
    if comp is not None:
        print("compare_omics should have returned None when passed a key of type int, but instead returned a {}".format(type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid key message above was expected.)")

    # Test two keys of type int
    comp = en.compare_omics(prot_name, acet_name, int_key, int_key)
    if comp is not None:
        print("compare_omics should have returned None when passed two keys of type int, but instead returned a {}".format(type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid key messages above were expected.)")

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_source_preservation():
    """Test that append_mutations_to_omics does not alter the dataframes it pulls data from."""
    print("Running test_append_mutations_source_preservation...")
    PASS = True

    # Load the source dataframes, set our variables
    mut = en.get_mutations()
    acet = en.get_acetylproteomics()
    acet_name = "acetylproteomics"
    mut_gene = 'TP53'
    mut_genes = ['TP53', 'PIK3R1']
    acet_gene = 'AAGAB'
    acet_genes = ['AACS', 'ZW10']

    # Copy the source dataframes, to compare at the end
    mut_copy = mut.copy()
    acet_copy = acet.copy()

    # Call append_mutations_to_omics a bunch of times
    en.append_mutations_to_omics(acet_name, mut_gene)
    en.append_mutations_to_omics(acet_name, mut_genes)
    en.append_mutations_to_omics(acet_name, mut_gene, acet_gene)
    en.append_mutations_to_omics(acet_name, mut_gene, acet_genes)
    en.append_mutations_to_omics(acet_name, mut_genes, acet_gene)
    en.append_mutations_to_omics(acet_name, mut_genes, acet_genes)
    en.append_mutations_to_omics(acet_name, mut_genes, acet_genes, show_location=False)

    # Check that the source dataframes weren't changed
    if not mut.equals(mut_copy):
        print("Mutations dataframe was altered by append_mutations_to_omics.")
        PASS = False

    if not acet.equals(acet_copy):
        print("Acetylproteomics dataframe was altered by append_mutations_to_omics.")
        PASS = False

    # Indicate whether the test passed
    print_test_result(PASS)

def test_append_mutations_one_mut_one_omics():
    """Test append_mutations_to_omics with one mutation gene and one omics gene."""
    print("Running test_append_mutations_one_mut_one_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics()
    phos_name = "phosphoproteomics"
    phos_gene = 'AAGAB'
    mut_gene = 'PIK3R1'

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(phos_name, mut_gene, omics_genes=phos_gene)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{} for {}, with somatic mutation data for {} gene'.format(phos.name, phos_gene, mut_gene)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Figure out which phosphoproteomics columns should have been grabbed
    phos_suffix = '-.*'
    phos_regex = build_omics_regex(phos_gene, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex) # Use the regex to get all matching columns

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos_cols.columns) + 4
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic_mutation dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_gene):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_one_mut_three_omics():
    """Test append_mutations_to_omics with one mutation gene and three omics genes."""
    print("Running test_append_mutations_one_mut_three_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics()
    phos_name = "phosphoproteomics"
    phos_genes = ['AAGAB', 'AACS', 'ZZZ3']
    mut_gene = 'PIK3R1'

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(phos_name, mut_gene, omics_genes=phos_genes)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{} for {} genes, with somatic mutation data for {} gene'.format(phos.name, len(phos_genes), mut_gene)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Figure out which phosphoproteomics columns should have been grabbed
    phos_suffix = '-.*'
    phos_regex = build_omics_regex(phos_genes, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex) # Use the regex to get all matching columns

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos_cols.columns) + 4
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic_mutation dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_gene):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_one_mut_all_omics():
    """Test append_mutations_to_omics with one mutation gene, and the default parameter of None for the omics gene, which should give the entire omics dataframe."""
    print("Running test_append_mutations_one_mut_all_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics()
    phos_name = "phosphoproteomics"
    mut_gene = 'PIK3R1'

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(phos_name, mut_gene)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{}, with somatic mutation data for {} gene'.format(phos.name, mut_gene)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos.columns) + 4
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic_mutation dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_gene):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_three_mut_one_omics():
    """Test append_mutations_to_omics with three mutation genes and one omics gene."""
    print("Running test_append_mutations_three_mut_one_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics()
    phos_name = "phosphoproteomics"
    phos_gene = 'AAGAB'
    mut_genes = ['PIK3R1', 'TP53', 'AURKA']

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(phos_name, mut_genes, omics_genes=phos_gene)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{} for {}, with somatic mutation data for {} genes'.format(phos.name, phos_gene, len(mut_genes))
    if not check_df_name(appended, exp_name):
        PASS = False

    # Figure out which phosphoproteomics columns should have been grabbed
    phos_suffix = '-.*'
    phos_regex = build_omics_regex(phos_gene, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex) # Use the regex to get all matching columns

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos_cols.columns) + 10
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic_mutation dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_genes):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_three_mut_three_omics():
    """Test append_mutations_to_omics with three mutation genes and three omics genes."""
    print("Running test_append_mutations_three_mut_three_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics()
    phos_name = "phosphoproteomics"
    phos_genes = ['AAGAB', 'AACS', 'ZZZ3']
    mut_genes = ['PIK3R1', 'TP53', 'AURKA']

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(phos_name, mut_genes, omics_genes=phos_genes)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{} for {} genes, with somatic mutation data for {} genes'.format(phos.name, len(phos_genes), len(mut_genes))
    if not check_df_name(appended, exp_name):
        PASS = False

    # Figure out which phosphoproteomics columns should have been grabbed
    phos_suffix = '-.*'
    phos_regex = build_omics_regex(phos_genes, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex) # Use the regex to get all matching columns

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos_cols.columns) + 10
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic_mutation dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_genes):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_three_mut_all_omics():
    """Test append_mutations_to_omics with three mutation genes, and the default parameter of None for the omics gene, which should give the entire omics dataframe."""
    print("Running test_append_mutations_three_mut_all_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics()
    phos_name = "phosphoproteomics"
    mut_genes = ['PIK3R1', 'TP53', 'AURKA']

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(phos_name, mut_genes)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{}, with somatic mutation data for {} genes'.format(phos.name, len(mut_genes))
    if not check_df_name(appended, exp_name):
        PASS = False

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos.columns) + 10
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic_mutation dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_genes):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_one_mut_one_omics_no_location():
    """Test append_mutations_to_omics with one mutation gene and one omics gene, and no location column."""
    print("Running test_append_mutations_one_mut_one_omics_no_location...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics()
    phos_name = "phosphoproteomics"
    phos_gene = 'AAGAB'
    mut_gene = 'PIK3R1'

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(phos_name, mut_gene, omics_genes=phos_gene, show_location=False)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{} for {}, with somatic mutation data for {} gene'.format(phos.name, phos_gene, mut_gene)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Figure out which phosphoproteomics columns should have been grabbed
    phos_suffix = '-.*'
    phos_regex = build_omics_regex(phos_gene, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex) # Use the regex to get all matching columns

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos_cols.columns) + 3
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic_mutation dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_gene, show_location=False):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_one_mut_three_omics_no_location():
    """Test append_mutations_to_omics with one mutation gene and three omics genes, and no location column."""
    print("Running test_append_mutations_one_mut_three_omics_no_location...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics()
    phos_name = "phosphoproteomics"
    phos_genes = ['AAGAB', 'AACS', 'ZZZ3']
    mut_gene = 'PIK3R1'

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(phos_name, mut_gene, omics_genes=phos_genes, show_location=False)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{} for {} genes, with somatic mutation data for {} gene'.format(phos.name, len(phos_genes), mut_gene)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Figure out which phosphoproteomics columns should have been grabbed
    phos_suffix = '-.*'
    phos_regex = build_omics_regex(phos_genes, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex) # Use the regex to get all matching columns

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos_cols.columns) + 3
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic_mutation dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_gene, show_location=False):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_one_mut_all_omics_no_location():
    """Test append_mutations_to_omics with one mutation gene, default of None for omics gene (to select all omics), and no location column."""
    print("Running test_append_mutations_one_mut_all_omics_no_location...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics()
    phos_name = "phosphoproteomics"
    mut_gene = 'PIK3R1'

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(phos_name, mut_gene, show_location=False)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{}, with somatic mutation data for {} gene'.format(phos.name, mut_gene)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos.columns) + 3
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic_mutation dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_gene, show_location=False):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_three_mut_one_omics_no_location():
    """Test append_mutations_to_omics with three mutation genes and one omics gene, and no location column."""
    print("Running test_append_mutations_three_mut_one_omics_no_location...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics()
    phos_name = "phosphoproteomics"
    phos_gene = 'AAGAB'
    mut_genes = ['PIK3R1', 'TP53', 'AURKA']

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(phos_name, mut_genes, omics_genes=phos_gene, show_location=False)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{} for {}, with somatic mutation data for {} genes'.format(phos.name, phos_gene, len(mut_genes))
    if not check_df_name(appended, exp_name):
        PASS = False

    # Figure out which phosphoproteomics columns should have been grabbed
    phos_suffix = '-.*'
    phos_regex = build_omics_regex(phos_gene, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex) # Use the regex to get all matching columns

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos_cols.columns) + 7
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic_mutation dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_genes, show_location=False):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_three_mut_three_omics_no_location():
    """Test append_mutations_to_omics with three mutation genes and three omics genes, and no location column."""
    print("Running test_append_mutations_three_mut_three_omics_no_location...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics()
    phos_name = "phosphoproteomics"
    phos_genes = ['AAGAB', 'AACS', 'ZZZ3']
    mut_genes = ['PIK3R1', 'TP53', 'AURKA']

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(phos_name, mut_genes, omics_genes=phos_genes, show_location=False)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{} for {} genes, with somatic mutation data for {} genes'.format(phos.name, len(phos_genes), len(mut_genes))
    if not check_df_name(appended, exp_name):
        PASS = False

    # Figure out which phosphoproteomics columns should have been grabbed
    phos_suffix = '-.*'
    phos_regex = build_omics_regex(phos_genes, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex) # Use the regex to get all matching columns

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos_cols.columns) + 7
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic_mutation dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_genes, show_location=False):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_three_mut_all_omics_no_location():
    """Test append_mutations_to_omics with three mutation genes, default of None for the omics gene (to select all omics), and no location column."""
    print("Running test_append_mutations_three_mut_all_omics_no_location...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics()
    phos_name = "phosphoproteomics"
    mut_genes = ['PIK3R1', 'TP53', 'AURKA']

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(phos_name, mut_genes, show_location=False)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{}, with somatic mutation data for {} genes'.format(phos.name, len(mut_genes))
    if not check_df_name(appended, exp_name):
        PASS = False

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos.columns) + 7
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic_mutation dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_genes, show_location=False):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_invalid_key():
    """Test that append_mutations_to_omics gracefully handles invalid mutation gene keys."""
    print("Running test_append_mutations_invalid_key...")
    PASS = True

    # Load our dataframe to test with, and set our keys to use
    acet = en.get_acetylproteomics()
    acet_name = "acetylproteomics"
    invalid = 'lorem ipsum'
    invalid_list = ['PIK3R1', 'TAF1', 'GP6', 'lorem ipsum']

    # Test one invalid key
    appended = en.append_mutations_to_omics(acet_name, invalid)
    if appended is not None:
        print("append_mutations_to_omics should have returned None when passed an invalid key, but instead returned a {}".format(type(appended)))
        PASS = False
    else:
        print("(NOTE: The invalid key message above was expected.)")

    # Test one invalid key in a list of valid keys
    appended = en.append_mutations_to_omics(acet_name, invalid_list)
    if appended is not None:
        print("append_mutations_to_omics should have returned None when passed a list of valid keys containing one invalid key, but instead returned a {}".format(type(appended)))
        PASS = False
    else:
        print("(NOTE: The invalid key message above was expected.)")

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_invalid_key_types():
    """Test that append_mutations_to_omics gracefully handles invalid mutation gene key types."""
    print("Running test_append_mutations_invalid_key_types...")
    PASS = True

    # Load our dataframe to test with, and set our keys to use
    prot_name = "proteomics"
    int_key = 44
    dict_key = {0:"TP53", 1:"PIK3R1", 2:"AURKA"}

    # Test a key of type int
    appended = en.append_mutations_to_omics(prot_name, int_key)
    if appended is not None:
        print("append_mutations_to_omics should have returned None when passed a key of invalid type int, but instead returned a {}".format(type(appended)))
        PASS = False
    else:
        print("(NOTE: The invalid key message above was expected.)")

    # Test a key of type dict
    appended = en.append_mutations_to_omics(prot_name, dict_key)
    if appended is not None:
        print("append_mutations_to_omics should have returned None when passed a key of invalid type dict, but instead returned a {}".format(type(appended)))
        PASS = False
    else:
        print("(NOTE: The invalid key message above was expected.)")

    # Print whether the test passed
    print_test_result(PASS)

def test_append_metadata_source_preservation():
    """Test that append_metadata_to_omics does not alter the dataframes it pulls data from."""
    print("Running test_append_metadata_source_preservation...")
    PASS = True

    # Load the source dataframes, and set our variables
    clin = en.get_clinical()
    derived_mol = en.get_derived_molecular()
    exp_setup = en.get_experimental_setup()
    phos = en.get_phosphoproteomics()
    clin_name = "clinical"
    derived_mol_name = "derived_molecular"
    exp_setup_name = "experimental_setup"
    phos_name = "phosphoproteomics"

    clin_col = "Country"
    clin_cols = ["Country", "tumor_Stage-Pathological", "LVSI"]
    derived_mol_col = "CIBERSORT_Eosinophils" 
    derived_mol_cols = ["CIBERSORT_Eosinophils", "Pathway_activity_JAK.STAT", "Progesterone_Receptor_%"]
    exp_setup_col = "RNAseq_R1_UUID"
    exp_setup_cols = ["RNAseq_R1_UUID", "Methylation_available", "WXS_tumor_UUID"]
    phos_col = "ZZZ3"
    phos_cols = ["AAAS", "AAED1"]

    # Copy the source dataframes, to compare at the end
    clin_copy = clin.copy()
    derived_mol_copy = derived_mol.copy()
    exp_setup_copy = exp_setup.copy()
    phos_copy = phos.copy()

    # Call append_metadata_to_omics a bunch of times
    en.append_metadata_to_omics(clin_name, phos_name)
    en.append_metadata_to_omics(clin_name, phos_name, metadata_cols=clin_col)
    en.append_metadata_to_omics(clin_name, phos_name, omics_genes=phos_col)
    en.append_metadata_to_omics(clin_name, phos_name, metadata_cols=clin_col, omics_genes=phos_col)
    en.append_metadata_to_omics(clin_name, phos_name, metadata_cols=clin_cols)
    en.append_metadata_to_omics(clin_name, phos_name, omics_genes=phos_cols)
    en.append_metadata_to_omics(clin_name, phos_name, metadata_cols=clin_cols, omics_genes=phos_col)
    en.append_metadata_to_omics(clin_name, phos_name, metadata_cols=clin_col, omics_genes=phos_cols)
    en.append_metadata_to_omics(clin_name, phos_name, metadata_cols=clin_cols, omics_genes=phos_cols)

    en.append_metadata_to_omics(derived_mol_name, phos_name)
    en.append_metadata_to_omics(derived_mol_name, phos_name, metadata_cols=derived_mol_col)
    en.append_metadata_to_omics(derived_mol_name, phos_name, omics_genes=phos_col)
    en.append_metadata_to_omics(derived_mol_name, phos_name, metadata_cols=derived_mol_col, omics_genes=phos_col)
    en.append_metadata_to_omics(derived_mol_name, phos_name, metadata_cols=derived_mol_cols)
    en.append_metadata_to_omics(derived_mol_name, phos_name, omics_genes=phos_cols)
    en.append_metadata_to_omics(derived_mol_name, phos_name, metadata_cols=derived_mol_cols, omics_genes=phos_col)
    en.append_metadata_to_omics(derived_mol_name, phos_name, metadata_cols=derived_mol_col, omics_genes=phos_cols)
    en.append_metadata_to_omics(derived_mol_name, phos_name, metadata_cols=derived_mol_cols, omics_genes=phos_cols)

    en.append_metadata_to_omics(exp_setup_name, phos_name)
    en.append_metadata_to_omics(exp_setup_name, phos_name, metadata_cols=exp_setup_col)
    en.append_metadata_to_omics(exp_setup_name, phos_name, omics_genes=phos_col)
    en.append_metadata_to_omics(exp_setup_name, phos_name, metadata_cols=exp_setup_col, omics_genes=phos_col)
    en.append_metadata_to_omics(exp_setup_name, phos_name, metadata_cols=exp_setup_cols)
    en.append_metadata_to_omics(exp_setup_name, phos_name, omics_genes=phos_cols)
    en.append_metadata_to_omics(exp_setup_name, phos_name, metadata_cols=exp_setup_cols, omics_genes=phos_col)
    en.append_metadata_to_omics(exp_setup_name, phos_name, metadata_cols=exp_setup_col, omics_genes=phos_cols)
    en.append_metadata_to_omics(exp_setup_name, phos_name, metadata_cols=exp_setup_cols, omics_genes=phos_cols)

    # Check that the source dataframes weren't changed
    if not clin.equals(clin_copy):
        print("clinical dataframe was altered by append_metadata_to_omics.")
        PASS = False

    if not derived_mol.equals(derived_mol_copy):
        print("derived_molecular dataframe was altered by append_metadata_to_omics.")
        PASS = False

    if not exp_setup.equals(exp_setup_copy):
        print("experimental_setup dataframe was altered by append_metadata_to_omics.")
        PASS = False

    if not phos.equals(phos_copy):
        print("phosphoproteomics dataframe was altered by append_metadata_to_omics.")
        PASS = False

    # Indicate whether the test passed
    print_test_result(PASS)

# One meta one omics
def test_append_metadata_one_meta_one_omics():
    """Test append_metadata_to_omics with one metadata column and one omics gene."""
    print("Running test_append_metadata_one_meta_one_omics...")
    PASS = True

    # Load the source dataframes
    derived_mol = en.get_derived_molecular()
    phos = en.get_phosphoproteomics()

    # Set the names
    derived_mol_name = "derived_molecular"
    phos_name = "phosphoproteomics"

    # Run the function, make sure it returned properly
    derived_mol_col = "Purity_Stroma"
    phos_gene = "USP36"
    appended = en.append_metadata_to_omics(derived_mol_name, phos_name, metadata_cols=derived_mol_col, omics_genes=phos_gene)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip other tests, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{} from {}, with {} for {}".format(derived_mol_col, derived_mol.name, phos.name, phos_gene)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Get the columns that should've been selected from phosphoproteomics
    phos_suffix = "-.*"
    phos_regex = build_omics_regex(phos_gene, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex)

    # Check dataframe shape
    exp_num_rows = len(derived_mol.index.intersection(phos.index))
    exp_num_cols = len(phos_cols.columns) + 1
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check column values
    if not check_appended_columns(derived_mol, appended, derived_mol_col, derived_mol_col):
        PASS = False

    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

# One meta three omics
def test_append_metadata_one_meta_three_omics():
    """Test append_metadata_to_omics with one metadata column and three omics genes."""
    print("Running test_append_metadata_one_meta_three_omics...")
    PASS = True

    # Load the source dataframes
    derived_mol = en.get_derived_molecular()
    phos = en.get_phosphoproteomics()

    # Set the names
    derived_mol_name = "derived_molecular"
    phos_name = "phosphoproteomics"

    # Run the function, make sure it returned properly
    derived_mol_col = "Purity_Stroma"
    phos_genes = ["USP36", "TMEM209", "STXBP5"]
    appended = en.append_metadata_to_omics(derived_mol_name, phos_name, metadata_cols=derived_mol_col, omics_genes=phos_genes)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip other tests, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{} from {}, with {} for {} genes".format(derived_mol_col, derived_mol.name, phos.name, len(phos_genes))
    if not check_df_name(appended, exp_name):
        PASS = False

    # Get the columns that should've been selected from phosphoproteomics
    phos_suffix = "-.*"
    phos_regex = build_omics_regex(phos_genes, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex)

    # Check dataframe shape
    exp_num_rows = len(derived_mol.index.intersection(phos.index))
    exp_num_cols = len(phos_cols.columns) + 1
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check column values
    if not check_appended_columns(derived_mol, appended, derived_mol_col, derived_mol_col):
        PASS = False

    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

# One meta all omics
def test_append_metadata_one_meta_all_omics():
    """Test append_metadata_to_omics with one metadata column, and the default parameter of None for the omics gene, which should cause the entire omics dataframe to be selected."""
    print("Running test_append_metadata_one_meta_all_omics...")
    PASS = True

    # Load the source dataframes
    derived_mol = en.get_derived_molecular()
    phos = en.get_phosphoproteomics()

    # Set the names
    derived_mol_name = "derived_molecular"
    phos_name = "phosphoproteomics"

    # Run the function, make sure it returned properly
    derived_mol_col = "Purity_Stroma"
    appended = en.append_metadata_to_omics(derived_mol_name, phos_name, metadata_cols=derived_mol_col)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip other tests, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{} from {}, with {}".format(derived_mol_col, derived_mol.name, phos.name)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Check dataframe shape
    exp_num_rows = len(derived_mol.index.intersection(phos.index))
    exp_num_cols = len(phos.columns) + 1
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check column values
    if not check_appended_columns(derived_mol, appended, derived_mol_col, derived_mol_col):
        PASS = False

    if not check_appended_columns(phos, appended, phos.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

# Three meta one omics
def test_append_metadata_three_meta_one_omics():
    """Test append_metadata_to_omics with three metadata columns and one omics gene."""
    print("Running test_append_metadata_three_meta_one_omics...")
    PASS = True

    # Load the source dataframes
    derived_mol = en.get_derived_molecular()
    phos = en.get_phosphoproteomics()

    # Set the names
    derived_mol_name = "derived_molecular"
    phos_name = "phosphoproteomics"

    # Run the function, make sure it returned properly
    derived_mol_cols = ["Purity_Stroma", "POLE_subtype", "CIBERSORT_T _cells _CD4 _memory _resting"]
    phos_gene = "USP36"
    appended = en.append_metadata_to_omics(derived_mol_name, phos_name, metadata_cols=derived_mol_cols, omics_genes=phos_gene)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip other tests, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{} columns from {}, with {} for {}".format(len(derived_mol_cols), derived_mol.name, phos.name, phos_gene)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Get the columns that should've been selected from phosphoproteomics
    phos_suffix = "-.*"
    phos_regex = build_omics_regex(phos_gene, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex)

    # Check dataframe shape
    exp_num_rows = len(derived_mol.index.intersection(phos.index))
    exp_num_cols = len(phos_cols.columns) + len(derived_mol_cols)
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check column values
    if not check_appended_columns(derived_mol, appended, derived_mol_cols, derived_mol_cols):
        PASS = False

    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

# Three meta three omics
def test_append_metadata_three_meta_three_omics():
    """Test append_metadata_to_omics with three metadata columns and three omics genes."""
    print("Running test_append_metadata_three_meta_three_omics...")
    PASS = True

    # Load the source dataframes
    derived_mol = en.get_derived_molecular()
    phos = en.get_phosphoproteomics()

    # Set the names
    derived_mol_name = "derived_molecular"
    phos_name = "phosphoproteomics"

    # Run the function, make sure it returned properly
    derived_mol_cols = ["Purity_Stroma", "POLE_subtype", "CIBERSORT_T _cells _CD4 _memory _resting"]
    phos_genes = ["USP36", "TMEM209", "STXBP5"]
    appended = en.append_metadata_to_omics(derived_mol_name, phos_name, metadata_cols=derived_mol_cols, omics_genes=phos_genes)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip other tests, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{} columns from {}, with {} for {} genes".format(len(derived_mol_cols), derived_mol.name, phos.name, len(phos_genes))
    if not check_df_name(appended, exp_name):
        PASS = False

    # Get the columns that should've been selected from phosphoproteomics
    phos_suffix = "-.*"
    phos_regex = build_omics_regex(phos_genes, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex)

    # Check dataframe shape
    exp_num_rows = len(derived_mol.index.intersection(phos.index))
    exp_num_cols = len(phos_cols.columns) + len(derived_mol_cols)
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check column values
    if not check_appended_columns(derived_mol, appended, derived_mol_cols, derived_mol_cols):
        PASS = False

    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

# Three meta all omics
def test_append_metadata_three_meta_all_omics():
    """Test append_metadata_to_omics with three metadata columns, and the default parameter of None for the omics gene, which should cause the entire omics dataframe to be selected."""
    print("Running test_append_metadata_three_meta_all_omics...")
    PASS = True

    # Load the source dataframes
    derived_mol = en.get_derived_molecular()
    phos = en.get_phosphoproteomics()

    # Set the names
    derived_mol_name = "derived_molecular"
    phos_name = "phosphoproteomics"

    # Run the function, make sure it returned properly
    derived_mol_cols = ["Purity_Stroma", "POLE_subtype", "CIBERSORT_T _cells _CD4 _memory _resting"]
    appended = en.append_metadata_to_omics(derived_mol_name, phos_name, metadata_cols=derived_mol_cols)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip other tests, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{} columns from {}, with {}".format(len(derived_mol_cols), derived_mol.name, phos.name)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Check dataframe shape
    exp_num_rows = len(derived_mol.index.intersection(phos.index))
    exp_num_cols = len(phos.columns) + len(derived_mol_cols)
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check column values
    if not check_appended_columns(derived_mol, appended, derived_mol_cols, derived_mol_cols):
        PASS = False

    if not check_appended_columns(phos, appended, phos.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

# All meta one omics
def test_append_metadata_all_meta_one_omics():
    """Test append_metadata_to_omics with the default of None for metadata_cols, which should select the entire dataframe, and one omics gene."""
    print("Running test_append_metadata_all_meta_one_omics...")
    PASS = True

    # Load the source dataframes
    derived_mol = en.get_derived_molecular()
    phos = en.get_phosphoproteomics()

    # Set the names
    derived_mol_name = "derived_molecular"
    phos_name = "phosphoproteomics"

    # Run the function, make sure it returned properly
    phos_gene = "USP36"
    appended = en.append_metadata_to_omics(derived_mol_name, phos_name, omics_genes=phos_gene)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip other tests, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{}, with {} for {}".format(derived_mol.name, phos.name, phos_gene)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Get the columns that should've been selected from phosphoproteomics
    phos_suffix = "-.*"
    phos_regex = build_omics_regex(phos_gene, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex)

    # Check dataframe shape
    exp_num_rows = len(derived_mol.index.intersection(phos.index))
    exp_num_cols = len(phos_cols.columns) + len(derived_mol.columns)
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check column values
    if not check_appended_columns(derived_mol, appended, derived_mol.columns, derived_mol.columns):
        PASS = False

    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

# All meta three omics
def test_append_metadata_all_meta_three_omics():
    """Test append_metadata_to_omics with the default of None for metadata_cols, which should select the entire dataframe, and three omics genes."""
    print("Running test_append_metadata_all_meta_three_omics...")
    PASS = True

    # Load the source dataframes
    derived_mol = en.get_derived_molecular()
    phos = en.get_phosphoproteomics()

    # Set the names
    derived_mol_name = "derived_molecular"
    phos_name = "phosphoproteomics"

    # Run the function, make sure it returned properly
    phos_genes = ["USP36", "TMEM209", "STXBP5"]
    appended = en.append_metadata_to_omics(derived_mol_name, phos_name, omics_genes=phos_genes)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip other tests, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{}, with {} for {} genes".format(derived_mol.name, phos.name, len(phos_genes))
    if not check_df_name(appended, exp_name):
        PASS = False

    # Get the columns that should've been selected from phosphoproteomics
    phos_suffix = "-.*"
    phos_regex = build_omics_regex(phos_genes, suffix=phos_suffix)
    phos_cols = phos.filter(regex=phos_regex)

    # Check dataframe shape
    exp_num_rows = len(derived_mol.index.intersection(phos.index))
    exp_num_cols = len(phos_cols.columns) + len(derived_mol.columns)
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check column values
    if not check_appended_columns(derived_mol, appended, derived_mol.columns, derived_mol.columns):
        PASS = False

    if not check_appended_columns(phos, appended, phos_cols.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

# All meta all omics (default parameters)
def test_append_metadata_default_parameters():
    """Test append_metadata_to_omics with the default parameters of None for metadata_cols and omics_genes, which should cause it to select the entire metadata and omics dataframes."""
    print("Running test_append_metadata_default_parameters...")
    PASS = True

    # Load the source dataframes
    derived_mol = en.get_derived_molecular()
    phos = en.get_phosphoproteomics()

    # Set the names
    derived_mol_name = "derived_molecular"
    phos_name = "phosphoproteomics"

    # Run the function, make sure it returned properly
    appended = en.append_metadata_to_omics(derived_mol_name, phos_name)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip other tests, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = "{}, with {}".format(derived_mol.name, phos.name)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Check dataframe shape
    exp_num_rows = len(derived_mol.index.intersection(phos.index))
    exp_num_cols = len(phos.columns) + len(derived_mol.columns)
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check column values
    if not check_appended_columns(derived_mol, appended, derived_mol.columns, derived_mol.columns):
        PASS = False

    if not check_appended_columns(phos, appended, phos.columns):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)


# All valid dfs (omics or meta)
def test_append_metadata_all_dfs():
    """Test that append_metadata_to_omics works with all dataframes that are valid for the function."""

    # Load our dataframes to test, and set the keys we'll use.
    clin = en.get_clinical()
    derived_mol = en.get_derived_molecular()
    exp_setup = en.get_exp_setup()

    acet = en.get_acetylproteomics()
    CNA = en.get_CNA()
    phosg = en.get_phosphoproteomics_gene()
    phoss = en.get_phosphoproteomics()
    prot = en.get_proteomics()
    tran = en.get_transcriptomics()

    # Call append_metadata_to_omics on the dataframes

    # Check the return values

    # Print whether the test passed

# Invalid dfs

# Invalid keys

# Invalid key types

# ADD ALL DFS TEST FOR MUT

en = cptac.Endometrial()

print("\nRunning tests:\n")
 
print("Testing getters...")
test_get_clinical()
test_get_derived_molecular()
test_get_experimental_setup()
test_get_acetylproteomics()
test_get_proteomics()
test_get_transcriptomics()
test_get_circular_RNA()
test_get_miRNA()
test_get_CNA()
test_get_phosphoproteomics()
test_get_phosphoproteomics_gene()
test_get_phosphosites()
test_get_mutations()
test_get_mutations_binary()

print("\nTesting compare and append functions...")
test_compare_omics_source_preservation()
test_compare_omics_default_parameters()
test_compare_omics_one_gene()
test_compare_omics_multiple_genes()
test_compare_omics_all_dfs()
test_compare_omics_invalid_dfs()
test_compare_omics_one_invalid_key()
test_compare_omics_both_invalid_keys()
test_compare_omics_one_list_with_invalid_key()
test_compare_omics_both_list_with_invalid_key()
test_compare_omics_invalid_key_types() 

test_append_metadata_source_preservation()
test_append_metadata_one_meta_one_omics()
test_append_metadata_one_meta_three_omics()
test_append_metadata_one_meta_all_omics()
test_append_metadata_three_meta_one_omics()
test_append_metadata_three_meta_three_omics()
test_append_metadata_three_meta_all_omics()
test_append_metadata_all_meta_one_omics()
test_append_metadata_all_meta_three_omics()
test_append_metadata_default_parameters()

test_append_mutations_source_preservation()
test_append_mutations_one_mut_one_omics()
test_append_mutations_one_mut_three_omics()
test_append_mutations_one_mut_all_omics()
test_append_mutations_three_mut_one_omics()
test_append_mutations_three_mut_three_omics()
test_append_mutations_three_mut_all_omics()
test_append_mutations_one_mut_all_omics_no_location()
test_append_mutations_one_mut_one_omics_no_location()
test_append_mutations_one_mut_three_omics_no_location()
test_append_mutations_three_mut_one_omics_no_location()
test_append_mutations_three_mut_three_omics_no_location()
test_append_mutations_three_mut_all_omics_no_location()
test_append_mutations_invalid_key()

print("Version:", cptac.version())
