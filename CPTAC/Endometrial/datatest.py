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
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import CPTAC.Endometrial as en
from utilities import *

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

def check_getter(df, exp_name, exp_dim, exp_headers, coordinates, values): # private
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
    """Builds a regex from a list of genes to grab all columns corresponding to those genes from the acetylproteomics or phosphoproteomics dataframes.

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

def check_appended_column(source_df, dest_df, source_header, dest_header): # private
    """Checks whether a column appended to a dataframe has the same values for each index as it has in the dataframe it was taken from.

    Parameters:
    source_df (pandas.core.frame.DataFrame): the dataframe the column was taken from.
    dest_df (pandas.core.frame.DataFrame): the dataframe the column was appended to (with the column appended to it).
    source_header (str): the column's header in the source dataframe.
    dest_header (str): the column's header in the dataframe it was appended to.

    Returns:
    bool: Indicates whether the column in the destination dataframe and the column in the source dataframe had the same values for each index.
    """
    PASS = True

    original_column = source_df[source_header]
    merged_column = dest_df[dest_header]

    # Sort columns. If they're in different orders in the dataframes, that's fine as long as the right samples are with the right data. But we do need them in the same order for series.equal.
    original_column = original_column.sort_index()
    merged_column = merged_column.sort_index()

    if not original_column.equals(merged_column):
        PASS = False
        for idx in merged_column.index:
            merged_value = merged_column.loc[idx]
            original_value = original_column.loc[idx]
            if (merged_value != original_value) and (pd.notna(merged_value) or pd.notna(original_value)):
                print("Merged dataframe had incorrect values.\n\tSample: {}\tColumn: {}\n\tExpected: {}\tActual: {}\n".format(idx, dest_header, original_value, merged_value))
    return PASS

def check_appended_columns(source_df, dest_df, headers):
    """Checks whether a column or list of columns appended to a dataframe have the same values for each index in that dataframe as they had in the dataframe they were taken from.

    Parameters:
    source_df (pandas.core.frame.DataFrame): The dataframe the columns were taken from.
    dest_df (pandas.core.frame.DataFrame): The dataframe the columns were appended to (with them appended to it).
    header (str or list or pandas.core.indexes.base.Index): The header(s) of the columns in source_df. str if single, list or index of str if multiple. Header(s) in dest_df will be constructed from this, and the name of source_df.
    
    Returns:
    bool: Indicates whether the specified column(s) in dest_df had the same values for each index as they did in source_df.
    """
    PASS = True
    if isinstance(headers, str):
        headers = [headers]
    for source_header in headers:
        dest_header = source_header + '_' + source_df.name
        if not check_appended_column(source_df, dest_df, source_header, dest_header):
            PASS = False
    return PASS

def check_mutation_columns_single_gene(mutations, merged_df, gene, multiple_mutations, show_location):
    """
    Parameters
    mutations (pandas.core.frame.DataFrame): The somatic mutations dataframe.
    merged_df (pandas.core.frame.DataFrame): The merged datframe.
    gene (str): The gene the mutation data was collected for.
    multiple_mutations (bool): Whether multiple mutations for one sample were included in the column.
    show_location (bool): Whether the location column was included in merged_df.

    Returns
    bool: Indicates whether the mutation data for that gene and each sample in the merged dataframe matched the data in the somatic mutations dataframe.
    """
    PASS = True

    # Set our column names for use later
    patient_col = 'Clinical_Patient_Key'
    gene_col = 'Gene'
    location_col = 'Location'
    mutation_col = 'Mutation'
    status_col = 'Sample_Status'
    merged_location_col = gene + '_' + location_col
    merged_mutation_col = gene + '_' + mutation_col

    # Load all the mutations for the gene
    gene_df = mutations.loc[mutations[gene_col] == gene] 
    gene_df = gene_df.set_index(patient_col)

    for sample in merged_df.index.values:
        # Get the rows for just this sample from our two dataframes
        sample_df = gene_df.loc[gene_df.index == sample]
        merged_sample_df = merged_df.loc[merged_df.index == sample]

        if len(sample_df.index) == 0: # There were no mutations for that gene in this sample
            original_location = 'No_mutation'
            if sample <= 'S104':
                original_mutation = 'Wildtype_Tumor'
            else:
                original_mutation = 'Wildtype_Normal'
            sample_dict = { # Create a prep dictionary with what the values for the different columns should be
                patient_col:[sample], # This will be our index
                location_col:original_location,
                mutation_col:original_mutation} 
            sample_df = pd.DataFrame(data=sample_dict, index=sample_dict[patient_col]) # Make that dict a dataframe, and set it as our sample_df

        elif len(sample_df.index) == 1: # There was one mutation for that gene in this sample
            pass # sample_df and merged_sample_df are already ready for the check_appended_column calls

        else: # There were multiple mutations for that gene in this sample
            if multiple_mutations:
                # Since there are multiple rows in our two dataframes, sort by location for consistency
                sample_df = sample_df.sort_values(by=location_col)
                merged_sample_df = merged_sample_df.sort_values(by=merged_location_col) 
            else:
                # Get just the mutation we care about
                sample_df = Utilities().add_mutation_hierarchy(sample_df)
                sample_df = sample_df.sort_values(by=[patient_col, 'Mutation_Hierarchy'], ascending=[True,False]) # Sorts by patient key, then by hierarchy so the duplicates will come with the lower number first
                sample_df = sample_df.iloc[[0]]

        # Test our location and mutation columns
        if show_location: # Only try to test the location if it was included haha
            if not check_appended_column(sample_df, merged_sample_df, location_col, merged_location_col):
                PASS = False
        if not check_appended_column(sample_df, merged_sample_df, mutation_col, merged_mutation_col):
            PASS = False

        # Test our Sample_Status column
        if sample <= 'S104': # Figure out what our sample status should be
            original_status = 'Tumor'
        else:
            original_status = 'Normal'
        sample_df = sample_df.assign(**{status_col:original_status}) # Append a Sample_Status column to our original values dataframe, with what the value should be, for comparison to the merged dataframe.
        if not check_appended_column(sample_df, merged_sample_df, status_col, status_col):
            PASS = False

    return PASS

def check_mutation_columns(mutations, merged_df, genes, multiple_mutations=False, show_location=True):
    """
    Parameters
    mutations (pandas.core.frame.DataFrame): The somatic mutations dataframe.
    merged_df (pandas.core.frame.DataFrame): The merged datframe.
    gene (str or list): The gene(s) the mutation data was collected for. str if one, list of str if multiple.
    multiple_mutations (bool, optional): Whether multiple mutations for one sample were included in the column. Default is False.
    show_location (bool, optional): Whether the location column was included in merged_df. Default is True.

    Returns
    bool: Indicates whether the mutation data for that gene and each sample in the merged dataframe matched the data in the somatic mutations dataframe.
    """
    PASS = True
    if isinstance(genes, str):
        genes = [genes]
    for gene in genes:
        if not check_mutation_columns_single_gene(mutations, merged_df, gene, multiple_mutations, show_location):
            PASS = False
    return PASS
    
# Test functions that get dataframes
def test_get_clinical_filtered():
    """Test get_clinical with the default parameter unfiltered=False."""

    print('Running test_get_clinical with the default parameter unfiltered=False...')

    df = en.get_clinical()
    name = "clinical"
    dimensions = (144, 26)
    headers = ['Proteomics_Participant_ID', 'Proteomics_Tumor_Normal', 'Country', 'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity', 'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site', 'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm', 'Num_full_term_pregnancies']
    test_coord = ((79, 16), (15, 25), (88, 2))
    test_vals = (77.0, '3', 'Poland')

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_clinical_unfiltered():
    """Test get_clinical with parameter unfiltered=True."""

    print('Running test_get_clinical with parameter unfiltered=True...')

    df = en.get_clinical(unfiltered=True)
    print("(NOTE: The unfiltered data warning above was expected.)") # To avoid confusion

    name = "clinical"
    dimensions = (153, 27)
    headers = ['Proteomics_Participant_ID', 'Case_excluded', 'Proteomics_Tumor_Normal', 'Country', 'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity', 'Path_Stage_Primary_Tumor-pT', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site', 'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm', 'Num_full_term_pregnancies']
    test_coord = ((23, 8), (151, 1), (32, 26))
    test_vals = ('Normal', 'No', '3')

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_derived_molecular_filtered():
    """Test get_derived_molecular with default parameter unfiltered=False."""

    print('Running test_get_derived_molecular with default parameter unfiltered=False...')

    df = en.get_derived_molecular()
    name = 'derived_molecular'
    dimensions = (144, 144) 
    headers = ['Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs', 'Proteomics_Aliquot_ID', 'Proteomics_OCT', 'Estrogen_Receptor', 'Estrogen_Receptor_%', 'Progesterone_Receptor', 'Progesterone_Receptor_%', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID', 'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality']
    test_coord = ((2, 3), (90, 143), (143, 4))
    test_vals = ('C3L-00032-01', 'PASS', 'CPT0230460002,CPT0230460003,CPT0230460004,CPT0230470002,CPT0230470003,CPT0230470004,CPT0230480002,CPT0230480003,CPT0230480004')

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_derived_molecular_unfiltered():
    """Test get_derived_molecular with parameter unfiltered=True."""

    print('Running test_get_derived_molecular with parameter unfiltered=True...')

    df = en.get_derived_molecular(unfiltered=True)
    print("(NOTE: The unfiltered data warning above was expected.)") # To avoid confusion

    name = 'derived_molecular'
    dimensions = (153, 144)
    headers = ['Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs', 'Proteomics_Aliquot_ID', 'Proteomics_OCT', 'Estrogen_Receptor', 'Estrogen_Receptor_%', 'Progesterone_Receptor', 'Progesterone_Receptor_%', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID', 'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality']
    test_coord = ((152, 2), (4, 143), (30, 60))
    test_vals = ('130N', 'PASS', -0.13)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_acetylproteomics_filtered():
    """Test get_acetylproteomics with default parameter unfiltered=False."""

    print('Running test_get_acetylproteomics with default parameter unfiltered=False...')

    df = en.get_acetylproteomics()
    name = 'acetylproteomics'
    dimensions = (144, 10862)
    headers = ['A2M-K1168', 'A2M-K1176', 'A2M-K135', 'A2M-K145', 'A2M-K516', 'A2M-K664', 'A2M-K682', 'AACS-K391', 'AAGAB-K290', 'AAK1-K201', 'ZSCAN31-K215', 'ZSCAN32-K659', 'ZW10-K634', 'ZYX-K24', 'ZYX-K25', 'ZYX-K265', 'ZYX-K272', 'ZYX-K279', 'ZYX-K533', 'ZZZ3-K117']
    test_coord = ((1, 1), (12, 10861), (90, 5849))
    test_vals = (0.47700000000000004, 0.16, 0.4098)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_acetylproteomics_unfiltered():
    """Test get_acetylproteomics with parameter unfiltered=True."""

    print('Running test_get_acetylproteomics with parameter unfiltered=True...')

    df = en.get_acetylproteomics(unfiltered=True)
    print("(NOTE: The unfiltered data warning above was expected.)") # To avoid confusion

    name = 'acetylproteomics'
    dimensions = (153, 10862)
    headers = ['A2M-K1168', 'A2M-K1176', 'A2M-K135', 'A2M-K145', 'A2M-K516', 'A2M-K664', 'A2M-K682', 'AACS-K391', 'AAGAB-K290', 'AAK1-K201', 'ZSCAN31-K215', 'ZSCAN32-K659', 'ZW10-K634', 'ZYX-K24', 'ZYX-K25', 'ZYX-K265', 'ZYX-K272', 'ZYX-K279', 'ZYX-K533', 'ZZZ3-K117']
    test_coord = ((1, 1), (15, 10861), (90, 4399))
    test_vals = (0.47700000000000004, 0.16, 0.6920000000000001)

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

def test_get_transcriptomics_linear():
    """Test get_transcriptomics_linear."""

    print('Running test_get_transcriptomics_linear...')

    df = en.get_transcriptomics_linear()
    name = "transcriptomics_linear"
    dimensions = (109, 28057)
    headers = ['A1BG', 'A1BG-AS1', 'A1CF', 'A2M', 'A2M-AS1', 'A2ML1', 'A2MP1', 'A3GALT2', 'A4GALT', 'A4GNT', 'ZWILCH', 'ZWINT', 'ZXDA', 'ZXDB', 'ZXDC', 'ZYG11A', 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3']
    test_coord = ((22, 25483), (108, 23), (101, 17748))
    test_vals = (0.82, 12.0, 6.19)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_transcriptomics_circular():
    """Test get_transcriptomics_circular."""

    print('Running test_get_transcriptomics_circular...')

    df = en.get_transcriptomics_circular()
    name = "transcriptomics_circular"
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

def test_get_cna():
    """Test get_cna."""

    print('Running test_get_cna...')

    df = en.get_cna()
    name = "cna"
    dimensions = (95, 28057)
    headers = ['MFSD14A', 'SASS6', 'TRMT13', 'LRRC39', 'DBT', 'RTCA-AS1', 'RTCA', 'MIR553', 'UBE4B', 'CDC14A', 'TSPY8', 'FAM197Y2', 'FAM197Y4', 'FAM197Y5', 'FAM197Y7', 'FAM197Y8', 'FAM197Y6', 'FAM197Y3', 'RBMY3AP', 'TTTY22']
    test_coord = ((12, 27865), (60, 8), (94, 15439))
    test_vals = (-0.07, 0.01, 0.03)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_phosphoproteomics_site():
    """Test get_phosphoproteomics_site."""

    print('Running test_get_phosphoproteomics_site...')

    df =  en.get_phosphoproteomics_site()
    name = "phosphoproteomics_site"
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
    name = 'phosphoproteomics_site for ' + gene
    dimensions = (144, 37)
    headers = ['AAK1-S14_phosphoproteomics_site', 'AAK1-S18_phosphoproteomics_site', 'AAK1-S20_phosphoproteomics_site', 'AAK1-S21_phosphoproteomics_site', 'AAK1-S26_phosphoproteomics_site', 'AAK1-S618_phosphoproteomics_site', 'AAK1-S623_phosphoproteomics_site', 'AAK1-S624_phosphoproteomics_site', 'AAK1-S637_phosphoproteomics_site', 'AAK1-S642_phosphoproteomics_site', 'AAK1-T448_phosphoproteomics_site', 'AAK1-T606_phosphoproteomics_site', 'AAK1-T620_phosphoproteomics_site', 'AAK1-T640_phosphoproteomics_site', 'AAK1-T653_phosphoproteomics_site', 'AAK1-T674_phosphoproteomics_site', 'AAK1-T681_phosphoproteomics_site', 'AAK1-T694_phosphoproteomics_site', 'AAK1-T848_phosphoproteomics_site', 'AAK1-Y687_phosphoproteomics_site']
    test_coord = ((5, 33), (64, 14), (128, 0))
    test_vals = (0.547, -0.5379999999999999, 0.1395)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_mutations_maf():
    """Test get_mutations_maf."""

    print('Running test_get_mutations_maf...')

    df = en.get_mutations_maf()
    name = "somatic MAF"
    dimensions = (52560, 5)
    headers = ['Clinical_Patient_Key', 'Patient_Id', 'Gene', 'Mutation', 'Location']
    test_coord = ((52000, 3), (12, 4), (34567, 0))
    test_vals = ('Missense_Mutation', 'p.T2121P', 'S059')

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_mutations_binary():
    """Test get_mutations_binary."""

    print('Running test_get_mutations_binary...')

    df = en.get_mutations_binary()
    name = "somatic binary"
    dimensions = (95, 51559)
    headers = ['A1BG_p.E298K', 'A1BG_p.S181N', 'A1CF_p.F487L', 'A1CF_p.S236Y', 'A2ML1_p.A8V', 'A2ML1_p.G1306D', 'A2ML1_p.L1347F', 'A2ML1_p.L82I', 'A2ML1_p.P712S', 'A2ML1_p.R443Q', 'ZYG11A_p.Q442H', 'ZYG11B_p.H315R', 'ZYG11B_p.R495M', 'ZYG11B_p.R728C', 'ZYX_p.C447Y', 'ZZEF1_p.A2723V', 'ZZEF1_p.D845Y', 'ZZEF1_p.K1251E', 'ZZEF1_p.K2387Sfs*40', 'ZZZ3_p.Y891C']
    test_coord = ((94, 51558), (0, 0), (45, 25436))
    test_vals = (0, 0, 0)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_mutations_unparsed():
    """Test get_mutations_unparsed."""

    print('Running test_get_mutations_unparsed...')

    df = en.get_mutations_unparsed()
    name = "somatic MAF unparsed"
    dimensions = (53101, 124)
    headers = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'ExAC_AC_AN_Adj', 'ExAC_AC_AN', 'ExAC_AC_AN_AFR', 'ExAC_AC_AN_AMR', 'ExAC_AC_AN_EAS', 'ExAC_AC_AN_FIN', 'ExAC_AC_AN_NFE', 'ExAC_AC_AN_OTH', 'ExAC_AC_AN_SAS', 'ExAC_FILTER']
    test_coord = ((52265, 45), (12, 70), (27658, 1))
    test_vals = ('strelkasnv-varssnv-mutectsnv', 'UPI0000167B91', 0)

    PASS = check_getter(df, name, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

# Test merging and appending functions
def test_compare_omics_source_preservation():
    """Tests that compare_omics does not alter the dataframes it pulls data from."""
    print("Running test_compare_omics_source_preservation...")
    PASS = True

    # Load the source dataframes
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics() # Acetylproteomics and phosphoproteomics have multiple columns for one gene. We use acetylproteomics to make sure compare_omics can grab all those values.

    # Copy the source dataframes so we can make sure later on that compare_omics doesn't alter them.
    prot_copy = prot.copy()
    acet_copy = acet.copy()

    # Call compare_omics on the dataframes, and make sure it doesn't return None.
    compared = en.compare_omics(prot, acet)
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
    """Tests compare_omics with default parameters cols1=None and cols2=None."""
    print("Running test_compare_omics_default_parameters...")
    PASS = True

    # Load the source dataframes
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics() # Acetylproteomics and phosphoproteomics have multiple columns for one gene. We use acetylproteomics to make sure compare_omics can grab all those values.

    # Run the function, make sure it returned properly
    compared = en.compare_omics(prot, acet) 
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

def test_compare_omics_single_gene():
    """Tests compare_omics with single genes for cols1 and cols2."""
    print("Running test_compare_omics_single_gene...")
    PASS = True

    # Load the source dataframes
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics() # Acetylproteomics and phosphoproteomics have multiple columns for one gene. We use acetylproteomics to make sure compare_omics can grab all those values.

    # Run the function, make sure it returned properly
    prot_gene = 'TP53'
    acet_gene = 'A2M'
    compared = en.compare_omics(prot, acet, prot_gene, acet_gene) 
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
    """Tests compare_omics with lists of genes for cols1 and cols2."""
    print("Running test_compare_omics_multiple_genes...")
    PASS = True

    # Load the source dataframes
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics() # Acetylproteomics and phosphoproteomics have multiple columns for one gene. We use acetylproteomics to make sure compare_omics can grab all those values.

    # Run the function, make sure it returned properly
    prot_genes = ['A4GALT', 'TP53', 'ZSCAN30']
    acet_genes = ['AAGAB', 'AACS', 'ZW10', 'ZYX']
    compared = en.compare_omics(prot, acet, prot_genes, acet_genes) 
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
    acet = en.get_acetylproteomics()
    cna = en.get_cna()
    phosg = en.get_phosphoproteomics_gene()
    phoss = en.get_phosphoproteomics_site()
    prot = en.get_proteomics()
    tran = en.get_transcriptomics_linear()
    gene1 = 'TP53'
    gene2 = 'AAGAB'

    # Call compare_omics on the dataframes
    acet_cna = en.compare_omics(acet, cna, gene1, gene2)
    phosg_phoss = en.compare_omics(phosg, phoss, gene1, gene2)
    prot_tran = en.compare_omics(prot, tran, gene1, gene2)

    # Check the return values
    if not check_returned_is_df(acet_cna):
        print("Dataframes compared: acetylproetomics and cna.")
        PASS = False
    if not check_returned_is_df(phosg_phoss):
        print("Dataframes compared: phosphoproteomics_gene and phosphoproteomics_site.")
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
    tran_cir = en.get_transcriptomics_circular() # Although transcriptomics_circular is an omics dataframe, it's of the wrong format to work with compare_omics

    # Test with one valid dataframe and one invalid one
    comp = en.compare_omics(prot, tran_cir)
    if comp is not None:
        print("compare_omics should have returned None when passed the {} dataframe, but instead returned a {}".format(tran_cir.name, type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid dataframe message above was expected.)")

    # Test with two invalid dataframes
    comp = en.compare_omics(clin, tran_cir)
    if comp is not None:
        print("compare_omics should have returned None when passed the {} and {} dataframes, but instead returned a {}".format(clin.name, tran_cir.name, type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid dataframe message above was expected.)")

    # Print whether the test passed
    print_test_result(PASS)

def test_compare_omics_invalid_keys():
    """Test that compare_omics will gracefully handle an invalid key, either alone or in a list of valid keys."""
    print("Running test_compare_omics_single_key_invalid...")
    PASS = True

    # Load dataframes to test with, and set keys to use
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics()
    invalid = 'gobbledegook'
    prot_invalid_list = ['A4GALT', 'TP53', 'ZSCAN30', 'gobbledegook']
    acet_valid = 'AACS' 
    acet_valid_list = ['AAGAB', 'AACS', 'ZW10', 'ZYX']
    acet_invalid_list = ['AAGAB', 'AACS', 'ZW10', 'ZYX', 'gobbledegook']

    # Test one invalid key and one valid key
    comp = en.compare_omics(prot, acet, invalid, acet_valid)
    if comp is not None:
        print("compare_omics should have returned None when passed one invalid key and one valid key, but instead returned a {}".format(type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid key message above was expected.)")

    # Test two invalid keys
    comp = en.compare_omics(prot, acet, invalid, invalid)
    if comp is not None:
        print("compare_omics should have returned None when passed two invalid keys, but instead returned a {}".format(type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid key messages above were expected.)")

    # Test one list of valid keys containing an invalid key, and one list of valid keys
    comp = en.compare_omics(prot, acet, prot_invalid_list, acet_valid_list)
    if comp is not None:
        print("compare_omics should have returned None when passed a list of valid keys containing one invalid key, but instead returned a {}".format(type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid key message above was expected.)")

    # Test two lists of valid keys containing one invalid key each
    comp = en.compare_omics(prot, acet, prot_invalid_list, acet_invalid_list)
    if comp is not None:
        print("compare_omics should have returned None when passed two lists of valid keys containing one invalid key each, but instead returned a {}".format(type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid key messages above were expected.)")

    # Print whether the test passed
    print_test_result(PASS)

def test_compare_omics_invalid_key_types():
    """Tests that compare_omics will gracefully handle a key of an invalid type."""
    print("Running test_compare_omics_invalid_key_types...")
    PASS = True

    # Load our dataframes
    prot = en.get_proteomics()
    acet = en.get_acetylproteomics()

    # Set our keys to use
    prot_valid = 'TP53'
    acet_valid = 'AACS'
    int_key = 100
    prot_valid_list = ['TP53', 'AURKA', 'PIK3CA']
    prot_dict = {0:'TP53', 1:'AURKA', 2:'PIK3CA'} # Create a prep dict for our series we'll use
    prot_series = pd.Series(prot_dict)
    acet_dict = {0:'AAGAB', 1:'AACS', 2:'ZW10'} # Create a prep dict for our series we'll use
    acet_series = pd.Series(acet_dict)

    # Test a key of type int
    comp = en.compare_omics(prot, acet, prot_valid, int_key)
    if comp is not None:
        print("compare_omics should have returned None when passed a key of type int, but instead returned a {}".format(type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid key message above was expected.)")

    # Test two keys of type int
    comp = en.compare_omics(prot, acet, int_key, int_key)
    if comp is not None:
        print("compare_omics should have returned None when passed two keys of type int, but instead returned a {}".format(type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid key messages above were expected.)")

    # Test a key of type pandas.core.series.Series
    comp = en.compare_omics(prot, acet, prot_valid_list, acet_series)
    if comp is not None:
        print("compare_omics should have returned None when passed a pandas.core.series.Series, but instead returned a {}".format(type(comp)))
        PASS = False
    else:
        print("(NOTE: The invalid key message above was expected.)")

    # Test two keys of type pandas.core.series.Series
    comp = en.compare_omics(prot, acet, prot_series, acet_series)
    if comp is not None:
        print("compare_omics should have returned None when passed two pandas.core.series.Series, but instead returned a {}".format(type(comp)))
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
    mut_gene = 'TP53'
    mut_genes = ['TP53', 'PIK3CA']
    acet_gene = 'AAGAB'
    acet_genes = ['AACS', 'ZW10']

    # Copy the source dataframes, to compare at the end
    mut_copy = mut.copy()
    acet_copy = acet.copy()

    # Call append_mutations_to_omics a bunch of times
    en.append_mutations_to_omics(mut_gene, acet)
    en.append_mutations_to_omics(mut_genes, acet)
    en.append_mutations_to_omics(mut_gene, acet, acet_gene)
    en.append_mutations_to_omics(mut_gene, acet, acet_genes)
    en.append_mutations_to_omics(mut_genes, acet, acet_gene)
    en.append_mutations_to_omics(mut_genes, acet, acet_genes)
    en.append_mutations_to_omics(mut_genes, acet, acet_genes, multiple_mutations=True)
    en.append_mutations_to_omics(mut_genes, acet, acet_genes, show_location=False)
    en.append_mutations_to_omics(mut_genes, acet, acet_genes, multiple_mutations=True, show_location=False)

    # Check that the source dataframes weren't changed
    if not mut.equals(mut_copy):
        print("Mutations dataframe was altered by append_mutations_to_omics.")
        PASS = False

    if not acet.equals(acet_copy):
        print("Acetylproteomics dataframe was altered by append_mutations_to_omics.")
        PASS = False

    # Indicate whether the test passed
    print_test_result(PASS)

def test_append_mutations_one_mut_all_omics():
    """Test append_mutations_to_omics with one mutation gene, and the default parameter of None for the omics gene, which should give the entire omics dataframe."""
    print("Running test_append_mutations_single_mut_all_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics_site()
    mut_gene = 'PIK3CA'

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(mut_gene, phos)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{}, with Somatic mutation data for {} gene'.format(phos.name, mut_gene)
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

    mutations = en.get_mutations() # Load the somatic mutations dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_gene):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_three_mut_all_omics():
    """Test append_mutations_to_omics with three mutation genes, and the default parameter of None for the omics gene, which should give the entire omics dataframe."""
    print("Running test_append_mutations_three_mut_all_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics_site()
    mut_genes = ['PIK3CA', 'TP53', 'AURKA']

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(mut_genes, phos)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{}, with Somatic mutation data for {} genes'.format(phos.name, len(mut_genes))
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

    mutations = en.get_mutations() # Load the somatic mutations dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_genes):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_one_mut_one_omics():
    """Test append_mutations_to_omics with one mutation gene and one omics gene."""
    print("Running test_append_mutations_one_mut_one_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics_site()
    phos_gene = 'AAGAB'
    mut_gene = 'PIK3CA'

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(mut_gene, phos, omics_genes=phos_gene)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{} for {}, with Somatic mutation data for {} gene'.format(phos.name, phos_gene, mut_gene)
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

    mutations = en.get_mutations() # Load the somatic mutations dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_gene):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_three_mut_one_omics():
    """Test append_mutations_to_omics with three mutation genes and one omics gene."""
    print("Running test_append_mutations_three_mut_one_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics_site()
    phos_gene = 'AAGAB'
    mut_genes = ['PIK3CA', 'TP53', 'AURKA']

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(mut_genes, phos, omics_genes=phos_gene)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{} for {}, with Somatic mutation data for {} genes'.format(phos.name, phos_gene, len(mut_genes))
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

    mutations = en.get_mutations() # Load the somatic mutations dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_genes):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_one_mut_three_omics():
    """Test append_mutations_to_omics with one mutation gene and three omics genes."""
    print("Running test_append_mutations_one_mut_three_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics_site()
    phos_genes = ['AAGAB', 'AACS', 'ZZZ3']
    mut_gene = 'PIK3CA'

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(mut_gene, phos, omics_genes=phos_genes)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{} for {} genes, with Somatic mutation data for {} gene'.format(phos.name, len(phos_genes), mut_gene)
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

    mutations = en.get_mutations() # Load the somatic mutations dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_gene):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_three_mut_three_omics():
    """Test append_mutations_to_omics with three mutation genes and three omics genes."""
    print("Running test_append_mutations_three_mut_three_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics_site()
    phos_genes = ['AAGAB', 'AACS', 'ZZZ3']
    mut_genes = ['PIK3CA', 'TP53', 'AURKA']

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(mut_genes, phos, omics_genes=phos_genes)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{} for {} genes, with Somatic mutation data for {} genes'.format(phos.name, len(phos_genes), len(mut_genes))
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

    mutations = en.get_mutations() # Load the somatic mutations dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_genes):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_one_mut_all_omics_no_location():
    """Test append_mutations_to_omics with one mutation gene, default of None for omics gene (to select all omics), and no location column."""
    print("Running test_append_mutations_single_mut_all_omics...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics_site()
    mut_gene = 'PIK3CA'

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(mut_gene, phos, show_location=False)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{}, with Somatic mutation data for {} gene'.format(phos.name, mut_gene)
    if not check_df_name(appended, exp_name):
        PASS = False

    # Check dataframe shape
    exp_num_rows = len(phos.index)
    exp_num_cols = len(phos.columns) + 2
    exp_shape = (exp_num_rows, exp_num_cols)
    if not check_df_shape(appended, exp_shape):
        PASS = False

    # Check values in columns
    if not check_appended_columns(phos, appended, phos.columns):
        PASS = False

    mutations = en.get_mutations() # Load the somatic mutations dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_gene, show_location=False):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

def test_append_mutations_three_mut_all_omics_no_location():
    """Test append_mutations_to_omics with three mutation genes, default of None for the omics gene (to select all omics), and no location column."""
    print("Running test_append_mutations_three_mut_all_omics_no_location...")
    PASS = True

    # Load the source dataframe and set our keys
    phos = en.get_phosphoproteomics_site()
    mut_genes = ['PIK3CA', 'TP53', 'AURKA']

    # Run the function, make sure it returned properly
    appended = en.append_mutations_to_omics(mut_genes, phos, show_location=False)
    if not check_returned_is_df(appended):
        PASS = False
        print_test_result(PASS)
        return # Skip remaining steps, since they won't work if it's not a dataframe.

    # Check dataframe name
    exp_name = '{}, with Somatic mutation data for {} genes'.format(phos.name, len(mut_genes))
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

    mutations = en.get_mutations() # Load the somatic mutations dataframe, which the mutation data was drawn from
    if not check_mutation_columns(mutations, appended, mut_genes, show_location=False):
        PASS = False

    # Print whether the test passed
    print_test_result(PASS)

# Multiple mutations, one mutation gene. Use PIK3CA

# Multiple mutations, three mutation genes. Use PIK3CA

# Show location, multiple mutations, one mutation gene. Use PIK3CA

# Show location, multiple mutations, three mutation genes, all omics. Use PIK3CA

# Show location, multiple mutations, three mutation genes, one omics. Use PIK3CA

# Show location, multiple mutations, three mutation genes, multiple omics. Use PIK3CA

# Invalid mutations keys (single and in list)
def test_append_mutations_invalid_key():
    """Test that append_mutations_to_omics gracefully handles invalid mutation gene keys."""
    print("Running test_append_mutations_invalid_key...")
    PASS = True

    # Load our dataframe to test with, and set our keys to use
    acet = en.get_acetylproteomics()
    invalid = 'lorem ipsum'
    invalid_list = ['PIK3CA', 'TAF1', 'GP6', 'lorem ipsum']

    # Test one invalid key
    appended = en.append_mutations_to_omics(invalid, acet)
    if appended is not None:
        print("append_mutations_to_omics should have returned None when passed an invalid key, but instead returned a {}".format(type(appended)))
        PASS = False
    else:
        print("(NOTE: The invalid key message above was expected.)")

    # Test one invalid key in a list of valid keys
    appended = en.append_mutations_to_omics(invalid_list, acet)
    if appended is not None:
        print("append_mutations_to_omics should have returned None when passed a list of valid keys containing one invalid key, but instead returned a {}".format(type(appended)))
        PASS = False
    else:
        print("(NOTE: The invalid key message above was expected.)")

    # Print whether the test passed
    print_test_result(PASS)

# Invalid mutation key types
def test_append_mutations_invalid_key_types():
    """Test that append_mutations_to_omics gracefully handles invalid mutation gene key types."""
    print("Running test_append_mutations_invalid_key_types...")

    # Load our dataframe to test with, and set our keys to use
    prot = en.get_proteomics()
    int_key = 44

    # Test a key of type int

    # Print whether the test passed

def evaluate_special_getters():
    print("Evaluating special getters...")
    results = []
    functions = {}
    results.append(en.get_clinical_cols()); functions[len(results)] = "clinical_cols"
    results.append(en.get_cohort_clinical(["Diabetes","BMI"])); functions[len(results)] = "cohort_meta"
    results.append(en.get_proteomics_quant(["S018","S100"])); functions[len(results)] = "proteomics_quant"
    results.append(en.get_proteomics_cols()); functions[len(results)] = "proteomics_cols"
    results.append(en.get_transcriptomics_cols()); functions[len(results)] = "transcriptomics_cols"
    results.append(en.get_cohort_proteomics(["A1BG","TP53"])); functions[len(results)] = "cohort_proteomics"
    results.append(en.get_cohort_transcriptomics(["A1BG","TP53"])); functions[len(results)] = "cohort_transcriptomics"
    results.append(en.get_cohort_cna(["SASS6","TTTY22"])); functions[len(results)] = "cohort_cna"
    results.append(en.get_cohort_phosphoproteomics(["TP53-S315","AAAS-S541"])); functions[len(results)] = "cohort_phosphoproteomics"
    results.append(en.get_patient_mutations("C3L-00157")); functions[len(results)] = "patient_mutations(Patient_Id)"
    results.append(en.get_patient_mutations("S013")); functions[len(results)] = "patient_mutations(Clinical_Patient_Key)"
    results.append(en.get_phosphosites("TP53")); functions[len(results)] = "phosphosites"
    PASS = True
    for x in range(0,len(results)):
        if results[x] is None:
            print("Error with get",functions[x+1], "function")
            PASS = False
    if PASS:
        print('\tPASS')
    else:
        print("\tFAIL\n")
class Stats:
    def __init__(self):
        pass
    def evaluate(data, trait):
        data_trait = en.compare_clinical(data, trait)
        threshold = .05 / len(data.columns)
        tscutoff = .5
        significantTests = []
        significantGenes = []
        for num in range(1,len(data_trait.columns)):
            gene = data_trait.columns[num]
            oneGene = data_trait[[trait, gene]]
            oneGene = oneGene.dropna(axis=0)
            spearmanrTest = stats.spearmanr(oneGene[trait], oneGene[gene])
            if (abs(spearmanrTest[0]) >= tscutoff) and (spearmanrTest[1] <= threshold):
                significantTests.append(spearmanrTest)
                significantGenes.append(gene)
        if len(significantGenes) > 0:
            return '\tPASS'
        else:
            return "\tFAIL\n"
class Plotter:
    def __init__(self):
        pass
    def plot(data, column1, column2, method):
        if method == "scatterplot":
            plot = sns.relplot(x = column1, y = column2, data = data)
        elif method == "barplot":
            plot = sns.barplot(x = column1, y = column2, data = data)
        elif method == "boxplot":
            plot = sns.boxplot(x = column1, y = column2, data = data)
        else:
            message = method + " not a recognized method"
            print(message)
            return ""
        plt.show()

print("\nRunning tests:\n")

#print("Testing getters...")
#test_get_clinical_filtered()
#test_get_clinical_unfiltered()
#test_get_derived_molecular_filtered()
#test_get_derived_molecular_unfiltered()
#test_get_acetylproteomics_filtered()
#test_get_acetylproteomics_unfiltered()
#test_get_proteomics()
#test_get_transcriptomics_linear()
#test_get_transcriptomics_circular()
#test_get_miRNA()
#test_get_cna()
#test_get_phosphoproteomics_site()
#test_get_phosphoproteomics_gene()
#test_get_phosphosites()
#test_get_mutations_maf()
#test_get_mutations_binary()
#test_get_mutations_unparsed()
#
#print("\nTesting compare and append functions...")
#test_compare_omics_source_preservation()
#test_compare_omics_default_parameters()
#test_compare_omics_single_gene()
#test_compare_omics_multiple_genes()
#test_compare_omics_all_dfs()
#test_compare_omics_invalid_dfs()
#test_compare_omics_invalid_keys()
#test_compare_omics_invalid_key_types()
#test_append_mutations_source_preservation()
#test_append_mutations_one_mut_all_omics()
#test_append_mutations_three_mut_all_omics()
#test_append_mutations_one_mut_one_omics()
#test_append_mutations_three_mut_one_omics()
#test_append_mutations_one_mut_three_omics()
#test_append_mutations_three_mut_three_omics()
#test_append_mutations_one_mut_all_omics_no_location()
#test_append_mutations_three_mut_all_omics_no_location()
test_append_mutations_invalid_key()

#evaluate_special_getters()

# The below tests are not so necessary anymore, now that we have better tests above.

#print("Plotting...")
#Plotter().plot(en.get_proteomics(), "A1BG","PTEN","scatterplot")
#Plotter().plot(en.get_clinical(), "Diabetes","BMI","barplot")
#Plotter().plot(en.get_clinical(), "Diabetes","BMI","boxplot")
#print('\tPASS')

#print("Running statistics...")
#message = Stats().evaluate(en.get_proteomics(), "Tumor_Size_cm")
#print(message)

print("Version:",en.version())
