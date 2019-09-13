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

def check_getter(df, exp_dim, exp_headers, coordinates, values): 
    """Test a dataframe's dimensions and headers, and three test values, then print whether it passed the test.

    Parameters
    df: the dataframe gotten by the getter we are testing
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

# Test functions that get dataframes
def test_get_clinical():
    """Test get_clinical."""

    print('Running test_get_clinical...')

    df = en.get_clinical()
    dimensions = (144, 26)
    headers = ['Patient_ID', 'Proteomics_Tumor_Normal', 'Country', 'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity', 'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site', 'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm', 'Num_full_term_pregnancies']
    test_coord = ((79, 16), (15, 25), (88, 2))
    test_vals = (77.0, '3', 'Poland')

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_derived_molecular():
    """Test get_derived_molecular."""

    print('Running test_get_derived_molecular...')

    df = en.get_derived_molecular()
    dimensions = (144, 125) 
    headers = ['Estrogen_Receptor', 'Estrogen_Receptor_%', 'Progesterone_Receptor', 'Progesterone_Receptor_%', 'MLH1', 'MLH2', 'MSH6', 'PMS2', 'p53', 'Other_IHC_specify', 'Log2_variant_total', 'Log2_SNP_total', 'Log2_INDEL_total', 'Genomics_subtype', 'Mutation_signature_C>A', 'Mutation_signature_C>G', 'Mutation_signature_C>T', 'Mutation_signature_T>C', 'Mutation_signature_T>A', 'Mutation_signature_T>G']
    test_coord = ((3, 4), (30, 117), (80, 52))
    test_vals = ('Intact nuclear expression', 5.459431619, -0.34)

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_experimental_design():
    """Test get_experimental_design."""

    print('Running test_get_experimental_design...')

    df = en.get_experimental_design()
    dimensions = (144, 26) 
    headers = ['Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs', 'Proteomics_Aliquot_ID', 'Proteomics_OCT', 'WXS_normal_sample_type', 'WXS_normal_filename', 'WXS_normal_UUID', 'WXS_tumor_sample_type', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID', 'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality']
    test_coord = ((2, 13), (143, 2), (67, 25))
    test_vals = ('a16b07d8-46c1-4fd9-8204-4f866aacfbec', '130N', 'PASS')

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_acetylproteomics():
    """Test get_acetylproteomics."""

    print('Running test_get_acetylproteomics...')

    df = en.get_acetylproteomics()
    dimensions = (144, 10862)
    headers = ['A2M-K1168', 'A2M-K1176', 'A2M-K135', 'A2M-K145', 'A2M-K516', 'A2M-K664', 'A2M-K682', 'AACS-K391', 'AAGAB-K290', 'AAK1-K201', 'ZSCAN31-K215', 'ZSCAN32-K659', 'ZW10-K634', 'ZYX-K24', 'ZYX-K25', 'ZYX-K265', 'ZYX-K272', 'ZYX-K279', 'ZYX-K533', 'ZZZ3-K117']
    test_coord = ((1, 1), (12, 10861), (90, 5849))
    test_vals = (0.47700000000000004, 0.16, 0.4098)

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_proteomics():
    """Test get_proteomics."""

    print('Running test_get_proteomics...')

    df = en.get_proteomics()
    dimensions = (144, 10999)
    headers = ['A1BG', 'A2M', 'A2ML1', 'A4GALT', 'AAAS', 'AACS', 'AADAT', 'AAED1', 'AAGAB', 'AAK1', 'ZSWIM8', 'ZSWIM9', 'ZW10', 'ZWILCH', 'ZWINT', 'ZXDC', 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3']
    test_coord = ((34, 6003), (99, 9544), (143, 32))
    test_vals = (0.0461, 1.68, 0.904)

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_transcriptomics():
    """Test get_transcriptomics."""

    print('Running test_get_transcriptomics...')

    df = en.get_transcriptomics()
    dimensions = (109, 28057)
    headers = ['A1BG', 'A1BG-AS1', 'A1CF', 'A2M', 'A2M-AS1', 'A2ML1', 'A2MP1', 'A3GALT2', 'A4GALT', 'A4GNT', 'ZWILCH', 'ZWINT', 'ZXDA', 'ZXDB', 'ZXDC', 'ZYG11A', 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3']
    test_coord = ((22, 25483), (108, 23), (101, 17748))
    test_vals = (0.82, 12.0, 6.19)

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_circular_RNA():
    """Test get_circular_RNA."""

    print('Running test_get_circular_RNA...')

    df = en.get_circular_RNA()
    dimensions = (109, 4945)
    headers = ['circ_chr10_100260218_100262063_CWF19L1', 'circ_chr10_100923975_100926019_SLF2', 'circ_chr10_100923978_100926019_SLF2', 'circ_chr10_100937402_100944128_SLF2', 'circ_chr10_100937402_100950753_SLF2', 'circ_chr10_101584602_101586156_POLL', 'circ_chr10_101667886_101676436_FBXW4', 'circ_chr10_101672915_101676436_FBXW4', 'circ_chr10_101792839_101807901_OGA', 'circ_chr10_101792839_101810314_OGA', 'circ_chrX_80288906_80310233_CHMP1B2P', 'circ_chrX_80289664_80310233_CHMP1B2P', 'circ_chrX_80707427_80719656_BRWD3', 'circ_chrX_80791854_80793772_BRWD3', 'circ_chrX_84096194_84164387_RPS6KA6', 'circ_chrX_84134782_84164387_RPS6KA6', 'circ_chrX_85067127_85074391_APOOL', 'circ_chrX_85978767_85981809_CHM', 'circ_chrX_91414904_91418871_PABPC5-AS1', 'circ_chrX_9691579_9693419_TBL1X']
    test_coord = ((108, 1), (30, 4935), (73, 2003))
    test_vals = (9.08, 6.56, 0.0)

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_miRNA():
    """Test get_miRNA."""

    print('Running test_get_miRNA...')

    df = en.get_miRNA()
    dimensions = (99, 2337)
    headers = ['hsa-let-7a-2-3p', 'hsa-let-7a-3p', 'hsa-let-7a-5p', 'hsa-let-7b-3p', 'hsa-let-7b-5p', 'hsa-let-7c-3p', 'hsa-let-7c-5p', 'hsa-let-7d-3p', 'hsa-let-7d-5p', 'hsa-let-7e-3p', 'hsa-miR-9901', 'hsa-miR-9902', 'hsa-miR-9903', 'hsa-miR-9983-3p', 'hsa-miR-9985', 'hsa-miR-9986', 'hsa-miR-99a-3p', 'hsa-miR-99a-5p', 'hsa-miR-99b-3p', 'hsa-miR-99b-5p']
    test_coord = ((5, 0), (98, 1597), (54, 2231))
    test_vals = (1.79, 1.36, 0.26)
    
    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_CNV():
    """Test get_CNV."""

    print('Running test_get_CNV...')

    df = en.get_CNV()
    dimensions = (95, 28057)
    headers = ['A1BG', 'A1BG-AS1', 'A1CF', 'A2M', 'A2M-AS1', 'A2ML1', 'A2MP1', 'A3GALT2', 'A4GALT', 'A4GNT', 'ZWILCH', 'ZWINT', 'ZXDA', 'ZXDB', 'ZXDC', 'ZYG11A', 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3']
    test_coord = ((12, 27865), (60, 8), (94, 15439))
    test_vals = (0.11, 0.01, -0.01)

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_phosphoproteomics():
    """Test get_phosphoproteomics."""

    print('Running test_get_phosphoproteomics...')

    df =  en.get_phosphoproteomics()
    dimensions = (144, 73212)
    headers = ['AAAS-S495', 'AAAS-S541', 'AAAS-Y485', 'AACS-S618', 'AAED1-S12', 'AAGAB-S310', 'AAGAB-S311', 'AAK1-S14', 'AAK1-S18', 'AAK1-S20', 'ZZZ3-S397', 'ZZZ3-S411', 'ZZZ3-S420', 'ZZZ3-S424', 'ZZZ3-S426', 'ZZZ3-S468', 'ZZZ3-S89', 'ZZZ3-T415', 'ZZZ3-T418', 'ZZZ3-Y399']
    test_coord = ((36, 46), (12, 72436), (96, 45361))
    test_vals = (0.579, 0.669, 0.156)

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_phosphoproteomics_gene():
    """Test get_phosphoproteomics_gene."""

    print('Running test_get_phosphoproteomics_gene...')

    df = en.get_phosphoproteomics_gene()
    dimensions = (144, 8466)
    headers = ['AAAS', 'AACS', 'AAED1', 'AAGAB', 'AAK1', 'AAMDC', 'AARS', 'AASDH', 'AATF', 'ABCA1', 'ZSCAN5C', 'ZSWIM3', 'ZSWIM8', 'ZUP1', 'ZW10', 'ZXDA', 'ZXDC', 'ZYX', 'ZZEF1', 'ZZZ3']
    test_coord =  ((2, 7999), (143, 1045), (71, 6543))
    test_vals = (-0.0879, 0.929, 0.153)

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_phosphosites():
    """Test get_phosphosites."""

    print('Running test_get_phosphosites...')

    gene = 'AAK1'
    df = en.get_phosphosites(gene)
    dimensions = (144, 37)
    headers = ['AAK1-S14_phosphoproteomics', 'AAK1-S18_phosphoproteomics', 'AAK1-S20_phosphoproteomics', 'AAK1-S21_phosphoproteomics', 'AAK1-S26_phosphoproteomics', 'AAK1-S618_phosphoproteomics', 'AAK1-S623_phosphoproteomics', 'AAK1-S624_phosphoproteomics', 'AAK1-S637_phosphoproteomics', 'AAK1-S642_phosphoproteomics', 'AAK1-T448_phosphoproteomics', 'AAK1-T606_phosphoproteomics', 'AAK1-T620_phosphoproteomics', 'AAK1-T640_phosphoproteomics', 'AAK1-T653_phosphoproteomics', 'AAK1-T674_phosphoproteomics', 'AAK1-T681_phosphoproteomics', 'AAK1-T694_phosphoproteomics', 'AAK1-T848_phosphoproteomics', 'AAK1-Y687_phosphoproteomics']
    test_coord = ((5, 33), (64, 14), (128, 0))
    test_vals = (0.547, -0.5379999999999999, 0.1395)

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_somatic_mutation():
    """Test get_somatic_mutation."""

    print('Running test_get_somatic_mutation...')

    df = en.get_somatic_mutation()
    dimensions = (52560, 3)
    headers = ['Gene', 'Mutation', 'Location']
    test_coord = ((52000, 2), (12, 0), (34567, 1))
    test_vals = ('p.T162I', 'ADGRA3', 'Missense_Mutation')

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

def test_get_somatic_mutation_binary():
    """Test get_somatic_mutation_binary."""

    print('Running test_get_somatic_mutation_binary...')

    df = en.get_somatic_mutation_binary()
    dimensions = (95, 51559)
    headers = ['A1BG_p.E298K', 'A1BG_p.S181N', 'A1CF_p.F487L', 'A1CF_p.S236Y', 'A2ML1_p.A8V', 'A2ML1_p.G1306D', 'A2ML1_p.L1347F', 'A2ML1_p.L82I', 'A2ML1_p.P712S', 'A2ML1_p.R443Q', 'ZYG11A_p.Q442H', 'ZYG11B_p.H315R', 'ZYG11B_p.R495M', 'ZYG11B_p.R728C', 'ZYX_p.C447Y', 'ZZEF1_p.A2723V', 'ZZEF1_p.D845Y', 'ZZEF1_p.K1251E', 'ZZEF1_p.K2387Sfs*40', 'ZZZ3_p.Y891C']
    test_coord = ((94, 51558), (0, 0), (45, 25436))
    test_vals = (0, 0, 0)

    PASS = check_getter(df, dimensions, headers, test_coord, test_vals)
    print_test_result(PASS)

en = cptac.Endometrial(version="latest")

print("\nRunning tests:\n")
 
print("Testing getters...")
test_get_clinical()
test_get_derived_molecular()
test_get_experimental_design()
test_get_acetylproteomics()
test_get_proteomics()
test_get_transcriptomics()
test_get_circular_RNA()
test_get_miRNA()
test_get_CNV()
test_get_phosphoproteomics()
test_get_phosphoproteomics_gene()
test_get_phosphosites()
test_get_somatic_mutation()
test_get_somatic_mutation_binary()

print("\nTesting compare and join functions...")
# Run tests!

print("Version:", cptac.version())
