# Tests for get_frequently_mutated from the utilities
# NOTE: change cptac.algorithms to cptac.utilities when the switch is made

import pandas as pd
import cptac

def test_get_frequently_mutated_en_default_cutoff():
    en = cptac.Endometrial()
    print('Running get_frequently_mutated...')
    df = cptac.algorithms.get_frequently_mutated(en)
    
    name = "frequently_mutated"
    dimensions = (232, 4)
    headers = ['Gene', 'Unique_Samples_Mut', 'Missense_Mut', 'Truncation_Mut']
    # test gene names
    test_coord_names = ((53, 0), (32, 0), (227, 0))
    test_vals_names = ('CTCF', 'CCDC168', 'ZNF536')
    
    total_tumors = 95
    # test missense and trucation don't equal the unique_sample_mutated 
    #(miss and trunc in same sample)
    test_coord_CTCF = ((53, 1), (53, 2), (53, 3)) 
    test_vals_CTCF = (27/total_tumors, 9/total_tumors, 23/total_tumors) 
    # testmissense and trucation values are equal
    test_coord_CCDC168 = ((32, 1),(32, 2),(32, 3))
    test_vals_CCDC168 = (16/total_tumors, 11/total_tumors, 11/total_tumors)
    # test no truncation type mutatations
    test_coord_ZNF536 = ((227, 1),(227, 2),(227, 3))
    test_vals_ZNF536 = (12/total_tumors, 12/total_tumors, 0/total_tumors)
    # test close to cutoff
    test_coord_DICER1 = ((61, 1),(61, 2),(61, 3))
    test_vals_DICER1 = (10/total_tumors, 10/total_tumors, 1/total_tumors)
    # common test
    test_coord_TP53 = ((205, 1),(205, 2),(205, 3))
    test_vals_TP53 = (21/total_tumors, 15/total_tumors, 7/total_tumors)

    test_coord_vals = [(test_coord_names, test_vals_names), (test_coord_CTCF, test_vals_CTCF), 
                       (test_coord_CCDC168, test_vals_CCDC168), (test_coord_ZNF536, test_vals_ZNF536),
                       (test_coord_DICER1, test_vals_DICER1), (test_coord_TP53, test_vals_TP53)]

    for coord, vals in test_coord_vals:
        PASS = check_getter(df, dimensions, headers, coord, vals)
    
    print_test_result(PASS)
    
    
def test_get_frequently_mutated_en_cutoff_20_cutoff():
    en = cptac.Endometrial()
    print('Running get_frequently_mutated...')
    df = cptac.algorithms.get_frequently_mutated(en, cutoff=0.2)
    
    dimensions = (10, 4)
    headers = ['Gene', 'Unique_Samples_Mut', 'Missense_Mut', 'Truncation_Mut']
    # test gene names
    test_coord_names = ((0, 0), (2, 0), (8, 0))
    test_vals_names = ('ARID1A', 'CTNNB1', 'TP53')
    
    total_tumors = 95
    # test missense and trucation don't equal the unique_samples_mutated 
    #(miss and trunc in same sample and counted in each category)
    test_coord_ARID1A = ((0, 1), (0, 2), (0, 3)) 
    test_vals_ARID1A = (43/total_tumors, 13/total_tumors, 38/total_tumors) 
    # test no truncation type mutatations
    test_coord_CTNNB1 = ((2, 1),(2, 2),(2, 3))
    test_vals_CTNNB1 = (29/total_tumors, 29/total_tumors, 0/total_tumors)
    # test close to the cutoff 
    test_coord_ZFHX3 = ((9, 1), (9, 2), (9, 3))
    test_vals_ZFHX3 = (21/total_tumors , 8/total_tumors , 16/total_tumors)
    # test miss and trunc almost equal
    test_coord_KMT2B = ((3, 1), (3, 2), (3, 3))
    test_vals_KMT2B = (23/total_tumors , 11/total_tumors , 12/total_tumors)
    # common test
    test_coord_TP53 = ((8, 1),(8, 2),(8, 3))
    test_vals_TP53 = (21/total_tumors, 15/total_tumors, 7/total_tumors)
    
    test_coord_vals = [(test_coord_names, test_vals_names), (test_coord_ARID1A, test_vals_ARID1A),
                        (test_coord_CTNNB1, test_vals_CTNNB1), (test_coord_ZFHX3, test_vals_ZFHX3),
                        (test_coord_KMT2B, test_vals_KMT2B), (test_coord_TP53, test_vals_TP53)]

    for coord, vals in test_coord_vals:
        PASS = check_getter(df, dimensions, headers, coord, vals)
        
    print_test_result(PASS)
    
    
    
def test_get_frequently_mutated_co_default_cutoff():
    co = cptac.Colon()
    print('Running get_frequently_mutated...')
    df = cptac.algorithms.get_frequently_mutated(co)
    
    dimensions = (612, 4)
    headers = ['Gene', 'Unique_Samples_Mut', 'Missense_Mut', 'Truncation_Mut']
    # test gene names
    test_coord_names = ((90, 0), (284, 0), (499, 0))
    test_vals_names = ('CASP5', 'KRAS', 'SPINK5')
    
    total_tumors = 97
    # test when there are no missense type mutatations
    test_coord_CASP5 = ((90, 1), (90, 2), (90, 3))
    test_vals_CASP5 = (19/total_tumors, 0/total_tumors, 19/total_tumors) 
    # test when there are no truncation type mutatations
    test_coord_KRAS = ((284, 1),(284, 2),(284, 3))
    test_vals_KRAS = (35/total_tumors, 35/total_tumors, 0/total_tumors)
    # test when missense and trucation don't add up to equal the fraction mutated
    #(miss and trunc in same sample)
    test_coord_ANK2 = ((34, 1),(34, 2),(34, 3)) 
    test_vals_ANK2 = (15/total_tumors, 13/total_tumors, 4/total_tumors) 
    # test when miss and trunc count are the same
    test_coord_ATM = ((56, 1),(56, 2),(56, 3)) 
    test_vals_ATM = (10/total_tumors, 7/total_tumors, 7/total_tumors) 
    # test close to the cutoff
    test_coord_SPINK5 = ((499, 1),(499, 2),(499,3))
    test_vals_SPINK5 = (10/total_tumors, 5/total_tumors, 7/total_tumors)
    # common test
    test_coord_TP53 = ((554, 1),(554, 2),(554, 3))
    test_vals_TP53 = (56/total_tumors, 38/total_tumors, 21/total_tumors)

    test_coord_vals = [(test_coord_names, test_vals_names), (test_coord_CASP5, test_vals_CASP5),
                        (test_coord_KRAS, test_vals_KRAS), (test_coord_ANK2, test_vals_ANK2),
                        (test_coord_ATM, test_vals_ATM), (test_coord_SPINK5, test_vals_SPINK5), 
                        (test_coord_TP53, test_vals_TP53)]

    for coord, vals in test_coord_vals:
        PASS = check_getter(df, dimensions, headers, coord, vals)
    
    print_test_result(PASS)
    
    
def test_get_frequently_mutated_co_15_cutoff():
    co = cptac.Colon()
    print('Running get_frequently_mutated...')
    df = cptac.algorithms.get_frequently_mutated(co, 0.15)
    
    dimensions = (138, 4)
    headers = ['Gene', 'Unique_Samples_Mut', 'Missense_Mut', 'Truncation_Mut']
    # test gene names
    test_coord_names = ((15, 0), (66, 0), (102, 0))
    test_vals_names = ('CASP5', 'KRAS', 'RYR2')
    
    total_tumors = 97
    # test no missense type mutatations 
    test_coord_CASP5 = ((15, 1), (15, 2), (15, 3))
    test_vals_CASP5 = (19/total_tumors, 0/total_tumors, 19/total_tumors) 
    # test no truncation type mutatations
    test_coord_KRAS = ((66, 1),(66, 2),(66, 3))
    test_vals_KRAS = (35/total_tumors, 35/total_tumors, 0/total_tumors)
    # test missense and truncation equal fraction mutated
    test_coord_PIK3CA = ((92, 1),(92, 2),(92, 3))
    test_vals_PIK3CA = (24/total_tumors, 23/total_tumors, 1/total_tumors)
    # test missense and trucation don't equal unique_samples_mutated (miss and trunc in same sample)
    test_coord_RYR2 = ((102, 1),(102, 2),(102, 3))
    test_vals_RYR2 = (21/total_tumors, 19/total_tumors, 7/total_tumors)
    # test close to the cutoff
    test_coord_ANK2 = ((6, 1),(6, 2),(6, 3)) 
    test_vals_ANK2 = (15/total_tumors, 13/total_tumors, 4/total_tumors) 
    # common test
    test_coord_TP53 = ((123, 1),(123, 2),(123, 3))
    test_vals_TP53 = (56/total_tumors, 38/total_tumors, 21/total_tumors)

    test_coord_vals = [(test_coord_names, test_vals_names), (test_coord_CASP5, test_vals_CASP5), 
                       (test_coord_KRAS, test_vals_KRAS), (test_coord_PIK3CA, test_vals_PIK3CA), 
                       (test_coord_RYR2, test_vals_RYR2), (test_coord_ANK2, test_vals_ANK2), 
                       (test_coord_TP53, test_vals_TP53)]

    for coord, vals in test_coord_vals:
        PASS = check_getter(df, dimensions, headers, coord, vals)
    
    print_test_result(PASS)
    
    
def test_get_frequently_mutated_ov_default_cutoff():
    ov = cptac.Ovarian()
    print('Running get_frequently_mutated...')
    df = cptac.algorithms.get_frequently_mutated(ov)
    
    dimensions = (16, 4)
    headers = ['Gene', 'Unique_Samples_Mut', 'Missense_Mut', 'Truncation_Mut']
    # test genes names
    test_coord_names = ((15, 0), (13, 0), (2, 0))
    test_vals_names = ('WDFY4', 'TP53', 'MT-CO1')
    
    total_tumors = 83
    #test missense and trucation not equal to unique_samples_mutated 
    #(miss and trunc in same sample)
    test_coord_WDFY4 = ((15, 1), (15, 2), (15, 3)) 
    test_vals_WDFY4 = (10/total_tumors, 8/total_tumors, 3/total_tumors) 
    # test miss and trunc equal to unique_samples_mutated
    test_coord_MUC4 = ((8, 1),(8, 2),(8, 3))
    test_vals_MUC4 = (27/total_tumors, 26/total_tumors, 1/total_tumors)
    # test no truncation mutations
    test_coord_MTCO1 = ((2, 1),(2, 2),(2, 3))
    test_vals_MTCO1 = (10/total_tumors, 10/total_tumors, 0/total_tumors)
    # test close to cutoff
    test_coord_FSIP2 = ((1, 1),(1, 2),(1, 3))
    test_vals_FSIP2 = (9/total_tumors, 8/total_tumors, 2/total_tumors)
    # common test and highest count
    test_coord_TP53 = ((13, 1),(13, 2),(13, 3))
    test_vals_TP53 = (77/total_tumors, 50/total_tumors, 27/total_tumors)
    
    
    test_coord_vals = [(test_coord_names, test_vals_names), (test_coord_WDFY4, test_vals_WDFY4),
                       (test_coord_MUC4, test_vals_MUC4), (test_coord_MTCO1, test_vals_MTCO1), 
                       (test_coord_FSIP2, test_vals_FSIP2), (test_coord_TP53, test_vals_TP53)]

    for coord, vals in test_coord_vals:
        PASS = check_getter(df, dimensions, headers, coord, vals)
    
    print_test_result(PASS)
    
    
def test_get_frequently_mutated_ov_05_cutoff():
    ov = cptac.Ovarian()
    print('Running get_frequently_mutated...')
    df = cptac.algorithms.get_frequently_mutated(ov, 0.05)
    
    dimensions = (142, 4)
    headers = ['Gene', 'Unique_Samples_Mut', 'Missense_Mut', 'Truncation_Mut']
    # test genes names
    test_coord_names = ((133, 0), (127, 0), (141, 0))
    test_vals_names = ('WDFY4', 'TP53', 'ZNF865')
    
    total_tumors = 83
    #test missense and trucation not equal to unique_samples_mutated 
    #(miss and trunc in same sample)
    test_coord_WDFY4 = ((133, 1), (133, 2), (133, 3)) 
    test_vals_WDFY4 = (10/total_tumors, 8/total_tumors, 3/total_tumors) 
    # test miss and trunc almost equal
    test_coord_CDK12 = ((11, 1),(11, 2),(11, 3))
    test_vals_CDK12 = (6/total_tumors, 4/total_tumors, 3/total_tumors)
    # test no truncation mutations
    test_coord_ZNF865 = ((141, 1),(141, 2),(141, 3))
    test_vals_ZNF865 = (5/total_tumors, 5/total_tumors, 0/total_tumors)
    # test close to cutoff 
    test_coord_SYNE1 = ((122, 1),(122, 2),(122, 3))
    test_vals_SYNE1 = (5/total_tumors, 5/total_tumors, 1/total_tumors)
    # common test and highest count
    test_coord_TP53 = ((127, 1),(127, 2),(127, 3))
    test_vals_TP53 = (77/total_tumors, 50/total_tumors, 27/total_tumors)
    
    #CHECK silent mut not counted
    
    test_coord_vals = [(test_coord_names, test_vals_names), (test_coord_WDFY4, test_vals_WDFY4),
                        (test_coord_CDK12, test_vals_CDK12), (test_coord_ZNF865, test_vals_ZNF865),  
                        (test_coord_SYNE1, test_vals_SYNE1), (test_coord_TP53, test_vals_TP53)]

    for coord, vals in test_coord_vals:
        PASS = check_getter(df, dimensions, headers, coord, vals)
    
    print_test_result(PASS)
    
def test_get_frequently_mutated_renal_default_cutoff():
    rc = cptac.RenalCcrcc()
    print('Running get_frequently_mutated...')
    df = cptac.algorithms.get_frequently_mutated(rc)
   
    dimensions = (6, 4)
    headers = ['Gene', 'Unique_Samples_Mut', 'Missense_Mut', 'Truncation_Mut']
    # test genes names
    test_coord_names = ((0, 0), (2, 0), (4, 0))
    test_vals_names = ('BAP1', 'PBRM1', 'TTN')
    
    total_tumors = 110
    # test miss and trunc equal to unique_samples_mutated
    test_coord_BAP1 = ((0, 1), (0, 2), (0, 3)) 
    test_vals_BAP1 = (17/total_tumors, 7/total_tumors, 10/total_tumors) 
    # test high truncation, low missense count
    test_coord_PBRM1 = ((2, 1),(2, 2),(2, 3))
    test_vals_PBRM1 = (44/total_tumors, 8/total_tumors, 37/total_tumors)
    # check that silent mutations are not counted (TTN has many silent mutations)
    # and missense and trucation not equal to unique_samples_mutated 
    test_coord_TTN = ((4, 1),(4, 2),(4, 3))
    test_vals_TTN = (13/total_tumors, 10/total_tumors, 4/total_tumors)
    # test close to cutoff
    test_coord_SETD2 = ((3, 1), (3, 2), (3, 3)) 
    test_vals_SETD2 = (15/total_tumors, 2/total_tumors, 13/total_tumors)               
    # common test and highest count
    test_coord_VHL = ((5, 1),(5, 2),(5, 3))
    test_vals_VHL = (82/total_tumors, 33/total_tumors, 49/total_tumors)
    
    
    test_coord_vals = [(test_coord_names, test_vals_names), (test_coord_BAP1, test_vals_BAP1),
                        (test_coord_PBRM1, test_vals_PBRM1), (test_coord_TTN, test_vals_TTN),
                        (test_coord_SETD2, test_vals_SETD2), (test_coord_VHL, test_vals_VHL)]

    for coord, vals in test_coord_vals:
        PASS = check_getter(df, dimensions, headers, coord, vals)
    
    print_test_result(PASS)
    
    
def test_get_frequently_mutated_renal_01_cutoff():
    rc = cptac.RenalCcrcc()
    print('Running get_frequently_mutated...')
    df = cptac.algorithms.get_frequently_mutated(rc, cutoff=0.01)
    
    dimensions = (1106, 4)
    headers = ['Gene', 'Unique_Samples_Mut', 'Missense_Mut', 'Truncation_Mut']
    # test genes names
    test_coord_names = ((11, 0), (992, 0), (1080, 0))
    test_vals_names = ('ABCC3', 'TTN', 'ZNF532')
    
    total_tumors = 110
    # test no missense 
    test_coord_ABCC3 = ((11, 1),(11, 2),(11, 3))
    test_vals_ABCC3 = (2/total_tumors, 0/total_tumors, 2/total_tumors)
    # test no truncation and close to cutoff
    test_coord_ZNF532 = ((1080, 1), (1080, 2), (1080, 3)) 
    test_vals_ZNF532 = (2/total_tumors, 2/total_tumors, 0/total_tumors)
    # test miss and trunc equal to unique_samples_mutated
    test_coord_NAV3 = ((611, 1), (611, 2), (611, 3)) 
    test_vals_NAV3 = (7/total_tumors, 5/total_tumors, 2/total_tumors) 
    # check that silent mutations are not counted (TTN has many silent mutations)
    # and missense and trucation not equal to unique_samples_mutated 
    test_coord_TTN = ((992, 1),(992, 2),(992, 3))
    test_vals_TTN = (13/total_tumors, 10/total_tumors, 4/total_tumors)
    # common test and highest count
    test_coord_VHL = ((1019, 1),(1019, 2),(1019, 3))
    test_vals_VHL = (82/total_tumors, 33/total_tumors, 49/total_tumors)
    
    
    test_coord_vals = [(test_coord_names, test_vals_names), (test_coord_ABCC3, test_vals_ABCC3),
                        (test_coord_ZNF532, test_vals_ZNF532), (test_coord_NAV3, test_vals_NAV3),
                        (test_coord_TTN, test_vals_TTN), (test_coord_VHL, test_vals_VHL)]

    for coord, val in test_coord_vals:
        PASS = check_getter(df, dimensions, headers, coord, val)
    
    print_test_result(PASS)
    
    
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
    

# Run tests
print("\nTesting get_frequently_mutated from algorithms...")
import cptac as cptac
test_get_frequently_mutated_en_default_cutoff()
test_get_frequently_mutated_co_default_cutoff()
test_get_frequently_mutated_ov_default_cutoff()
test_get_frequently_mutated_renal_default_cutoff()


test_get_frequently_mutated_en_cutoff_20_cutoff()
test_get_frequently_mutated_co_cutoff_15_cutoff()
test_get_frequently_mutated_ov_05_cutoff()
test_get_frequently_mutated_renal_01_cutoff()