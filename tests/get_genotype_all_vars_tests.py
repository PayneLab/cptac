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
import cptac.utils as ut

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
    
    if not isinstance(returned, pd.DataFrame):
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


def test_genotype_ccrcc_KRAS():
    
    # test when there is no data in the somatic mutations df
    print('Running get_genotype_all_vars...')
    df = k.get_genotype_all_vars('KRAS')
    
    dimensions = (110, 2)
    headers = ['KRAS', 'Mutation']
    
    # get index (int) of patient_ID
    index_1 = df.index.get_loc('C3L-00010') # Test No_Mutation
    index_2 = df.index.get_loc('C3L-01560')
    index_3 = df.index.get_loc('C3N-00646')
    index_4 = df.index.get_loc('C3L-00800') # No del vals (test more No_Mutation) 
    index_5 = df.index.get_loc('C3L-01281')
    index_6 = df.index.get_loc('C3N-00154') 
    index_7 = df.index.get_loc('C3N-00492') # Test Amp
    index_8 = df.index.get_loc('C3L-01287')
    index_9 = df.index.get_loc('C3N-00852')

    # Test No_Mutation
    test_coord_1 = ((index_1, 1), (index_2, 1), (index_3, 1)) # C3N-01515
    test_vals_1 = ('No_Mutation', 'No_Mutation', 'No_Mutation') 
    test_coord_2 = ((index_4, 1),(index_5, 1),(index_6, 1))
    test_vals_2 = ('No_Mutation', 'No_Mutation', 'No_Mutation')
    # Test Amp 
    test_coord_3 = ((index_7, 1), (index_8, 1), (index_9, 1))
    test_vals_3 = ('Amplification', 'Amplification', 'Amplification') 

    test_coord_vals = [(test_coord_1, test_vals_1), (test_coord_2, test_vals_2), 
                       (test_coord_3, test_vals_3)]

    for coord, vals in test_coord_vals:
        PASS = check_getter(df, dimensions, headers, coord, vals)
    
    print_test_result(PASS)
    
    
    
def test_genotype_gbm_KRAS():
    
    # test when there is no data in the somatic mutations df
    print('Running get_genotype_all_vars...')
    df = g.get_genotype_all_vars('KRAS')
    
    dimensions = (98, 2)
    headers = ['KRAS', 'Mutation']
    
    # get index (int) of patient_ID
    index_1 = df.index.get_loc('C3N-03473') # Test No_Mutation
    index_2 = df.index.get_loc('C3N-03183')
    index_3 = df.index.get_loc('C3N-01515')
    index_4 = df.index.get_loc('C3L-01049') # Test Del (only 2)
    index_5 = df.index.get_loc('C3L-02708')
    index_6 = df.index.get_loc('C3N-02256') 
    index_7 = df.index.get_loc('C3N-01816') # Test Amp
    index_8 = df.index.get_loc('C3N-02769')
    index_9 = df.index.get_loc('C3N-02784')

    # Test No_Mutation
    test_coord_1 = ((index_1, 1), (index_2, 1), (index_3, 1)) # C3N-01515
    test_vals_1 = ('No_Mutation', 'No_Mutation', 'No_Mutation') 
    
    # Test Del (only 2 del)
    test_coord_2 = ((index_4, 1),(index_5, 1),(index_6, 1))
    test_vals_2 = ('Deletion', 'Deletion', 'No_Mutation')
    # Test Amp 
    test_coord_3 = ((index_7, 1), (index_8, 1), (index_9, 1))
    test_vals_3 = ('Amplification', 'Amplification', 'Amplification') 

    test_coord_vals = [(test_coord_1, test_vals_1), (test_coord_2, test_vals_2), 
                       (test_coord_3, test_vals_3)]

    for coord, vals in test_coord_vals:
        PASS = check_getter(df, dimensions, headers, coord, vals)
    
    print_test_result(PASS)
    
    
def test_genotype_hnscc_KRAS():
    # test when there is no data in the somatic mutations df
    print('Running get_genotype_all_vars...')
    df = h.get_genotype_all_vars('KRAS')
    
    dimensions = (109, 2)
    headers = ['KRAS', 'Mutation']
    
    # get index (int) of patient_ID
    index_1 = df.index.get_loc('C3L-00999') # Test No_Mutation
    index_2 = df.index.get_loc('C3N-01946')
    index_3 = df.index.get_loc('C3N-03487')
    index_4 = df.index.get_loc('C3N-01337') # Test Del 
    index_5 = df.index.get_loc('C3N-03012')
    index_6 = df.index.get_loc('C3N-03785') 
    index_7 = df.index.get_loc('C3L-04844') # Test Amp
    index_8 = df.index.get_loc('C3L-00987')
    index_9 = df.index.get_loc('C3N-03488')

    # Test No_Mutation
    test_coord_1 = ((index_1, 1), (index_2, 1), (index_3, 1)) # C3N-01515
    test_vals_1 = ('No_Mutation', 'No_Mutation', 'No_Mutation') 
    
    # Test Del 
    test_coord_2 = ((index_4, 1),(index_5, 1),(index_6, 1))
    test_vals_2 = ('Deletion', 'Deletion', 'Deletion')
    # Test Amp 
    test_coord_3 = ((index_7, 1), (index_8, 1), (index_9, 1))
    test_vals_3 = ('Amplification', 'Amplification', 'Amplification') 

    test_coord_vals = [(test_coord_1, test_vals_1), (test_coord_2, test_vals_2), 
                       (test_coord_3, test_vals_3)]

    for coord, vals in test_coord_vals:
        PASS = check_getter(df, dimensions, headers, coord, vals)
    
    print_test_result(PASS)

print("\nRunning tests:\n")
test_genotype_ccrcc_KRAS()
test_genotype_gbm_KRAS()
test_genotype_hnscc_KRAS()

print("Version:", cptac.version())
