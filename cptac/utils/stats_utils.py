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
import scipy.stats
import statsmodels.stats.multitest
import operator

from cptac.exceptions import InvalidParameterError

'''
@Param df:
    A dataframe containing the label column, and one or more real valued comparison columns.

@Param label_column:
    The name of the label column. This column must be in the dataframe, and must contain exactly 2 unique values.

@Param comparison_columns (default - will use all in dataframe):
    A list of columns on which t-tests will be performed. Each column must be in the dataframe, and must be real valued.
    If no value is specified, by default it will use every column in the dataframe, aside from the specified label column.

@Param alpha (default = .05):
    Significance level. Will be adjusted using parameter correction_method if more than 1 comparison is done.

@Param equal_var (default = True):
    Whether the variances of the two groups are equal. If True, will perform Student's t test. If False, will perform
    Welch's t test.

@Param return_all (default = False):
    Boolean. If true, will return a dataframe containing all comparisons and p-values, regardless of significance.
    If false, will only return significant comparisons and p-values in the dataframe, or None if no significant comparisons.

@Param correction_method (default = 'bonferroni'):
    String. Specifies method of adjustment for multiple testing. See -
    https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
    - for documentation and available methods.

@Param mincount (default=3):
    The minimum number of samples that must have a recorded value, in order for a ttest to be performed. The default is 3,
    which means that there must be at least 3 values recorded for both the mutation type (e.g missense_mutation) and the wildtype.

@Param pval_return_corrected (default=True):
    Returns corrected pvalues if True,
    Returns uncorrected pvalues if False 

@Return:
    A pandas dataframe of column names and corresponding p-values which were determined to be significant in
    the comparison, sorted by significance (smallest p-values at the head). The 2 columns of the dataframe are
    'Comparison' and 'P_Value'.
    Returns None if dataframe was not formatted properly, or if no comparison was significant.

This method takes as a parameter a dataframe. Must be formatted in the following way. 1 column declared as the label column, with
the name of this column passed in as the second parameter. The Label column must contain exactly 2 unique entries,
and every row in the dataframe must have one of these 2 values in this column. The remaining columns will be real
valued columns on which t-tests will be done. A list of real valued columns on which to do t-tests will be passed in
as the third parameter. No t-test will be done on columns not included in this list.

The wrap_ttest method will then compare the two groups, as partitioned by the two values in the Label column, and
perform t-tests for each real valued column in the passed in list, generating a p-value.
The resulting p-values will be corrected for multiple testing, using a specified 'correction_method', and a dataframe with
the significant results will be returned as a dataframe, sorted by p-value.
'''

def wrap_ttest(df, label_column, comparison_columns=None, alpha=.05, equal_var=True, return_all=False, correction_method='bonferroni', mincount=3, pval_return_corrected=True):
    try:
        '''Verify precondition that label column exists and has exactly 2 unique values'''
        label_values = df[label_column].unique()
        if len(label_values) != 2:
            print("Incorrectly Formatted Dataframe! Label column must have exactly 2 unique values.")
            return None

        '''Partition dataframe into two sets, one for each of the two unique values from the label column'''
        partition1 = df.loc[df[label_column] == label_values[0]]
        partition2 = df.loc[df[label_column] == label_values[1]]

        '''If no comparison columns specified, use all columns except the specified labed column'''
        if not comparison_columns:
            comparison_columns = list(df.columns)
            comparison_columns.remove(label_column)

        '''Determine the number of real valued columns on which we will do t-tests'''
        number_of_comparisons = len(comparison_columns)

        '''Store comparisons and p-values in two arrays'''
        comparisons = []
        pvals = []

        '''Loop through each comparison column, perform the t-test, and record the p-val'''

        for column in comparison_columns:
            if len(partition1[column].dropna(axis=0)) <= mincount:
                continue
            elif len(partition2[column].dropna(axis=0)) <= mincount:
                continue
            else:
                stat, pval = scipy.stats.ttest_ind(
                    a=partition1[column].dropna(axis=0),
                    b=partition2[column].dropna(axis=0), 
                    equal_var=equal_var
                )

                comparisons.append(column)
                pvals.append(pval)

        if len(pvals) == 0: # None of the groups had enough members to pass the mincount
            raise InvalidParameterError("No groups had enough members to pass mincount; no tests run.")

        '''Correct for multiple testing to determine if each comparison meets the new cutoff'''
        results = statsmodels.stats.multitest.multipletests(pvals=pvals, alpha=alpha, method=correction_method)
        reject = results[0]

        '''Format results in a pandas dataframe'''
        results_df = pd.DataFrame(columns=['Comparison','P_Value'])

        '''If return all, add all comparisons and p-values to dataframe'''
        if return_all:
            if pval_return_corrected:
                results_df['Comparison'] = comparisons
                results_df['P_Value'] = results[1]

            else:
                results_df['Comparison'] = comparisons
                results_df['P_Value'] = pvals

            '''Else only add significant comparisons'''
        else:
            for i in range(0, len(reject)):
                if reject[i]:
                    if pval_return_corrected:
                        results_df = results_df.append({'Comparison':comparisons[i],'P_Value':results[1][i]}, ignore_index=True)
                    else:
                        results_df = results_df.append({'Comparison':comparisons[i],'P_Value':pvals[i]}, ignore_index=True)


        '''Sort dataframe by ascending p-value'''
        results_df = results_df.sort_values(by='P_Value', ascending=True)
        results_df = results_df.reset_index(drop=True)

        '''If results df is not empty, return it, else return None'''
        if len(results_df) > 0:
            return results_df
        else:
            return None


    except:
        print("Incorrectly Formatted Dataframe!")
        return None


'''
@Param df: Dataframe.Each column is a different gene/ comparison. Rows contains numeric values (such as proteomics) for correlation test
@Param label_column: String. Name of column that will be your x axis and will be compared to all values in df unless otherwise specified.
@Param alpha: significant level
@Param comparison_columns: columns that will be looped through and used as y axis for correlation test.
All other columns beside label column unless specified here.
@Param correction_method: String. Specifies method of adjustment for multiple testing. See -
https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
    - for documentation and available methods.
This function will return a data frame with the columns comparison, the correlation coefficient, and the p value.
'''
def wrap_pearson_corr(df,label_column, alpha=.05,comparison_columns=None,correction_method='bonferroni',return_all = True):


    #df = df.dropna(axis=1, how="all")

    '''If no comparison columns specified, use all columns except the specified labed column'''
    if not comparison_columns:
        comparison_columns = list(df.columns)
        comparison_columns.remove(label_column)
    '''Store comparisons,p-values, correlation in their own array'''
    comparisons = []
    pvals = []
    correlation=[]


    '''Format results in a pandas dataframe'''
    newdf = pd.DataFrame(columns=['Comparison','Correlation','P_value'])
    for gene in comparison_columns:
        #create subset df with interacting gene/ gene (otherwise drop NaN drops everything)
        df_subset = df[[label_column,gene]]
        #do a linear regression to see if it's a meaningful association
        #dropna will remove rows with nan
        df_subset = df_subset.dropna(axis=0, how="any")
        count_row = df_subset.shape[0]
        if count_row > 20:
            x1 = df_subset[[label_column]].values
            y1 = df_subset[[gene]].values
            x1 = x1[:,0]
            y1 = y1[:,0]
            corr, pval = scipy.stats.pearsonr(x1,y1)

            comparisons.append(gene)
            pvals.append(pval)
            correlation.append(corr)


    '''Correct for multiple testing to determine if each comparison meets the new cutoff'''
    results = statsmodels.stats.multitest.multipletests(pvals=pvals, alpha=alpha, method=correction_method)
    reject = results[0]

    if return_all:
        for i in range(0,len(comparisons)):
            newdf = newdf.append({'Comparison': comparisons[i],"Correlation": correlation[i],'P_value': pvals[i]}, ignore_index=True)

    '''Else only add significant comparisons'''
    if (return_all == False):
            for i in range(0, len(reject)):
                if reject[i]:
                    newdf = newdf.append({'Comparison': comparisons[i],"Correlation": correlation[i],'P_value': pvals[i]}, ignore_index=True)

    '''Sort dataframe by ascending p-value'''
    newdf = newdf.sort_values(by='P_value', ascending=True)
    '''If results df is not empty, return it, else return None'''
    return newdf

def permutation_test_means(group1, group2, num_permutations, paired=False):
    """Use permutation testing to calculate a P value for the difference between the means of two groups. You would use this instead of a Student's t-test if your data do not follow a normal distribution. Note that permutation tests are still subject to the assumption of the Student's t-test that if you want to see if the means of the two groups are different, they need to have the same variance.

    Parameters:
    group1 (pandas.Series): The first group of samples. NaNs will be dropped before any analysis. If doing a paired test, samples from the same subject must have the same index value in both group1 and group2.
    group2 (pandas.Series): The second group of samples. NaNs will be dropped before any analysis. If doing a paired test, samples from the same subject must have the same index value in both group1 and group2.
    num_permutations (int): The number of permutations to perform
    paired (bool, optional): Whether to do a paired test. Default is False.

    Returns:
    float: The difference between the means.
    float: The P value for the null hypothesis that the two groups have the same mean.
    """
    # Drop NaN values
    group1 = group1.dropna()
    group2 = group2.dropna()

    extreme_count = 0

    # Create an independent pseudo-random number generator
    generator = np.random.RandomState(0)

    if paired:

        # Calculate the paired differences, and the paired difference in the means
        paired_diffs = group1.combine(group2, operator.sub).dropna().values
        actual_diff = np.mean(paired_diffs)
        abs_actual_diff = abs(actual_diff)

        for i in range(num_permutations):

            # Randomly flip the signs of the differences and recalculate the mean
            random_signs = generator.choice([1, -1], size=paired_diffs.size)
            diffs_signs_perm = random_signs * paired_diffs
            perm_diff = np.mean(diffs_signs_perm)

            # Keep count of how many are as or more extreme than the actual difference
            if abs(perm_diff) >= abs_actual_diff: # We compare the absolute values for a two-tailed test
                extreme_count += 1

    else:
        # Concatenate the series
        both = group1.append(group2)

        # Calculate the actual difference in the means
        actual_diff = np.mean(group1) - np.mean(group2)
        abs_actual_diff = abs(actual_diff)

        # Set some local references to speed up lookups in the loop
        # We tried this optimization for other loops, but this is the only one where it made a difference
        generator_permutation = generator.permutation
        np_mean = np.mean
        both_values = both.values
        group1_size = group1.size

        for i in range(num_permutations):

            # Permute values
            perm_array = generator_permutation(both_values)

            # Split out permuted groups
            perm_group1 = perm_array[:group1_size]
            perm_group2 = perm_array[group1_size:]

            # Calculate the permutation's difference in the means
            perm_diff = np_mean(perm_group1) - np_mean(perm_group2)

            # Keep count of how many are as or more extreme than the actual difference
            if abs(perm_diff) >= abs_actual_diff: # We compare the absolute values for a two-tailed test
                extreme_count += 1

    # Calculate the P value
    P_val = extreme_count / num_permutations # Don't need to multiply by 2 because we compared the absolute values of difference between means.

    return actual_diff, P_val


def permutation_test_corr(data, num_permutations):
    """Use permutation testing to calculate a P value for the linear correlation coefficient between two variables in several samples. You would use this if your distribution didn't follow the Pearson correlation test's assumption of being bivariate normal.

    Parameters:
    data (pandas.DataFrame): A dataframe where the rows are samples, and the columns are the two variables we're testing correlation between.

    Returns:
    float: The linear correlation coefficient for the two variables.
    float: The P value for the null hypothesis that the correlation coefficient is zero.
    """

    # Check the table dimensions
    if data.shape[1] != 2:
        raise InvalidParameterError(f"Expected 2 columns in dataframe. Found {data.shape[1]}.")

    # Drop NaN values
    data = data.dropna()

    # Extract the values
    var1 = data.iloc[:, 0].values
    var2 = data.iloc[:, 1].values

    # Create an independent pseudo-random number generator
    generator = np.random.RandomState(0)

    # Calculate the actual correlation coefficient
    actual_coef = np.corrcoef(var1, var2)[0, 1]

    extreme_count = 0

    for i in range(num_permutations):
        var1_perm = generator.permutation(var1)
        perm_coef = np.corrcoef(var1_perm, var2)[0, 1]

        # Keep count of how many are as or more extreme than our coefficient
        if abs(perm_coef) >= abs(actual_coef): # We compare the absolute values for a two-tailed test
            extreme_count += 1

    # Calculate the P value
    P_val = extreme_count / num_permutations # Don't need to multiply by 2 because we compared the absolute values of coefficients.

    return actual_coef, P_val
