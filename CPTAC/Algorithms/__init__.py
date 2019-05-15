import pandas as pd
import numpy as np
import scipy.stats
import re
import sys
import urllib3
import json
import operator
import collections

'''
@Param df:
    A dataframe containing the label column, and one or more real valued comparison columns.
    
@Param label_column:
    The name of the label column. This column must be in the dataframe, and must contain exactly 2 unique values.
    
@Param comparison_columns:
    A list of columns on which t-tests will be performed. Each column must be in the dataframe, and must be real valued.

@Param alpha (default = .05):
    Significance level. Will be adjusted using bonferroni correction if more than 1 comparison is done.
    
@Param verbose (default = False):
    Boolean. If true, will print p-value of every comparison, whether or not it meets significance cutoff.

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
The cutoff for significance will be determined using a bonferroni correction, and the significant columns, 
with their p-values, will be returned as a dataframe, sorted by p-value.
'''

def wrap_ttest(df, label_column, comparison_columns, alpha=.05, verbose=False):
    try:
        '''Verify precondition that label column exists and has exactly 2 unique values'''
        label_values = df[label_column].unique()
        if len(label_values) != 2:
            print("Incorrectly Formatted Dataframe! Label column must have exactly 2 unique values.")
            return None
        
        '''Partition dataframe into two sets, one for each of the two unique values from the label column'''
        partition1 = df.loc[df[label_column] == label_values[0]]
        partition2 = df.loc[df[label_column] == label_values[1]]
        
        '''Determine the number of real valued columns on which we will do t-tests'''
        number_of_comparisons = len(comparison_columns)
        
        '''Use a bonferroni correction to adjust for multiple testing by altering the p-value needed for acceptance'''
        bonferroni_cutoff = alpha/number_of_comparisons
        
        '''Store significant comparisons with their p-values in a dictionary'''
        significant_comparisons = {}
        
        '''Loop through each comparison column, perform the t-test, and determine whether it meets the significance cutoff'''
        for column in comparison_columns:
            stat, pval = scipy.stats.ttest_ind(partition1[column].dropna(axis=0), partition2[column].dropna(axis=0))
            if verbose:
                print(column, ": ", pval)
            if pval <= bonferroni_cutoff:
                significant_comparisons[column] = pval
        
        '''If no comparison met the significance cutoff, notify that no comparison was signficant, and return None'''
        if len(significant_comparisons) == 0:
            print("No significant comparisons.")
            return None
        
            '''If one or more comparison did meet the significance cutoff, sort the dictionary by significance and return it to the caller'''
        else:
            '''Sort dictionary to list smallest p-values first'''
            sorted_significant_comparisons = sorted(significant_comparisons.items(), key=lambda kv: kv[1])
            '''Format as a dataframe and return to caller'''
            sorted_significant_comparisons_df = pd.DataFrame.from_dict(sorted_significant_comparisons)
            sorted_significant_comparisons_df.columns = ['Comparison', 'P_Value']
            return sorted_significant_comparisons_df
    
    except:
        print("Incorrectly Formatted Dataframe!")
        return None

'''
@Param protein:
    The name of the protein that you want to generate a list of interacting proteins for.

@Param number (default=25):
    The number of interacting proteins that you want to get.
    
@Return:
    A list of proteins known by the String api to be interacting partners with the specified protein.
    Returns None if specified protein isn't found in String database, or connection to String api fails.
    
    
This method takes as a parameter the name of a protein. It then accesses the STRING database, through
a call to their public API, and generates a list of proteins known to be interacting partners with the specified
protein. Optional second parameter is number (which by default is 25), which specifies in the API call how many
interacting partners to retrieve from the database. The list of interacting proteins is returned to the caller
as a python list.
'''

def get_interacting_proteins_string(protein, number=25):
    '''Use urllib3 to access the string database api, gather list of interacting proteins'''
    urllib3.disable_warnings()
    string_api_url = "https://string-db.org/api"
    output_format = "json"
    method = "network"

    '''Use the specified gene and homo sapiens species code'''
    my_protein = [protein]
    species = "9606"

    '''Format the api request to collect the appropriate information'''
    request_url = string_api_url + "/" + output_format + "/" + method + "?"
    request_url += "identifiers=%s" % "%0d".join(my_protein)
    request_url += "&" + "species=" + species
    request_url += "&" + "limit=" + str(number)

    '''Send a request to the API, print the response status'''
    try:
        http = urllib3.PoolManager()
        response = http.request('GET',request_url)
        '''Catch exception if it fails while accessing the api'''
    except urllib3.HTTPError as err:
        error_message = err.read()
        print("Error accessing STRING api, " , error_message)
        sys.exit()
    
    '''Get the data from the api response'''
    interacting_proteins = []
    if response.status == 200: 
        '''Get the data from the API's response'''
        data = response.data
        y = json.loads(data)

        '''Make a list of the resulting interacting proteins'''
        for entry in y:
            if entry["preferredName_A"] not in interacting_proteins:
                interacting_proteins.append(entry["preferredName_A"])
            if entry["preferredName_B"] not in interacting_proteins:
                interacting_proteins.append(entry["preferredName_B"])
        
        return interacting_proteins
        
        '''If we didnt get a successful response from the api, notify the caller and return None'''
    else:
        print("\nSpecified gene was not found in String database, double check that you have it correctly!")
        return None


'''
@Param protein:
    The name of the protein that you want to generate a list of interacting proteins for.

@Param number (default=25):
    The number of interacting proteins that you want to get.
    
@Return:
    A list of proteins known by the biogrid api to be interacting partners with the specified protein.
    Returns None if specified protein isn't found in biogrid database, or connection to biogrid api fails.
    
    
This method takes as a parameter the name of a protein. It then accesses the biogrid database, through
a call to their public API, and generates a list of proteins known to be interacting partners with the specified
protein. Optional second parameter is number (which by default is 25), which specifies in the API call how many
interacting partners to retrieve from the database. The list of interacting proteins is returned to the caller
as a python list.
'''
def get_interacting_proteins_biogrid(protein, number=25):
    '''Store interacting proteins in a list'''
    interacting_proteins = []
    urllib3.disable_warnings()
    
    '''Configure url for request'''
    request_url = "https://webservice.thebiogrid.org/interactions/?searchNames=true&geneList=" + protein +"&includeInteractors=true&format=json&taxId=9606&start=0&max=" + str(number) + "&accesskey=0ff59dcf3511928e78aad499688381c9"
    try:
        '''Send request, get response'''
        http = urllib3.PoolManager()
        response = http.request('GET',request_url)
        
        '''If response was successful'''
        if response.status == 200: 
            '''Get the data from the API's response'''
            data = response.data
            y = json.loads(data)
            
            '''Add name of each protein to list of interacting proteins'''
            for entry in y:
                if y[entry]['OFFICIAL_SYMBOL_A'] not in interacting_proteins:
                    interacting_proteins.append(y[entry]['OFFICIAL_SYMBOL_A'])
            
            '''Return this list to caller'''
            return interacting_proteins
        
        else:
            '''If response was not successful, notify caller of error, return None'''
            print("Error accessing api!")
            return None
        
        '''Catch exception, notify caller of error, return None'''
    except Exception as err:
        print("Error accessing api, " , err)
        return None


'''
@Param protein:
    The name of the protein that you want to generate a list of interacting proteins for.

@Param number (default=25):
    The number of interacting proteins that you want to get from both STRING and BioGrid(used by uniprot). This 
    number of proteins will be generated by both String and BioGrid, and the two will be combined. The actual number of 
    proteins in the list returned by this method will be between the number specified and 2 times the number specified, 
    depending on how many of the interacting proteins the two APIs 'agree' on.
    
@Return:
    A list of proteins known by the String and BioGrid APIs to be interacting partners with the specified protein.
    Returns None if specified protein isn't found in either database, or both API calls fail.
    
    
This method takes as a parameter the name of a protein. It then accesses the STRING and BioGrid databases, through
a call to their public API, and generates a list of proteins known to be interacting partners with the specified
protein. Optional second parameter is number (which by default is 25), which specifies in the API call how many
interacting partners to retrieve from the database. The list of interacting proteins is returned to the caller
as a python list.
'''
def get_interacting_proteins(protein, number=25):
    string_list = get_interacting_proteins_string(protein, number)
    biogrid_list = get_interacting_proteins_biogrid(protein, number)
    
    if string_list == None and biogrid_list == None:
        return None
    
    else:
        interacting_proteins = []
        for prot in string_list:
            if prot not in interacting_proteins:
                interacting_proteins.append(prot)
        for prot in biogrid_list:
            if prot not in interacting_proteins:
                interacting_proteins.append(prot)

        return interacting_proteins


def get_frequently_mutated(somatic_df, omics_mutations_df, cutoff=.1, show_percentage=False):  
     
    """take DataFrames of somatic mutations and omics_mutations to determine the percent of 
    mutated genes for all tumors. Frequently mutated genes are greater than the cutoff.

    Parameters:
    somatic_df (pandas.core.frame.DataFrame): Somatic mutations dataframe.
    omics_mutations_df (pandas.core.frame.DataFrame): merged dataframe of any gene and proteomics dataframe
    (used to find total_tumor_patients)
    cutoff (float): used as comparison to determine status of gene mutation frequency

    Returns:
    freq_mutated (list): list of frequently mutated genes passing the cutoff"""
    
    unique_genes = somatic_df['Gene'].unique() # Get series of all mutated genes
    freq_mutated = []
    gene_and_freq_d = {}
    
    # Get total tumor patients
    tumors = omics_mutations_df.loc[omics_mutations_df['Sample_Status'] == 'Tumor']
    total_tumor_patients = len(tumors)
    
    # Find sample (col or index) in gene_mutated (samples represent patients)
    # gene_mutated: Endometrial and Colon samples found in column 0, Ovarian found in index
    gene = 'PTEN'
    gene_mutated = somatic_df.loc[somatic_df['Gene'] == gene]
    ovarian = False
    if gene_mutated.columns[0] == 'Gene':
        ovarian = True
    
    # Find percentage of gene mutation and add frequently mutated genes to dictionary
    if ovarian == True:
        print('ovarian')
        for gene in unique_genes:
            gene_mutated = somatic_df.loc[somatic_df['Gene'] == gene].index
            num_gene_mutated = len(gene_mutated.unique())
            percentage = (num_gene_mutated / total_tumor_patients)
            if percentage > cutoff:
                gene_and_freq_d[gene] = percentage
    
    else:
        for gene in unique_genes:
            gene_mutated = somatic_df.loc[somatic_df['Gene'] == gene].iloc[:, 0]
            gene_mutated.drop_duplicates(keep='first',inplace=True)
            num_gene_mutated = len(gene_mutated)
            percentage = (num_gene_mutated / total_tumor_patients)
            if percentage > cutoff:
                gene_and_freq_d[gene] = percentage

    # Sort dictionary descending order based on percent mutated
    sorted_d = sorted(gene_and_freq_d.items(), key=operator.itemgetter(1), reverse=True)  
    
    # Add frequently mutated gene to list. Option to include percentage.
    for i in range(0,len(sorted_d)):
        certain_tuple = sorted_d[i]
        gene, percent_mutated = certain_tuple
        if show_percentage == True:
            string_gene_percent = gene + ': %' + str('%.2f'%percent_mutated)
            freq_mutated.append(string_gene_percent)
        else:
            freq_mutated.append(gene)  
                   
    return freq_mutated