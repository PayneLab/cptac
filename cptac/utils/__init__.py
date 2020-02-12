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
import re
import sys
import urllib3
import json
import operator
import collections
import os


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
    
@Param return_all (default = False):
    Boolean. If true, will return a dataframe containing all comparisons and p-values, regardless of significance.
    If false, will only return significant comparisons and p-values in the dataframe, or None if no significant comparisons.

@Param correction_method (default = 'bonferroni')
    String. Specifies method of adjustment for multiple testing. See -
    https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
    - for documentation and available methods.
    
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

def wrap_ttest(df, label_column, comparison_columns=None, alpha=.05, return_all=False, correction_method='bonferroni'):
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
            stat, pval = scipy.stats.ttest_ind(partition1[column].dropna(axis=0), partition2[column].dropna(axis=0))
            comparisons.append(column)
            pvals.append(pval)
            
        '''Correct for multiple testing to determine if each comparison meets the new cutoff'''
        results = statsmodels.stats.multitest.multipletests(pvals=pvals, alpha=alpha, method=correction_method)
        reject = results[0]
        
        '''Format results in a pandas dataframe'''
        results_df = pd.DataFrame(columns=['Comparison','P_Value'])

        '''If return all, add all comparisons and p-values to dataframe'''
        if return_all:
            results_df['Comparison'] = comparisons
            results_df['P_Value'] = pvals
        
            '''Else only add significant comparisons'''
        else:
            for i in range(0, len(reject)):
                if reject[i]:
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

        if protein not in interacting_proteins:
            interacting_proteins.append(protein)
        
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
 

'''
@Param protein:
    The name of the protein that you want to generate a list of interacting proteins for.
    
@Return:
    A list of proteins which are interacting partners with the specified protein, according to the bioplex data table.
    Returns None if specified protein isn't found, or no interacting partners are found.
    
This method takes as a parameter the name of a protein. It then accesses the bioplex data table and returns a list of any protein found to be an interacting partner to the given protein.
'''
        
def get_interacting_proteins_bioplex(protein, secondary_interactions=False):
    path_here = os.path.abspath(os.path.dirname(__file__))
    file_name = "BioPlex_interactionList_v4a.tsv"
    file_path = os.path.join(path_here, file_name)
    
    bioplex_interactions = pd.read_csv(file_path, sep='\t')
    
    A_df = bioplex_interactions.loc[bioplex_interactions['SymbolA'] == protein]
    B_df = bioplex_interactions.loc[bioplex_interactions['SymbolB'] == protein]
    
    A_interactions = list(A_df['SymbolB'])
    B_interactions = list(B_df['SymbolA'])
    
    all_interactions = list(set(A_interactions + B_interactions))
    
    if secondary_interactions:
        secondary_interactions_list = []
        for interaction in all_interactions:
            secondary = get_interacting_proteins_bioplex(interaction, False)
            for si in secondary:
                secondary_interactions_list.append(si)
                
        for asi in secondary_interactions_list:
            if asi not in all_interactions:
                all_interactions.append(asi)
    
    if len(all_interactions) > 0:
        return all_interactions
    else:
        return None
    
def get_frequently_mutated(cancer_object, cutoff = 0.1):  
    """
    Takes a cancer object and find the frequently 
    mutated genes (in the tumor samples) compared to the cutoff.
    
    Parameters:
    cancer_object (object): cancer type from cptac module 
    cutoff (float): used as a comparison to determine the 
                    status of gene mutation frequency
    Returns:
    freq_mutated_df (pd.DataFrame): DataFrame of frequently 
        mutated genes passing the cutoff. Columns contain the 
        fractions of total unique mutations, missense type 
        mutations, and truncation type mutations per gene.
    
    The Missense_Mut column includes: 
        In_Frame_Del, In_Frame_Ins, Missense_Mutation
   
    The Truncation_Mut column includes: 
        Frame_Shift_Del, Frame_Shift_Ins, Splice_Site, 
        Nonsense_Mutation, Nonstop_Mutation
        
    These columns count multiple mutations of one gene in the 
    same sample, so fractions in the last two columns may 
    exceed the Unique_Samples_Mut column which only counts if 
    the gene was mutated once per sample.""" 
   
    # Get total tumors/patients
    omics_and_mutations = cancer_object.join_omics_to_mutations(
        mutations_genes = 'TP53', omics_df_name = 'proteomics', omics_genes = 'TP53')
    tumors = omics_and_mutations.Sample_Status

    if isinstance(tumors, pd.core.frame.DataFrame): # This would happen if our proteomics dataframe has a column multiindex, which leads to a joined df with a column multiindex, and causes our selection to be a dataframe instead of a series.
        tumors = tumors.iloc[:, 0]
        tumors.name = "Sample_Status"

    v = tumors.value_counts()
    total_tumors = v['Tumor']
    total_tumor_count = int(total_tumors)
    
    # Get mutations data frame
    somatic_mutations = cancer_object.get_somatic_mutation() 

    # Drop silent mutations for Ovarian and RenalCcrcc dataset, and synonymous SNV (i.e. silent) mutations in HNSCC
    if 'Silent' in somatic_mutations['Mutation'].unique():
        origin_df = somatic_mutations.loc[somatic_mutations['Mutation'] != 'Silent'].reset_index()
    elif 'synonymous SNV' in somatic_mutations['Mutation'].unique():
        origin_df = somatic_mutations.loc[somatic_mutations['Mutation'] != 'synonymous SNV'].reset_index()
    else:
        origin_df = somatic_mutations.reset_index() #prepare to count unique samples
        
    # Create two categories in Mutation column - 'M': Missense, 'T': Truncation
    if cancer_object.get_cancer_type() in ('colon', 'hnscc'):
        missense_truncation_groups = {'frameshift substitution': 'T', 
            'frameshift deletion': 'T', 'frameshift insertion': 'T', 
            'stopgain': 'T', 'stoploss': 'T', 'nonsynonymous SNV': 'M',
            'nonframeshift insertion': 'M','nonframeshift deletion': 'M', 
            'nonframeshift substitution': 'M'}
    else: 
        missense_truncation_groups = {'In_Frame_Del': 'M', 'In_Frame_Ins': 'M',
            'Missense_Mutation': 'M', 'Frame_Shift_Del': 'T','Nonsense_Mutation': 'T', 
            'Splice_Site': 'T', 'Frame_Shift_Ins': 'T','Nonstop_Mutation':'T'}
    
    mutations_replaced_M_T = origin_df.replace(missense_truncation_groups)
    
    # replace non_coding mutations for Gbm
    unique_mutations = len(mutations_replaced_M_T['Mutation'].unique())
    gbm = False
    if cancer_object.get_cancer_type() == 'gbm':
        gbm = True
        non_coding = {'Intron': 'NC', 'RNA': 'NC', "5'Flank": 'NC', "3'Flank": 'NC', 
            "5'UTR": 'NC', "3'UTR": 'NC', 'Splice_Region' : 'NC'}
        mutations_replaced_M_T = mutations_replaced_M_T.replace(non_coding)
        unique_mutations_2 = len(mutations_replaced_M_T['Mutation'].unique())
        
    elif unique_mutations != 2: # Check that all mutation names are catagorized
        print('Warning: New mutation name not classified. Counts will be affected.')
    
    # Find frequently mutated genes (total fraction > cutoff)
    # Same steps will be repeated for finding the missense and truncation mutation frequencies
    # Step 1 - group by gene and count unique samples
    # Step 2 - format
    # Step 3 - filter using the cutoff and create fraction 
    count_mutations = origin_df.groupby(['Gene']).nunique()
    count_mutations = count_mutations.rename(columns={"Patient_ID": "Unique_Samples_Mut"}) # Step 2 
    count_mutations = count_mutations.drop(['Gene', 'Mutation', 'Location'], axis = 1)
    fraction_mutated = count_mutations.apply(lambda x: x / total_tumor_count) # Step 3 
    fraction_greater_than_cutoff = fraction_mutated.where(lambda x: x > cutoff) #na used when not > cutoff
    filtered_gene_df = fraction_greater_than_cutoff.dropna() # drop genes below cutoff
    
    # Create and join Missense column (following similar steps as seen above)
    miss = mutations_replaced_M_T.loc[mutations_replaced_M_T['Mutation'] == 'M']
    count_miss = miss.groupby(['Gene']).nunique()
    missense_df = count_miss.rename(columns={"Patient_ID": "Missense_Mut"})
    missense_df = missense_df.drop(['Gene', 'Mutation', 'Location'], axis = 1)
    fraction_missense = missense_df.apply(lambda x: x / total_tumor_count)
    freq_mutated_df = filtered_gene_df.join(fraction_missense, how='left').fillna(0)
    
    # Create and join Truncation column (following similar steps as seen above)
    trunc = mutations_replaced_M_T.loc[mutations_replaced_M_T['Mutation'] == 'T']
    count_trunc = trunc.groupby(['Gene']).nunique()
    truncation_df = count_trunc.rename(columns={"Patient_ID": "Truncation_Mut"})
    truncation_df = truncation_df.drop(['Gene', 'Mutation', 'Location'], axis = 1)
    fraction_truncation = truncation_df.apply(lambda x: x / total_tumor_count)
    freq_mutated_df = freq_mutated_df.join(fraction_truncation, how='left').fillna(0)
    
    
    if gbm == True:
        # Create and join non-coding column (following similar steps as seen above)
        nc = mutations_replaced_M_T.loc[mutations_replaced_M_T['Mutation'] == 'NC']
        count_nc = nc.groupby(['Gene']).nunique()
        nc_df = count_nc.rename(columns={"Patient_ID": "Non-Coding"})
        nc_df = nc_df.drop(['Gene', 'Mutation', 'Location'], axis = 1)
        fraction_nc = nc_df.apply(lambda x: x / total_tumor_count)
        freq_mutated_df = freq_mutated_df.join(fraction_nc, how='left').fillna(0)
        
    freq_mutated_df = freq_mutated_df.reset_index() #move genes to their own column
    
    return freq_mutated_df

def parse_hotspot(path, mut_df):
    '''
    @Param path:
        (String) The path to the cluster output file that is on your computer after running the Hotspot analysis
    
    @Param mut_df:
        (Dataframe) The dataframe that is obtained by performing the .get_somatic_mutation() function of cptac
        
    @Return:
        There will be four outputs for this function:
        
        vis_hs_df: 
            visualize hotspot dataframe
            
            A small dataframe which will allow quick visualization regarding the number of cancer patients that contain hotspot mutations
        
        bin_hs_df: 
            binary hotspot dataframe
            
            A larger dataframe that contains boolean values for each patient and their relationship with the hotspot(True = patient has a hotspot mutation, False = patient does not have a hotspot mutation)
        
        det_hs_df: 
            detailed hotspot dataframe
            
            A larger dataframe that contains nonbinary values for each patient and their relationship with the hotspot(No = no mutation, Yes = mutation but not in the hotspot, Yes_HS = mutation in the hotspot)
        
        mut_dict:
            mutations dictionary
            
            A dictionary that contains the hotspot gene as the key, and a list of mutations that make up that hotspot
            
    This function will take two parameters (cluster file path and mutations dataframe) and use them to parse the Hotspot3D program output. It creates a cluster dataframe from the Hotspot3D output, and identifies the patients who contain hotspot mutations. The outputs of this function can be used to run further statistical analysis and exploration on the cancer datasets.
    '''
    #Importing the desired cluster file from the specified path on the computer
    cluster_df = pd.read_csv(path, sep = '\t')
    
    #Creating a list of all the identified hotspot clusters
    cluster_list_initial = (cluster_df.Cluster.unique()).tolist()
    cluster_list = list()
    
    #Checking each cluster to make sure that only clusters containing 2 or more mutations are looked at ('clusters' with only 1 mutation are technically just frequently mutated)
    for value in cluster_list_initial:
        length = len(cluster_df[cluster_df['Cluster'] == value])
        if length >= 2:
            cluster_list.append(value)
    
    #Sorting the list numerically
    cluster_list.sort()
    
    #If there are no clusters that have more than one mutation, the function ends and returns the statement below
    if len(cluster_list) == 0:
        print('There are no hotspot clusters that contain more than one mutation.')
        return None
    
    #creating the multiple dictionaries that are used to compile hotspots and corresponding mutations
    gene_dict = {}
    mut_dict = {}
    rev_mut_dict = {}
    hs_count = {}
    
    #This loop contructs a reverse dictionary to be used to classify patients' mutations as well as the mutation dictionary output
    for value in cluster_list:
        gene_dict[value] = cluster_df.loc[cluster_df['Cluster'] == value, 'Gene/Drug'].values[0]
        mut_list = cluster_df[cluster_df['Cluster'] == value]['Mutation/Gene'].values.tolist()
        if str(value).endswith('0'):
            mut_dict[gene_dict[value]] = mut_list
            hs_count[gene_dict[value]] = 0
        else:
            mut_dict[str(gene_dict[value]) + '_' + str(value)[-1]] = mut_list
            hs_count[str(gene_dict[value]) + '_' + str(value)[-1]] = 0
    
    #This loop finalizes the reverse dictionary
    for hs in mut_dict.keys():
        for mutation in mut_dict[hs]:
            rev_mut_dict[mutation] = hs
    
    #The three dataframe outputs are initialized
    vis_hs_df = pd.DataFrame()
    vis_hs_df['hotspot_id'] = mut_dict.keys()
    
    bin_hs_df = pd.DataFrame()
    bin_hs_df['sample_id'] = mut_df.index.unique()
    
    det_hs_df = pd.DataFrame()
    det_hs_df['sample_id'] = mut_df.index.unique()
    
    #This loop populates default values for each patient and hotspot
    for hs in mut_dict.keys():
        bin_hs_df[hs] = False
        det_hs_df[hs] = 'No'
    
    #This loop iterates through each individual mutation and then properly identifies the mutation in the different dataframes
    for row in mut_df.iterrows():
        info = list(row[1])
        gene = info[0]
        location = info[2]
        if str(location)[0] != 'p':
            location = 'p.'+str(location)
        sample_id = row[0]
        
        #This statement checks to see if the mutation is one of the hotspot mutations
        if location in rev_mut_dict.keys():
            hs = rev_mut_dict[location]
            hs_count[hs] += 1
            
            bin_hs_df.loc[bin_hs_df['sample_id'] == sample_id, hs] = True
            det_hs_df.loc[det_hs_df['sample_id'] == sample_id, hs] = 'Yes_HS'
        
        #This statement is used if the mutation is not a hotspot mutation, but if it still on one of the proteins that contains a hotspot
        elif gene in mut_dict.keys():
            det_hs_df.loc[det_hs_df['sample_id'] == sample_id, hs] = 'Yes'
    
    #This loop adds the patient count for each hotspot to the small visualize hotspot dataframe
    for hs in hs_count.keys():
        vis_hs_df.loc[vis_hs_df['hotspot_id'] == hs, 'patients_within'] = hs_count[hs]
        
    bin_hs_df = bin_hs_df.set_index('sample_id')
    det_hs_df = det_hs_df.set_index('sample_id')
    
    #Return of the three dataframes and mutation dictionary
    return(vis_hs_df, bin_hs_df, det_hs_df, mut_dict)
