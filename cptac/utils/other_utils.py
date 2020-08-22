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

import os
import pandas as pd
import numpy as np
import requests
import webbrowser
import warnings

from cptac.exceptions import DropFromSingleIndexError, InvalidParameterError, MissingFileError, DuplicateColumnHeaderWarning, FileNotUpdatedWarning, FlattenSingleIndexWarning

def get_corum_protein_lists(update=True):
    """Reads file from CORUM and returns a dictionary where the keys are protein complex names, and the values are lists of proteins that are members of those complexes. Data is downloaded from the CORUM website (https://mips.helmholtz-muenchen.de/corum/#). We also provide get_hgnc_protein_lists to get similar data from HGNC. The CORUM data has more specific subgroups than the HGNC data, but the HGNC data is more comprehensive than the CORUM data--it contains proteins that aren't included in the CORUM data.
    Parameters:
    update (bool, optional): Whether to download the latest version of the file from CORUM. Default True. Otherwise, uses a previously downloaded copy (if one exists).

    Returns:
    dict: Keys are complex names; values are lists of proteins in each complex.
    """

    # Set the paths we need
    path_here = os.path.abspath(os.path.dirname(__file__))
    data_files_path = os.path.join(path_here, "data")
    corum_file_path = os.path.join(data_files_path, 'corum_protein_complexes.tsv.zip')

    if update:

        corum_url = "https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip"

        try:
            response = requests.get(corum_url)
            response.raise_for_status() # Raises a requests.HTTPError if the response code was unsuccessful

        except requests.RequestException: # Parent class for all exceptions in the requests module
            warnings.warn("Insufficient internet to update data file. Data from most recent download will be used.", FileNotUpdatedWarning, stacklevel=2)

        else:
            # Check that the data directory exists, create if it doesn't
            if not os.path.isdir(data_files_path):
                os.mkdir(data_files_path)

            # Save the file
            with open(corum_file_path, 'wb') as dest:
                dest.write(response.content)

    # Make sure the file exists
    if not os.path.isfile(corum_file_path):
        raise MissingFileError("CORUM data file has not been downloaded previously, and either you passed False to the 'update' parameter, or there wasn't sufficient internet to update the file (in which case a warning has been issued telling you that). Depending on which of these situations is currently yours, either pass True to the 'update' parameter, or try again when you have a better internet connection.")

    member_proteins = pd.read_csv(corum_file_path, sep='\t')
    member_proteins = member_proteins.loc[member_proteins['Organism'] == 'Human']
    member_proteins = member_proteins.set_index("ComplexName")

    # Select the member proteins column and split the proteins in each complex into values of a list
    member_proteins = member_proteins['subunits(Gene name)'].str.split(';')

    # For complexes with multiple entries (due to different samples), combine the lists
    member_proteins = member_proteins.groupby(member_proteins.index).agg(sum) # Sum will concatenate lists
    member_proteins = member_proteins.apply(set).apply(sorted) # Get rid of duplicates by converting to set. Then go back to list.
    member_proteins = member_proteins.to_dict()

    return member_proteins

def get_hgnc_protein_lists(update=True):
    """Reads file from the HGNC gene family dataset and returns a dictionary where the keys are protein complex names, and the values are lists of proteins that are members of those complexes. Data downloaded from the HGNC BioMart server (https://biomart.genenames.org/). We also provide get_corum_protein_lists to get similar data from CORUM. The HGNC data is more comprehensive than the CORUM data--it contains proteins that aren't included in the CORUM data. Additionally, the CORUM data has more specific subgroups than HGNC, so the HGNC data is easier to query when we just want all proteins associated with a particular structure.

    Parameters:
    update (bool, optional): Whether to download the latest version of the file from HGNC. Default True. Otherwise, uses a previously downloaded copy (if one exists).

    Returns:
    dict: Keys are complex names; values are lists of proteins in each complex.
    """
    # Set the paths we need
    path_here = os.path.abspath(os.path.dirname(__file__))
    data_files_path = os.path.join(path_here, "data")
    hgnc_file_path = os.path.join(data_files_path, 'hgnc_protein_families.tsv')

    if update:

        xml_query = '<!DOCTYPE Query><Query client="biomartclient" processor="TSV" limit="-1" header="1"><Dataset name="hgnc_family_mart" config="family_config"><Attribute name="family__name_103"/><Attribute name="family__gene__symbol_104"/></Dataset></Query>'
        hgnc_biomart_url = "https://biomart.genenames.org/martservice/results"

        try:
            response = requests.post(hgnc_biomart_url, data={"query": xml_query})
            response.raise_for_status() # Raises a requests.HTTPError if the response code was unsuccessful

        except requests.RequestException: # Parent class for all exceptions in the requests module
            warnings.warn("Insufficient internet to update data file. Data from most recent download will be used.", FileNotUpdatedWarning, stacklevel=2)

        else:
            # Check that the data directory exists, create if it doesn't
            if not os.path.isdir(data_files_path):
                os.mkdir(data_files_path)

            # Save the file
            with open(hgnc_file_path, 'wb') as dest:
                dest.write(response.content)

    # Make sure the file exists
    if not os.path.isfile(hgnc_file_path):
        raise MissingFileError("HGNC data file has not been downloaded previously, and either you passed False to the 'update' parameter, or there wasn't sufficient internet to update the file (in which case a warning has been issued telling you that). Depending on which of these situations is currently yours, either pass True to the 'update' parameter, or try again when you have a better internet connection.")

    # Read the file
    member_proteins = pd.read_csv(hgnc_file_path, sep='\t')
    member_proteins = member_proteins.set_index("Family name")
    member_proteins = member_proteins["Approved symbol"]

    # Combine multiple rows per family into one row with a list
    member_proteins = member_proteins.groupby(member_proteins.index).agg(list)
    member_proteins = member_proteins.apply(set).apply(sorted) # Get rid of duplicates by converting to set. Then go back to list.
    member_proteins = member_proteins.to_dict()

    return member_proteins

def search(term):
    """Search for a term in a web browser.

    Parameters:
    term (str): term to be searched

    Returns: None
    """
    url = "https://www.google.com/search?q=" + term
    message = f"Searching for {term} in web browser..."
    print(message, end='\r')
    webbrowser.open(url)
    print(" " * len(message), end='\r') # Erase the message

def reduce_multiindex(df, levels_to_drop=None, flatten=False, sep='_', tuples=False, quiet=False):
    """Drop levels from and/or flatten the column axis of a dataframe with a column multiindex.

    Parameters:
    df (pandas.DataFrame): The dataframe to make the changes to.
    levels_to_drop (str, int, or list or array-like of str or int, optional): Levels, or indices of levels, to drop from the dataframe's column multiindex. These must match the names or indices of actual levels of the multiindex. Must be either all strings, or all ints. Default of None will drop no levels.
    flatten (bool, optional): Whether or not to flatten the multiindex. Default of False will not flatten. Cannot be used if tuples=True.
    sep (str, optional): String to use to separate index levels when flattening. Default is underscore. Only relevant if flatten=True.
    tuples (bool, optional): Whether to return the multiindex as a single-level index of tuples. Cannot be used if flatten=True. Default False.
    quiet (bool, optional): Whether to suppress warnings if duplicate column headers being created when column index levels are dropped, or if you tried to flatten or tuple-ify an index with only one level. Default False.

    Returns:
    pandas.DataFrame: The dataframe, with the desired column index changes made.
    """
    # Parameter check
    if flatten and tuples:
        raise InvalidParameterError("You passed 'True' for both 'flatten' and 'tuples'. This is an invalid combination of arguments. Either pass 'True' to 'flatten' to combine index levels and make a single-level index of strings, or pass 'True' to 'tuples' to return a single-level index of tuples; but just pick one or the other.")

    # Make a copy, so the original dataframe is preserved
    df = df.copy(deep=True)

    if levels_to_drop is not None:
        if df.columns.nlevels < 2:
            raise DropFromSingleIndexError("You attempted to drop level(s) from an index with only one level.")

        if isinstance(levels_to_drop, (str, int)):
            levels_to_drop = [levels_to_drop]
        elif not isinstance(levels_to_drop, (list, pd.Series, pd.Index)):
            raise InvalidParameterError(f"Parameter 'levels_to_drop' is of invalid type {type(levels_to_drop)}. Valid types: str, int, list or array-like of str or int, or NoneType.")

        # Check that they're not trying to drop too many columns
        existing_len = len(df.columns.names)
        to_drop_len = len(levels_to_drop)
        if to_drop_len >= existing_len:
            raise InvalidParameterError(f"You tried to drop too many levels from the dataframe column index. The most levels you can drop is one less than however many exist. {existing_len} levels exist; you tried to drop {to_drop_len}.")

        # Check that the levels they want to drop all exist
        to_drop_set = set(levels_to_drop)
        if all(isinstance(level, int) for level in to_drop_set):
            existing_set_indices = set(range(len(df.columns.names)))
            if not to_drop_set <= existing_set_indices:
                raise InvalidParameterError(f"Some level indices in {levels_to_drop} do not exist in dataframe column index, so they cannot be dropped. Existing column level indices: {list(range(len(df.columns.names)))}")
        else:
            existing_set = set(df.columns.names)
            if not to_drop_set <= existing_set:
                raise InvalidParameterError(f"Some levels in {levels_to_drop} do not exist in dataframe column index, so they cannot be dropped. Existing column levels: {df.columns.names}")

        df.columns = df.columns.droplevel(levels_to_drop)

        num_dups = df.columns.duplicated(keep=False).sum()
        if num_dups > 0 and not quiet:
            warnings.warn(f"Due to dropping the specified levels, dataframe now has {num_dups} duplicated column headers.", DuplicateColumnHeaderWarning, stacklevel=2)

    if flatten:
        if df.columns.nlevels < 2 and not quiet:
            warnings.warn("You tried to flatten a column index that didn't have multiple levels, so we didn't actually change anything.", FlattenSingleIndexWarning, stacklevel=2)
            return df

        tuples = df.columns.to_flat_index() # Converts multiindex to an index of tuples
        no_nan = tuples.map(lambda x: [item for item in x if pd.notnull(item) and item != ""]) # Cut any NaNs and empty strings out of tuples
        joined = no_nan.map(lambda x: sep.join(x)) # Join each tuple
        df.columns = joined
        df.columns.name = "Name" # For consistency
    elif tuples:
        if df.columns.nlevels < 2 and not quiet:
            warnings.warn("You tried to turn a column index into tuples, but it didn't have multiple levels so we didn't actually change anything.", FlattenSingleIndexWarning, stacklevel=2)
            return df

        df.columns = df.columns.to_flat_index()

    return df

"""
Takes a cancer object and find the frequently
mutated genes (in the tumor samples) compared to the cutoff.

@Param cancer_object:
    Cancer dataset object from the cptac module.

@Param cutoff:
    Float. Used as a comparison to determine the status of
    gene mutation frequency.

@Return:
    DataFrame of frequently mutated genes passing the cutoff.
    Columns contain the fractions of total unique mutations,
    missense type mutations, and truncation type mutations per gene.

The Missense_Mut column includes:
    In_Frame_Del, In_Frame_Ins, Missense_Mutation

The Truncation_Mut column includes:
    Frame_Shift_Del, Frame_Shift_Ins, Splice_Site,
    Nonsense_Mutation, Nonstop_Mutation

These columns count multiple mutations of one gene in the
same sample, so fractions in the last two columns may
exceed the Unique_Samples_Mut column which only counts if
the gene was mutated once per sample."""

def get_frequently_mutated(cancer_object, cutoff = 0.1):  
    # Get total tumor count
    clinical_df = cancer_object.get_clinical()
    tumor_status = clinical_df[['Sample_Tumor_Normal']]
    tumor = tumor_status.loc[tumor_status['Sample_Tumor_Normal'] == 'Tumor']
    total_tumor_count = float(len(tumor))
    
    # Get mutations data frame
    somatic_mutations = cancer_object.get_somatic_mutation() 

    # Drop silent mutations for Hnscc, Ovarian, and Ccrcc dataset, and synonymous SNV (i.e. silent) mutations in HNSCC
    if 'Silent' in somatic_mutations['Mutation'].unique():
        somatic_mutations = somatic_mutations.loc[somatic_mutations['Mutation'] != 'Silent']
    if 'RNA' in somatic_mutations['Mutation'].unique():
        somatic_mutations = somatic_mutations.loc[somatic_mutations['Mutation'] != 'RNA'] #ignore RNA in LSCC
    if 'synonymous SNV' in somatic_mutations['Mutation'].unique():
        somatic_mutations = somatic_mutations.loc[somatic_mutations['Mutation'] != 'synonymous SNV']
        
    origin_df = somatic_mutations.reset_index() #prepare to count unique samples
        
    # Create two categories in Mutation column - 'M': Missense, 'T': Truncation
    if cancer_object.get_cancer_type() in ('hnscc') and cancer_object.version() == '0.1':
        dif_mut_names = True
    elif cancer_object.get_cancer_type() in ('colon'):
        dif_mut_names = True
    else: 
        dif_mut_names = False
        
    if dif_mut_names == True:
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
        print(mutations_replaced_M_T['Mutation'].unique())
    
    # Find frequently mutated genes (total fraction > cutoff)
    # Same steps will be repeated for finding the missense and truncation mutation frequencies
    # Step 1 - group by gene and count unique samples
    # Step 2 - format
    # Step 3 - filter using the cutoff and create fraction 
    count_mutations = origin_df.groupby(['Gene']).nunique()
    count_mutations = count_mutations.rename(columns={"Patient_ID": "Unique_Samples_Mut"}) # Step 2 
    count_mutations = count_mutations.drop(['Mutation', 'Location'], axis = 1)
    fraction_mutated = count_mutations.apply(lambda x: x / total_tumor_count) # Step 3 
    fraction_greater_than_cutoff = fraction_mutated.where(lambda x: x > cutoff) #na used when not > cutoff
    filtered_gene_df = fraction_greater_than_cutoff.dropna() # drop genes below cutoff
    
    
    # Create and join Missense column (following similar steps as seen above) *Counts missense once in sample
    miss = mutations_replaced_M_T.loc[mutations_replaced_M_T['Mutation'] == 'M']
    count_miss = miss.groupby(['Gene']).nunique()
    missense_df = count_miss.rename(columns={"Patient_ID": "Missense_Mut"})
    missense_df = missense_df.drop(['Mutation', 'Location'], axis = 1)
    fraction_missense = missense_df.apply(lambda x: x / total_tumor_count)
    freq_mutated_df = filtered_gene_df.join(fraction_missense, how='left').fillna(0)
    
    # Create and join Truncation column (following similar steps as seen above)
    trunc = mutations_replaced_M_T.loc[mutations_replaced_M_T['Mutation'] == 'T']
    count_trunc = trunc.groupby(['Gene']).nunique()
    truncation_df = count_trunc.rename(columns={"Patient_ID": "Truncation_Mut"})
    truncation_df = truncation_df.drop(['Mutation', 'Location'], axis = 1)
    fraction_truncation = truncation_df.apply(lambda x: x / total_tumor_count)
    freq_mutated_df = freq_mutated_df.join(fraction_truncation, how='left').fillna(0)
    
    if gbm == True:
        # Create and join non-coding column (following similar steps as seen above)
        nc = mutations_replaced_M_T.loc[mutations_replaced_M_T['Mutation'] == 'NC']
        count_nc = nc.groupby(['Gene']).nunique()
        nc_df = count_nc.rename(columns={"Patient_ID": "Non-Coding"})
        nc_df = nc_df.drop(['Mutation', 'Location'], axis = 1)
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
