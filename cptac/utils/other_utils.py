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

import json
import os
import pandas as pd
import requests
import webbrowser
import warnings

from cptac.exceptions import DropFromSingleIndexError, InvalidParameterError, MissingFileError, DuplicateColumnHeaderWarning, FileNotUpdatedWarning, FlattenSingleIndexWarning

def df_to_tree(df):
    df = df.\
    assign(Datatypes=df["Datatypes"].str.split("\ *,\ *", expand=False, regex=True)).\
    explode("Datatypes").\
    reset_index(drop=True)
    # Print our dataframe as a pretty tree structure
    info = {}
    for row in df.set_index(["Cancers", "Sources", "Datatypes"]).index.values:
        if row[0] not in info.keys():
            info[row[0]] = {}
        if row[1] not in info[row[0]].keys():
            info[row[0]][row[1]] = []
        info[row[0]][row[1]].append(row[2])

    df_tree = _tree(info)
    return df_tree

def _tree(nest, prepend=""):
    """Recursively build a formatted string to represent a dictionary"""
    tree_str = ""
    if isinstance(nest, dict):
        for i, (k, v) in enumerate(nest.items()):
            if i == len(nest.keys()) - 1:
                branch = "└"
                newprepend = prepend + "    "
            else:
                branch = "├"
                newprepend = prepend + "│   "
            tree_str += f"{prepend}{branch}── {k}\n"
            tree_str += _tree(nest=v, prepend=newprepend)
    elif isinstance(nest, list):
        for i, v in enumerate(nest):
            if i == len(nest) - 1:
                branch = "└"
            else:
                branch = "├"
            tree_str += f"{prepend}{branch}── {v}\n"
    else:
        raise ValueError(f"Unexpected type '{type(nest)}'")

    return tree_str



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
    
def get_boxnote_text(filepath):
    """Convert a boxnote to plain text.

    Parameters:
    filepath (str): the path to the boxfile

    Returns: the text of the boxfile as a string
    """

    f = open(filepath, encoding="utf8")
    text = json.loads(f.read())["atext"]["text"]
    f.close()
    return text

def validate_levels_to_drop(levels_to_drop, column_names):
    if isinstance(levels_to_drop, (str, int)):
        levels_to_drop = [levels_to_drop]

    if not isinstance(levels_to_drop, (list, pd.Series, pd.Index)):
        raise TypeError(f"Parameter 'levels_to_drop' is of invalid type {type(levels_to_drop)}. Valid types: str, int, list or array-like of str or int, or NoneType.")
    
    levels_to_drop_set = set(levels_to_drop)
    if all(isinstance(level, int) for level in levels_to_drop_set):
        if not levels_to_drop_set.issubset(set(range(len(column_names)))):
            raise ValueError(f"Some level indices in {levels_to_drop} do not exist in dataframe column index.")
    elif not levels_to_drop_set.issubset(set(column_names)):
        raise ValueError(f"Some levels in {levels_to_drop} do not exist in dataframe column index.")

    return levels_to_drop

def drop_levels_from_multiindex(df, levels_to_drop):
    df.columns = df.columns.droplevel(levels_to_drop)
    if df.columns.duplicated(keep=False).sum() > 0:
        warnings("Due to dropping the specified levels, dataframe now has duplicated column headers.", stacklevel=2)

def flatten_multiindex(df, sep='_'):
    if df.columns.nlevels < 2:
        warnings("You tried to flatten a column index that didn't have multiple levels, so we didn't actually change anything.", stacklevel=2)
        return df

    df.columns = df.columns.to_flat_index().map(lambda x: sep.join([item for item in x if pd.notnull(item) and item != ""]))
    df.columns.name = "Name"
    return df

def tupleify_multiindex(df):
    if df.columns.nlevels < 2:
        warnings("You tried to turn a column index into tuples, but it didn't have multiple levels so we didn't actually change anything.", stacklevel=2)
        return df

    df.columns = df.columns.to_flat_index()
    return df

def reduce_multiindex(df, levels_to_drop=None, flatten=False, sep='_', tuples=False):
    if flatten and tuples:
        raise ValueError("Either pass 'True' to 'flatten' or to 'tuples', but not both.")
    
    df = df.copy(deep=True)
    if levels_to_drop is not None:
        levels_to_drop = validate_levels_to_drop(levels_to_drop, df.columns.names)
        drop_levels_from_multiindex(df, levels_to_drop)

    return flatten_multiindex(df, sep) if flatten else tupleify_multiindex(df) if tuples else df

"""
Takes a cancer object and finds the frequently
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


def parse_hotspot(cluster_file_path, mutation_df):
    '''
    Parses the Hotspot3D program output.

    Parameters
    ----------
    cluster_file_path: str
        Path to the cluster output file generated by Hotspot analysis.
    mutation_df: pandas.DataFrame
        Dataframe obtained by performing the .get_somatic_mutation() function of cptac.

    Returns
    -------
    Tuple[pandas.DataFrame, pandas.DataFrame, pandas.DataFrame, Dict]
        Tuple containing the following four items:
            * vis_hs_df: pandas.DataFrame
                Small dataframe for quick visualization of the number of cancer patients that contain hotspot mutations.
            * bin_hs_df: pandas.DataFrame
                Larger dataframe with boolean values representing each patient and their relationship with the hotspot.
            * det_hs_df: pandas.DataFrame
                Larger dataframe with nonbinary values for each patient and their relationship with the hotspot.
            * mut_dict: Dict
                A dictionary with the hotspot gene as the key, and a list of mutations that make up that hotspot as the value.
    '''

    # Load cluster dataframe
    cluster_df = pd.read_csv(cluster_file_path, sep='\t')

    # Identify clusters with 2 or more mutations
    valid_clusters = cluster_df['Cluster'].value_counts() >= 2
    if not valid_clusters.any():
        print('There are no hotspot clusters that contain more than one mutation.')
        return None

    valid_cluster_list = valid_clusters[valid_clusters].index.sort_values()

    # Construct mutation and hotspot dictionaries
    mutation_to_hotspot = {}
    hotspot_to_mutations = {}
    hotspot_counts = {}

    for cluster in valid_cluster_list:
        gene = cluster_df.loc[cluster_df['Cluster'] == cluster, 'Gene/Drug'].values[0]
        mutations = cluster_df[cluster_df['Cluster'] == cluster]['Mutation/Gene'].values.tolist()

        if cluster % 10 == 0:
            hotspot_key = gene
        else:
            hotspot_key = f"{gene}_{cluster % 10}"

        hotspot_to_mutations[hotspot_key] = mutations
        hotspot_counts[hotspot_key] = 0

        for mutation in mutations:
            mutation_to_hotspot[mutation] = hotspot_key

    # Initialize output dataframes
    vis_hs_df = pd.DataFrame(hotspot_to_mutations.keys(), columns=['hotspot_id'])
    bin_hs_df = pd.DataFrame(mutation_df.index.unique(), columns=['sample_id'])
    det_hs_df = bin_hs_df.copy()

    for hotspot_key in hotspot_to_mutations.keys():
        bin_hs_df[hotspot_key] = False
        det_hs_df[hotspot_key] = 'No'

    # Update dataframes based on mutation data
    for row in mutation_df.itertuples():
        sample_id, gene, location = row.Index, row[1], row[3]
        location = 'p.' + location if str(location)[0] != 'p' else location

        if location in mutation_to_hotspot:
            hotspot = mutation_to_hotspot[location]
            hotspot_counts[hotspot] += 1

            bin_hs_df.loc[bin_hs_df['sample_id'] == sample_id, hotspot] = True
            det_hs_df.loc[det_hs_df['sample_id'] == sample_id, hotspot] = 'Yes_HS'
        elif gene in hotspot_to_mutations:
            det_hs_df.loc[det_hs_df['sample_id'] == sample_id, hotspot] = 'Yes'

    vis_hs_df['patients_within'] = hotspot_counts.values()

    bin_hs_df.set_index('sample_id', inplace=True)
    det_hs_df.set_index('sample_id', inplace=True)

    return vis_hs_df, bin_hs_df, det_hs_df, hotspot_to_mutations
