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
import warnings
from cptac.exceptions import CptacDevError, ReindexMapError, FailedReindexWarning
from contextlib import contextmanager
import sys, os
import mygene


@contextmanager
def suppress_stdout():
    """Temporarily stops output from being printed. Use before a function when you don't want its output to
    be printed.
    Example: with suppress_stdout():
                 function()
    """
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


def rename_duplicate_labels(df, label_type='columns'):
    """Returns a df with unique labels for columns or indices
    Parameters:
    df (pandas.DataFrame): The df containing duplicate labels.
    label_type (str): use 'columns' to make unique column labels and 'index' to make unique indices.

    Returns:
    pandas.DataFrame: df with with unique labels. ".i" is appended to every duplicate, increasing incrementally.
    """
    if label_type == 'columns':
        labs = pd.Series(df.columns[:])

        for dup in labs[labs.duplicated()].unique():
            labs[labs[labs == dup].index.values.tolist()] = [dup + '.' + str(i) if i != 0 else dup for i in range(sum(labs == dup))]

        # rename the columns with the labs list
        df.columns=labs

    elif label_type == 'index':
        labs = pd.Series(df.index[:])

        for dup in labs[labs.duplicated()].unique():
            labs[labs[labs == dup].index.values.tolist()] = [dup + '.' + str(i) if i != 0 else dup for i in range(sum(labs == dup))]

        # rename the columns with the labs list
        df.index=labs
    return df

def average_replicates(df, id_list = [], normal_identifier = '.N', common = '\.', to_drop = '\.\d$'):
    """Returns a df with one row for each patient_ID (all replicates for a patient are averaged)
    Parameters:
    df (pandas.DataFrame): The df containing replicates (duplicate entries for the same tissue_type).
    id_list: list of IDs with replicates (use the ID format that is common between replicates so the
            list can be used to slice out all replicates for that ID). Make sure the IDs
            in the list include the symbol that distinguishes normal samples ('.N' or '-N').
    common: regex string that is common between replicates (identifies duplicate entries)
    to_drop: regex string to drop to find each patient_ID that has replicates (used to slice out all replicates)

    Returns:
    pandas.DataFrame: df with with replicate rows averaged and one row for each patient_ID.
    """
    # If no list of replicate IDs is given, make list from common regex
    if len(id_list) == 0:
        replicate_df = df[df.index.str.contains(common)]
        patient_ids = pd.Series(replicate_df.index) # create series of replicate IDs to prep removing appended ".i"
        ids = patient_ids.replace(to_drop, '', regex=True)
        id_list = list(set(ids)) #id_list contains only patient_IDs of replicates (without #s)

    new_df = df.copy()
    for patient_ID in id_list:
        # Can slice only normals with patient_ID because of normal identifier (won't slice out tumors)
        if normal_identifier in patient_ID:
            id_df = df[df.index.str.contains(patient_ID, regex = True)] # slice out replicates for a single patient
        # If tumor, need to slice out normals
        else:
            if normal_identifier == '.N':
                norm_regex = '\.N' # prep for regex use
            else:
                norm_regex = normal_identifier
            id_df = df[df.index.str.contains(patient_ID, regex = True) & \
                       ~ df.index.str.contains(norm_regex, regex = True)] # don't include normals
        #print(id_df.index.to_list())
        vals = list(id_df.mean(axis=0))
        new_df = new_df.drop(id_df.index.to_list(), axis = 'index') # drop unaveraged rows
        new_df.loc[patient_ID] = vals # add averaged row

    return new_df

def get_reindex_map(series):
    """Generate a reindexing map from a series where the index is the new indices, and the values are the old indices.

    Parameters:
    series (pandas Series): The series to generate the reindex map from.

    Returns:
    pandas Series: The reindexing map, with the old inidicies as the index, and the new indices as the values.
    """
    old_index_name = series.name
    new_index_name = series.index.name
    if new_index_name is None:
        new_index_name = "index" # This is the default name pandas gives unnamed indices when they're made columns

    # Check that the mapping is one to one
    series = series.dropna()
    if len(series) != len(series.drop_duplicates()):
        raise ReindexMapError("Reindex map is not one to one.")

    # Make the index the values, and the values the index
    df = series.reset_index()
    df = df.set_index(old_index_name)
    reindex_map = df[new_index_name]
    return reindex_map

def reindex_dataframe(df, reindex_map, new_index_name, keep_old):
    """Reindex a dataframe based on a mapping of the old index values to new ones.

    Parameters:
    df (pandas.DataFrame): The dataframe to reindex.
    reindex_map (dict or pandas Series): A dictionary or pandas Series with the old index values as the keys or index, and the new ones as the values.
    new_index_name (str): The desired name for the new index.
    keep_old (bool): Whether to retain the old index in the dataframe as a column.

    Returns:
    pandas.DataFrame: A copy of the given dataframe, with the new index.
    """
    if not df.index.isin(reindex_map.keys()).all(): # This works for a dict or a pandas Series, because Series have a .keys() attribute that's an alias for the index
        not_in = df.index[~(df.index.isin(reindex_map.keys()))]
        raise ReindexMapError(not_in)

    new_index = df.index.map(reindex_map.get)

    if keep_old:
        df = df.reset_index()
    df.index = new_index
    df.index.name = new_index_name
    df = df.sort_index()
    return df

def reindex_all_sample_id_to_patient_id(data_dict, reindex_map, additional_to_keep_col=[], skip=[]):
    """Reindex all the dataframes with Patient_IDs instead of Sample_IDs

    Parameters:
    data_dict (keys are str, values are pandas.DataFrame): The data dictionary to reindex
    reindex_map (dict or pandas Series): A dictionary or pandas Series with the old index values (Sample_IDs) as the keys or index, and the new ones (Patient_IDs) as the values. Must map for all existing index values in the entire dataset.
    additional_to_keep_col (list of str, optional): The function will already keep the old index as a column in the clinical dataframe. If you want it to do this for any other dataframes, put their names in this list.
    skip (str or list of str, optional): A list of dataframes to skip when reindexing.

    Returns:
    dict: The data dictionary, with all dataframes reindexed with Patient_IDs
    """
    if isinstance(skip, str): # If it's a single dataframe name, make it a list so we can treat everything the same
        skip = [skip]

    dfs_to_delete = []
    dfs_to_keep_col = ["clinical"] + additional_to_keep_col

    for name in data_dict.keys(): # Only loop over keys, to avoid changing the structure of the object we're looping over

        # Skip any specified to skip
        if name in skip:
            continue

        df = data_dict[name]
        df.index.name = "Sample_ID" # So that it's labeled properly when we keep it as a column in the clinical dataframe.
        keep_old = name in dfs_to_keep_col # Keep the old Patient_ID index as a column in the clinical dataframe (and any additionally specified dataframes), so we have a record of it.

        try:
            df = reindex_dataframe(df, reindex_map, "Patient_ID", keep_old)
        except ReindexMapError:
            warnings.warn(f"Error reindexing {name} dataframe. At least one Sample_ID did not have corresponding Patient_ID mapped in clinical dataframe. {name} dataframe not loaded.", FailedReindexWarning, stacklevel=3) # stacklevel=3 ensures that the warning is registered as originating from the file that called the __init__ function, instead of from here directly, because the former is more useful information.
            dfs_to_delete.append(name)
            continue

        data_dict[name] = df

    for name in dfs_to_delete: # Delete any dataframes that had issues reindexing
        del data_dict[name]

    return data_dict

def reformat_normal_patient_ids(data_dict, existing_identifier=None, existing_identifier_location=None):
    """Reformat the patient IDs for normal samples to be marked by an appended ".N"

    Parameters:
    data_dict (dict): The data dictionary for a dataset. All the tables must be indexed by Patient IDs.
    existing_identifier (str, optional): A normal sample identifier that already exists on the normal samples' patient IDs, which we will remove before adding the new identifier. Default of None will cause nothing to be removed.
    existing_identifier_location (str, optional): Either "start" or "end": Indicates whether the existing identifier is at the beginning or end of the normal samples' patient IDs, so we know which end to remove it from. Optional if nothing is passed to the existing_identifier parameter.

    Returns:
    dict: The data dictionary for the dataset, with normal samples' patient IDs reformatted in the specified dataframes.
    """

    # Check parameters
    if (existing_identifier is None and existing_identifier_location is not None) or (existing_identifier is not None and existing_identifier_location is None):
        raise CptacDevError("Parameters existing_identifier and existing_identifier_location must either both be None, or both not be None.")

    sample_statuses_old_index = data_dict["clinical"]["Sample_Tumor_Normal"] # We'll need this every time, and we need it with the un-reformatted index

    for name in data_dict.keys(): # Loop over the keys so we can edit the values without any issues

        df = data_dict[name]

        # Add in the tumor/normal statuses for these samples, if they aren't already in the table
        added_sample_statuses = False # So we can keep track of whether to drop the column when we're done
        if "Sample_Tumor_Normal" not in df.columns:
            df = join_col_to_dataframe(df, sample_statuses_old_index)
            added_sample_statuses = True

        df.index.name = "Patient_ID" # So we can easily access it when we've made it into a column
        df = df.reset_index() # This makes the Patient_ID index a column, so we can edit it.

        if (existing_identifier is not None) and (existing_identifier_location is not None): # There's an existing normal sample identifier to remove
            existing_length = len(existing_identifier)

            if existing_identifier_location == "start":
                df["Patient_ID"] = df["Patient_ID"].where(
                    cond=(~((df["Sample_Tumor_Normal"] == "Normal") & (df["Patient_ID"].str[0:existing_length] == existing_identifier))),
                    other=df["Patient_ID"].str[existing_length:]
                )

            elif existing_identifier_location == "end":
                df["Patient_ID"] = df["Patient_ID"].where(
                    cond=(~((df["Sample_Tumor_Normal"] == "Normal") & (df["Patient_ID"].str[-existing_length:] == existing_identifier))),
                    other=df["Patient_ID"].str[:-existing_length] # Note that we use the negative of the existing length, since we're working with the end of the string
                )

            else:
                raise CptacDevError("existing_identifier_location parameter must be either 'start' or 'end'")

        # Append ".N" to the patient IDs of normal samples
        df["Patient_ID"] = df["Patient_ID"].where(
            cond=(~(df["Sample_Tumor_Normal"] == "Normal")),
            other=df["Patient_ID"] + ".N"
        )

        # Set the index to the reformatted Patient IDs
        df = df.set_index("Patient_ID")

        # If we added the Sample_Tumor_Normal column, drop it
        if added_sample_statuses:
            if isinstance(df.columns, pd.MultiIndex):
                df = df.drop(columns="Sample_Tumor_Normal", level=0) # level=0 prevents a PerformanceWarning
            else:
                df = df.drop(columns="Sample_Tumor_Normal")

        # Put the dataframe with reformatted patient IDs back into the data dictionary
        data_dict[name] = df

    return data_dict

def join_col_to_dataframe(df, col):
    """Join a sample status column into a dataframe, automatically accounting for whether the dataframe has a column multiindex or not.

    Parameters:
    df (pandas.DataFrame): The dataframe to join the column into
    col (pandas Series): The column to join into the dataframe, with a matching index.

    Returns:
    pandas.DataFrame: The dataframe with the column joined in.
    """
    col_df = col.to_frame().copy(deep=True)

    # Make sure the columns axes all have the same name
    df.columns.name = "Name"
    col_df.columns.name = "Name"

    # If df has a column multiindex, edit the col_df column index to match, so we can join them
    if col_df.columns.names != df.columns.names:
        col_df.columns = add_index_levels(to=col_df.columns, source=df.columns)

    if col_df.columns.names != df.columns.names: # Just to make sure
        raise CptacDevError("col_df's column axes had levels not found in the dataframe's columns.")

    df = df.join(col_df, how="left") # We do a left join because we only want rows that exist in our dataframe

    return df

def standardize_axes_names(df):
    """For all dataframes in the given dictionary, sets the name of the index axes to "Patient_ID", because that's what they all are by that point, and sets the name of the column axes to "Name".

    Parameters:
    data_dict (dict): The dataframe dictionary of the dataset.

    Returns:
    dict: The dataframe dictionary, with the dataframe axes' names standardized. Keys are str of dataframe names, values are pandas.DataFrame
    """
    df.index.name = "Patient_ID"
    df.columns.name = "Name"

def add_index_levels(to, source, fill=""):
    """Add levels to the "to" index so it has all levels in the "source" index. The possible levels are, in this order: "Name", "Site", "Peptide", "Database_ID"

    Parameters:
    to (pandas.Index or pandas.MultiIndex): The index to add levels to.
    source (pandas.Index or pandas.MultiIndex): The index to match the levels of.
    fill (optional): Value to fill empty levels with. Default is an empty string, which allows us to select a column with just the first level. This is useful for boolean filters.

    Returns:
    pandas.MultiIndex: The levels of "to", with any levels from "source" that "to" didn't have originally.
    """
    to_set = set(to.names)
    source_set = set(source.names)
    if source_set <= to_set:
        return to # Because otherwise we'd just end up constructing a duplicate of "to", and who would want to do that?

    all_names = ["Name", "Site", "Peptide", "Database_ID"]
    levels = {}

    for name in all_names:
        if name in to.names:
            levels[name] = to.get_level_values(name)
        elif name in source.names:
            levels[name] = [fill for i in range(to.size)]

    new_columns = pd.MultiIndex.from_arrays(list(levels.values()), names=list(levels.keys()))
    return new_columns
