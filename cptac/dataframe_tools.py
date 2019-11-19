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
from .exceptions import CptacDevError, ReindexMapError

def unionize_indices(dataset):
    """Return a union of all indices in a dataset, without duplicates.

    Parameters:
    dataset (dict of str: pandas DataFrame): The data dictionary containing the dataset.

    Returns:
    pandas Index: Union of all indices in the dataset, without duplicates.
    """
    indices = [df.index for df in dataset.values()]
    master_index = pd.Index([])
    for index in indices:
        master_index = master_index.union(index)
        master_index = master_index.drop_duplicates()
    return master_index

def generate_sample_status_col(df, normal_test):
    """Create a sample status column, called Sample_Tumor_Normal, for a dataframe.

    Parameters:
    df (pandas DataFrame): The dataframe to create a Sample_Status column for, indexed with Patient_IDs.
    normal_test (function): A function that takes a given Patient_ID and returns a bool indicating whether it corresponds to a normal sample.

    Returns:
    pandas Series: A sample status column for the dataframe.
    """
    sample_status_array = np.where(df.index.map(normal_test), "Normal", "Tumor")
    sample_status_col = pd.Series(data=sample_status_array, index=df.index.copy())
    return sample_status_col

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

def generate_sample_id_map(master_index):
    """Generate sample ids for all samples in a dataset, and map them to the existing index.

    Parameters:
    master_index (pandas Index): Union of all the indices in the dataset, without duplicates.

    Returns:
    dict: A dictionary with the old index values as the keys, and their corresponding sample ids as the values.
    """
    sample_id_dict = {}
    for i in range(len(master_index)):
        patient_id = master_index[i]
        sample_id_dict[patient_id] = "S{:0>3}".format(i + 1) # Use string formatter to give each sample id the format S*** filled with zeroes, e.g. S001, S023, or S112
    return sample_id_dict

def reindex_dataframe(df, reindex_map, new_index_name, keep_old):
    """Reindex a dataframe based on a mapping of the old index values to new ones.

    Parameters:
    df (pandas DataFrame): The dataframe to reindex.
    reindex_map (dict or pandas Series): A dictionary or pandas Series with the old index values as the keys or index, and the new ones as the values.
    new_index_name (str): The desired name for the new index.
    keep_old (bool): Whether to retain the old index in the dataframe as a column.

    Returns:
    pandas DataFrame: A copy of the given dataframe, with the new index.
    """
    if not df.index.isin(reindex_map.keys()).all(): # This works for a dict or a pandas Series, because Series have a .keys() attribute that's an alias for the index
        not_in = df.index[~(df.index.isin(reindex_map.keys()))]
        raise ReindexMapError(not_in)

    new_index = df.index.map(reindex_map.get)

    if keep_old:
        df = df.reset_index() # This gives the dataframe a numerical index and makes the old index a column, so it's not dropped when we set the new index.
    df.index = new_index
    df.index.name = new_index_name
    df = df.sort_index()
    return df

def reindex_all(data_dict, master_index, additional_to_keep_col=[]):
    """Reindex all the dataframes with Sample_IDs instead of Patient_IDs

    Parameters:
    data_dict (keys are str, values are pandas DataFrames): The data dictionary to reindex
    master_index (pandas Index): An index of all patient IDs that exist in the dataset

    Returns:
    dict: The data dictionary, with all dataframes reindexed with Sample_IDs
    """
    dfs_to_delete = []
    dfs_to_keep_col = ["clinical"] + additional_to_keep_col
    sample_id_dict = generate_sample_id_map(master_index) # Generate a sample ID for each patient ID

    for name in data_dict.keys(): # Only loop over keys, to avoid changing the structure of the object we're looping over
        df = data_dict[name]
        df.index.name = "Patient_ID"
        keep_old = name in dfs_to_keep_col # Keep the old Patient_ID index as a column in the clinical dataframe (and any additionally specified dataframes), so we have a record of it.
        try:
            df = reindex_dataframe(df, sample_id_dict, "Sample_ID", keep_old)
        except ReindexMapError:
            warnings.warn(f"Error mapping sample ids in {name} dataframe. At least one Patient_ID did not have corresponding Sample_ID mapped in clinical dataframe. {name} dataframe not loaded.", FailedReindexWarning, stacklevel=3) # stacklevel=3 ensures that the warning is registered as originating from the file that called the __init__ function, instead of from here directly, because the former is more useful information.
            dfs_to_delete.append(name)
            continue

        data_dict[name] = df

    for name in dfs_to_delete: # Delete any dataframes that had issues reindexing
        del data_dict[name]

    return data_dict

def reformat_normal_patient_ids(data_dict, existing_identifier=None, existing_identifier_location=None, additional_dfs_to_reformat=None):
    """WARNING: USE THIS FUNCTION ONLY AFTER ALL DATAFRAMES IN THE DATASET HAVE BEEN REINDEXED WITH SAMPLE IDS.
    Reformat the patient IDs for normal samples to be marked by a prepended "N."

    Parameters:
    data_dict (dict): The data dictionary for a dataset
    existing_identifier (str, optional): A normal sample identifier that already exists on the normal samples' patient IDs, which we will remove before adding the new identifier. Default of None will cause nothing to be removed.
    existing_identifier_location (str, optional): Either "start" or "end": Indicates whether the existing identifier is at the beginning or end of the normal samples' patient IDs, so we know which end to remove it from. Optional if nothing is passed to the existing_identifier parameter.
    additional_dfs_to_reformat (str or list of str, optional): A single name or list of the names of additional dataframes within the data dictionary in which to reformat the patient ID column. The column is always reformatted in the clinical dataframe, so that's what will be done if this is left at the default of None.

    Returns:
    dict: The data dictionary for the dataset, with normal samples' patient IDs reformatted in the specified dataframes.
    """
    # Check parameters
    if (existing_identifier is None and existing_identifier_location is not None) or (existing_identifier is not None and existing_identifier_location is None):
        raise CptacDevError("Parameters existing_identifier and existing_identifier_location must either both be None, or both be not None.")

    dfs_to_reformat = ["clinical"] # We will always at least reformat the patient IDs in the clinical dataframe

    if additional_dfs_to_reformat is not None:
        if isinstance(additional_dfs_to_reformat, str):
            additional_dfs_to_reformat = [additional_dfs_to_reformat]
            
        dfs_to_reformat += additional_dfs_to_reformat
        dfs_to_reformat = set(dfs_to_reformat) # In case they duplicated the clinical dataframe in their input list

    clinical = data_dict["clinical"] # We'll need this every time for the Sample_Tumor_Normal column

    for name in dfs_to_reformat:
        if name not in data_dict.keys():
            raise CptacDevError("Invalid dataframe name passed to additional_dfs_to_reformat parameter.")

        df = data_dict[name]

        if (existing_identifier is not None) and (existing_identifier_location is not None): # Remove the existing normal sample identifier
        # We do this in a separate step because not all datasets have an existing identifier, and some datasets could have some patient IDs with and existing identifier but others without one. This would only remove it from the ones that have it.
            existing_length = len(existing_identifier)

            if existing_identifier_location == "start":
                df.loc[(clinical.loc[df.index, "Sample_Tumor_Normal"] == "Normal") & (df["Patient_ID"].str[0:existing_length] == existing_identifier), "Patient_ID"] = df["Patient_ID"].str[existing_length:]

            elif existing_identifier_location == "end":
                df.loc[(clinical.loc[df.index, "Sample_Tumor_Normal"] == "Normal") & (df["Patient_ID"].str[-existing_length:] == existing_identifier), "Patient_ID"] = df["Patient_ID"].str[:-existing_length] # Note that we use the negative of the existing length, since we're working with the end of the string

            else:
                raise CptacDevError("existing_identifier_location parameter must be either 'start' or 'end'")

        # Prepend "N." to the patient IDs of normal samples
        df.loc[clinical.loc[df.index, "Sample_Tumor_Normal"] == "Normal", "Patient_ID"] = "N." + df["Patient_ID"]

        # Put the dataframe with reformatted patient IDs back into the data dictionary
        data_dict[name] = df

    print("Yayyy!")
    return data_dict


def standardize_axes_names(data_dict):
    """For all dataframes in the given dictionary, sets the name of the index axes to "Sample_ID", because that's what they all are by that point, and sets the name of the column axes to "Name".

    Parameters:
    data_dict (dict): The dataframe dictionary of the dataset.

    Returns:
    dict: The dataframe dictionary, with the dataframe axes' names standardized. Keys are str of dataframe names, values are pandas DataFrames
    """
    # Rename indices to "Sample_ID", since that's what they all are.
    for name in data_dict.keys(): # Loop over the keys so we can alter the values without any issues
        df_rename_index = data_dict[name]
        df_rename_index.index.name = "Sample_ID"
        data_dict[name] = df_rename_index

    # Set name of column axis to "Name" for all dataframes
    for name in data_dict.keys(): # Loop over the keys so we can alter the values without any issues
        df = data_dict[name]
        df.columns.name = "Name"
        data_dict[name] = df

    return data_dict
