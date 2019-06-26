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

import hashlib
import os

def get_dataset_path(dataset):
    """Get the path to the main directory for a dataset.

    Parameters:
    dataset (str): The path to get the directory for. Must be all lowercase.

    Returns:
    str: The path to the main directory of the specified dataset.
    """
    path_here = os.path.abspath(os.path.dirname(__file__))
    if dataset == "gbm":
        dataset_dir = "mockup_data_gbm"
    else:
        dataset_dir = f"data_{dataset}"
    dataset_path = os.path.join(path_here, dataset_dir)
    if os.path.isdir(dataset_path):
        return dataset_path
    else:
        print(f"{dataset} is not a valid dataset.")
        return None

def validate_version(version, dataset, dataset_path, index, use_context):
    """Check that a given version number is valid. If version is "latest", check that index and installed latest match.

    Parameters:
    version (str): The version number to validate.
    dataset (str): The name of the dataset we're validating the version for.
    dataset_path (str): The path to the dataset the version is for.
    index (dict): The parsed index for the dataset.
    use_context (str): Either "sync" or "load", depending on whether the function is being called as part of a data sync or a dataset loading. Allows for more detailed error messages. Pass None if you don't want more detailed error messages.

    Returns:
    str: The version number, if valid input. If they passed the "latest" keyword but their latest installed didn't match the latest in the index, then return "ambiguous_latest". Else return None.
    """
    # See what the highest version in the index is
    index_latest = max(index.keys(), key=float)

    # Parse and validate the version they passed
    if version in index.keys():
        if float(version) < float(index_latest): # Print a warning if they're using an old version
            print(f"WARNING: You are using an old data version. Latest is {index_latest}. You are using {version}.")
        return version

    elif version.lower() == "latest":
        latest_installed = get_latest_installed(dataset_path)
        if (index_latest == latest_installed) or (latest_installed is None):
            return index_latest
        else: # If their latest installed version is different from the latest version recorded in the index, then we don't know which one they meant when they passed "latest".
            if use_context == "sync":
                message = f"You requested to sync the {dataset} dataset. Latest version is {index_latest}, which is not installed locally. To download it, run \"cptac.sync(dataset='{dataset}', version='{index_latest}')\". To instead sync the older version that is already installed, run \"cptac.sync(dataset='{dataset}', version='{latest_installed}')\"."
            elif use_context == "load":
                message = f"You requested to load the {dataset} dataset. Latest version is {index_latest}, which is not installed locally. To download it, run \"cptac.sync(dataset='{dataset}', version='{index_latest}')\". You will then be able to load the latest version by calling \"cptac.{dataset.title()}()\". Or, to instead load the older version that is already installed, call \"cptac.{dataset.title()}(version='{latest_installed}')\"."
            else:
                message = f"You requested the latest version. Latest version is {index_latest}, which is not installed locally. To download it, run \"cptac.sync(dataset='{dataset}', version='{index_latest}')\"."

            print(message)
            return None
    else:
        print(f"{version} is an invalid version for the {dataset} dataset. Valid versions: {', '.join(index.keys())}")
        return None

def get_latest_installed(dataset_path):
    """Return the latest version number installed in a dataset directory.

    Parameters:
    dataset_path (str): The path to the dataset of interest.

    Returns:
    str: The latest version installed locally.
    """
    dirs = [dir.strip() for dir in os.listdir(dataset_path)
                if os.path.isdir(os.path.join(dataset_path, dir))]

    dataset_dir = dataset_path.split(os.sep)[-1]
    dataset_name = dataset_dir.split('_')[1] # For example, we get 'endometrial' from 'data_endometrial'
    version_dir_prefix = dataset_name + '_v'
    versions = [dir.replace(version_dir_prefix, '') for dir in dirs
                    if dir.startswith(version_dir_prefix)]
    if len(versions) == 0:
        return None
    latest_installed = max(versions, key=float)
    return latest_installed

def get_index(dataset_path):
    """Get the index for a dataset, as a nested dictionary

    Parameters:
    dataset_path (str): The path to the dataset you want the index of.

    Returns:
    dict: The index, as a nested dictionary.
    """
    index_file = "index.txt"
    index_path = os.path.join(dataset_path, index_file)
    with open(index_path, 'r') as index_file:
        index_lines = index_file.readlines()

    index = {}
    version = None
    for line in index_lines:
        line = line.strip()
        if line.startswith('#'):
            version = line[1:]
            index[version] = {}
        else:
            line_list = line.split('\t')
            file_name = line_list[0]
            file_hash = line_list[1]
            file_url = line_list[2]
            index[version][file_name] = {}
            index[version][file_name]["hash"] = file_hash
            index[version][file_name]["url"] = file_url
    return index

def parse_tsv_dict(path):
    """Read in a dictionary from the given two column tsv file.

    Parameters:
    path (str): The path to the two column tsv file.

    Returns:
    dict: The tsv file read into a dictionary.
    """
    with open(path, 'r') as data_file:
        lines = data_file.readlines()

    data_dict = {}
    for line in lines:
        line_list = line.strip().split("\t")
        key = line_list[0]
        value = line_list[1]
        data_dict[key] = value

    return data_dict

def hash_file(path):
    """Return the md5 hash for the file at the given path. Return None if file doesn't exist.

    Parameters:
    path (str): The absolute path to the file to hash.

    Returns:
    str: The hash for the file.
    """
    if os.path.isfile(path):
        with open(path, 'rb') as file_obj:
            hash = hash_bytes(file_obj.read())
        return hash
    else:
        return None

def hash_bytes(bytes):
    """Hash the given bytes.

    Parameters:
    bytes (bytes): The bytes to has.

    Returns:
    str: The hash for the bytes.
    """
    hasher = hashlib.md5()
    hasher.update(bytes)
    hash = hasher.hexdigest()
    return hash
