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

def validate_version(version, dataset_path, index):
    """Check that a given version number is valid. If version is "latest", check that index and installed latest match.

    Parameters:
    version (str): The version number to validate.
    dataset_path (str): The path to the dataset the version is for.
    index (dict): The parsed index for the dataset.

    Returns:
    str: The version number, if valid input. Else None.
    """
    if version in index.keys():
        return version
    elif version.lower() == "latest":
        index_latest = max(index.keys(), key=float) # See what the highest version in the index is
        latest_installed = get_latest_installed(dataset_path)
        if (index_latest == latest_installed) or (latest_installed is None):
            return index_latest
        else:
            print(f"Ambiguous request for latest version. Latest version in index is {index_latest}, but latest version installed locally is {latest_installed}. To download the latest version in the index, run cptac.sync with '{index_latest}' as the version parameter. To run your function with the latest version that is installed locally, run the function with '{latest_installed}' as the version parameter.")
            return None
    else:
        print(f"{version} is an invalid version for this dataset. Valid versions: {', '.join(index.keys())}")
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
