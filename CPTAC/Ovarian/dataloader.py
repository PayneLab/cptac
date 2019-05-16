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
import os
import glob
import textwrap

def create_dataframe(path):
    """Create a dataframe for the data a file, parsed based on the data type.

    Parameters:
    path (str): The path from the current directory to the file to load the data from.

    Returns:
    pandas DataFrame: Dataframe of data in file, parsed depending on the data type.
    """
    path_elements = path.split(os.sep) # Get a list of all the levels of the path
    file_name = path_elements[-1] # The last element will be the name of the file
    df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

    print("Loading {} data...".format(df_name))
    if df_name == "proteomics":
        df = pd.read_csv(path,sep="\t", index_col = 0)
        df = df[df["hgnc_symbol"].notnull()] # Drops all nan values in hgnc_symbol column
        df = df.set_index("hgnc_symbol")
        df = df.sort_index()
        df = df.transpose()
        c_index = df.index[0:83].str[1:] # Drops letter off all indices with "C"
        index = c_index.append(df.index[83:])
        df = df.set_index(index)

        # Drop all OV_QC* samples--they're quality control samples not relevant for data analysis
        idx = df.index.values.tolist()
        idx_to_drop = [id for id in idx if id.startswith('OV_QC')]
        df = df.drop(idx_to_drop) 

        # Name and return the dataframe
        df.name = df_name
        return df
    elif df_name == "clinical":
        df = pd.read_csv(path, sep="\t")
        df = df.set_index("PPID")
        df = df[~df.index.duplicated(keep="first")]
        df.name = df_name
        return df
    elif df_name == "phosphoproteomics":
        df = pd.read_csv(path, sep = "\t",index_col = 0)
        df = df[df["site"].notnull()] # Drops all nan values in site column
        df = df.drop(["refseq_peptide","Peptide"],axis=1)
        df = df.set_index("site")
        df = df.sort_index()
        df = df.transpose()
        c_index = df.index[0:83].str[1:] # Drops letter off all indices with "C"
        index = c_index.append(df.index[83:])
        df = df.set_index(index)
        
        # Drop all OV_QC* samples--they're quality control samples not relevant for data analysis
        idx = df.index.values.tolist()
        idx_to_drop = [id for id in idx if id.startswith('OV_QC')]
        df = df.drop(idx_to_drop) 

        # Name and return the dataframe
        df.name = df_name
        return df
    elif df_name == "transcriptomics":
        df = pd.read_csv(path, sep="\t", index_col=0)
        df = df.sort_index()
        df = df.transpose()
        df = df.sort_index()
        df = df.drop(columns = df.columns[0:23]) # Drop all date values until new data is uploaded
        df.name = df_name
        return df
    elif df_name == "cnv":
        df = pd.read_csv(path, sep="\t", index_col=0)
        df = df.sort_index()
        df = df.transpose()
        df = df.sort_index()
        df.name = "CNV"
        return df
    elif df_name == "somatic_38":
        df = pd.read_csv(path, sep = "\t")
        if "Tumor_Sample_Barcode" in df.columns:
            split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n = 1, expand = True)
            df["Tumor_Sample_Barcode"] = split_barcode[0]
        parsedDf = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
        parsedDf = parsedDf.rename({"Tumor_Sample_Barcode":"Patient_Id","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
        parsedDf = parsedDf.set_index("Patient_Id")
        parsedDf.name = 'somatic_mutation'
        return parsedDf
    else:
        print("Error reading", path)

def warning():
    print("\n","******PLEASE READ******")
    warning = "WARNING: This data is under a publication embargo until June 1, 2019. CPTAC is a community resource project and data are made available rapidly after generation for community research use. The embargo allows exploring and utilizing the data, but the data may not be in a publication until June 1, 2019. Please see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter embargo() to open the webpage for more details."
    wrapped_list = textwrap.wrap(warning)
    for line in wrapped_list:
        print(line)

def set_sample_id_index(df, sample_id_dict, drop_patient_ids): # private
    """Replaces a patient ID index with a sample ID index for a dataframe.

    Parameters:
    df (pandas DataFrame): The dataframe to change the index of.
    sample_id_dict (dict): A dictionary where the keys are the patient ids, and the values are the corresponding sample ids. Every value in the original dataframe's index must match a key in this dictionary.
    drop_patient_ids (bool): Whether to drop the patient id column after it is no longer the index.

    Returns:
    pandas DataFrame: The original dataframe, with the patient id index replaced with a sample id index.
    """
    sample_id_column = []
    for row in df.index:
        if row in sample_id_dict.keys():
            sample_id_column.append(sample_id_dict[row]) 
        else:
            print("Error mapping sample ids in {0} dataframe. Patient_ID {1} did not have corresponding Sample_ID mapped in clinical dataframe. {0} dataframe not loaded.".format(df.name, row))
            return
    return_df = df.assign(Sample_ID=sample_id_column)
    if not drop_patient_ids:
        old_index_name = df.index.name
        if old_index_name is None:
            old_index_name = 'index'
        return_df = return_df.reset_index() # This gives the dataframe a default numerical index and makes the old index a column, which prevents it from being dropped when we set Sample_ID as the index.
        return_df = return_df.rename(columns={old_index_name:'Patient_ID'}) # Rename the old index as Patient_ID
    return_df = return_df.set_index('Sample_ID') # Make the Sample_ID column the index
    return_df.name = df.name
    return return_df

def get_dataframes():
    """Load all of the Ovarian dataframes, and format them properly.

    Parameters: None

    Returns:
    dict of pandas DataFrame: A dictionary containing all of the endometrial dataframes as values, with their names as keys.
    """
    dir_path = os.path.dirname(os.path.realpath(__file__)) #gets path to CPTAC package
    data_directory = dir_path + os.sep + "Data" + os.sep #appends Data to path
    path = data_directory + os.sep + "*.*" #appends "*.*" to path, which looks for all files
    files = glob.glob(path) #puts all files into iterable variable

    # Load the dataframes
    print("Loading CPTAC Ovarian data:")
    data = {}
    for file in files: #loops through files variable
        df = create_dataframe(file)
        data[df.name] = df #maps dataframe name to dataframe

    # Get a union of all dataframes' indicies, with duplicates removed
    indicies = [df.index for df in data.values()]
    master_index = pd.Index([])
    for index in indicies:
        master_index = master_index.union(index)
        master_index = master_index.drop_duplicates()

    # Generate a sample ID for each patient ID
    sample_id_dict = {}
    for i in range(len(master_index)):
        patient_id = master_index[i]
        sample_id_dict[patient_id] = "S{:0>3}".format(i + 1) # Use string formatter to give each sample id the format S*** filled with zeroes, e.g. S001 or S023

    # Put a mapping in the clinical dataframe of all patient ids to their sample ids, including patient ids for samples not originally in the clinical dataframe. Then, give the clinical dataframe a sample id index.
    master_df = pd.DataFrame(index=master_index)
    master_clinical = data['clinical'].join(master_df, how='outer') # Do an outer join with the clinical dataframe, so that clinical has a row for every sample in the dataset
    master_clinical.name = data["clinical"].name
    master_clinical = set_sample_id_index(master_clinical, sample_id_dict, drop_patient_ids=False) # Replace the patient id index with a sample id index in the clinical dataframe. Keep the patient ids so we can maps sample ids to their patient ids.

    # Add a column, Sample_Tumor_Normal, indicating whether each sample was a tumor or normal sample. Normal samples have a Patient_ID that begins with 'N'.
    clinical_status_col = []
    for sample in master_clinical["Patient_ID"]:
        if sample[0] == 'N':
            clinical_status_col.append("Normal")
        else:
            clinical_status_col.append("Tumor")
    master_clinical.insert(1, "Sample_Tumor_Normal", clinical_status_col)
    data['clinical'] = master_clinical # Replace the clinical dataframe in the data dictionary with our new and improved version!

    # Give the other dataframes Sample_ID indicies, but don't keep the old index, since we have a mapping in the clinical dataframe of all sample ids to their patient ids.
    for name in data.keys(): # Only loop over keys, to avoid changing the structure of the object we're looping over
        if name != 'clinical':
            data[name] = set_sample_id_index(data[name], sample_id_dict, drop_patient_ids=True)

    warning() #displays warning

    return data
