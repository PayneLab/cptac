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

import numpy as np
import pandas as pd
import os
import re
import glob

def create_dataframe(path):
    """Create a dataframe for the data a file, parsed based on the data type.

    Parameters:
    path (str): The path from the current directory to the file to load the data from.

    Returns:
    pandas DataFrame: Dataframe of data in file, parsed depending on the data type.
    """
    path_elements = path.split(os.sep) # Get a list of the levels of the path
    file_name = path_elements[-1] # The last element will be the name of the file
    df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

    df = None
    print("Loading {} data...".format(df_name))
    if df_name == "clinical" or df_name.split("_")[0] == "proteomics" or df_name.split("_")[0] == "transcriptomics":
        df = pd.read_csv(path, sep="\t",index_col=0)
        df = df.transpose()
        df.name = df_name
    elif df_name.split("_")[0] == "phosphoproteomics": #column names are SLC12A8_T485_A0AV02:T485, should this be changed?
        df = pd.read_csv(path, sep="\t",index_col=0)
        df = df.transpose()
        df.name = df_name
    elif df_name == "miRNA": #column names are hsa-let-7a-2-3p, should this be changed?
        df = pd.read_csv(path, sep="\t",index_col=0)
        df = df.transpose()
        df.name = df_name
    elif df_name.split("_")[0] == "mutation":
        if file_name.split(".")[1] == "cbt":
            df = pd.read_csv(path, sep="\t",index_col=0)
            df = df.transpose()
            df.name = "somatic_" + df_name
        elif file_name.split(".")[1] == "txt":
            df = pd.read_csv(path, sep="\t")
            df = df.sort_values(by="SampleID")
            df = df[["SampleID","Gene","Variant_Type","Protein_Change"]]
            df = df.rename({"Variant_Type":"Mutation","Protein_Change":"Location"},axis="columns")
            df.name = "somatic_" + df_name
    else:
        error_message = "Error reading " + path
        raise IOError(error_message)
    return df

def set_sample_id_index(df, sample_id_dict, drop_patient_ids):
    """Replaces a patient ID index with a sample ID index for a dataframe. Preserves dataframe name.

    Parameters:
    df (pandas DataFrame): The dataframe to change the index of.
    sample_id_dict (dict): A dictionary where the keys are the patient ids, and the values are the corresponding sample ids. Every value in the original dataframe's index must match a key in this dictionary.
    drop_patient_ids (bool): Whether to drop the patient id column after it is no longer the index.

    Returns:
    pandas DataFrame: The original dataframe, with the patient id index replaced with a sample id index, and the dataframe name preserved.
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
    return_df = return_df.sort_index() # Get everything in order
    return_df.name = df.name
    return return_df

def get_dataframes():
    """Load all of the colon dataframes, and format them properly.

    Parameters: None

    Returns:
    dict of pandas DataFrame: A dictionary containing all the colon dataframes as values, with their names as keys.
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    data_directory = dir_path + os.sep + "data" + os.sep
    path = data_directory + "*.*"
    files = glob.glob(path) # Puts all files into iterable variable
    data = {}
    print("Loading cptac colon data:")
    for file in files: # Loops through files variable
        try:
            print(file)
            df = create_dataframe(file)
            data[df.name] = df # Maps dataframe name to dataframe
        except IOError:
            print("Error reading ", file)
            print("Check that all file names coincide with dataframe.py specs.")

    # Separate clinical and derived molecular dataframes
    all_clinical_data = data.get("clinical")
    clinical_df = all_clinical_data.drop(columns=['StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'TumorPurity','immuneSubtype', 'CIN', 'Integrated.Phenotype'])
    clinical_df.name = "clinical"
    derived_molecular_df = all_clinical_data[['StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'TumorPurity', 'immuneSubtype', 'CIN', 'Integrated.Phenotype']]
    derived_molecular_df.name = "derived_molecular"

    # Put them in our data dictionary
    data["clinical"] = clinical_df # Replaces original clinical dataframe
    data["derived_molecular"] = derived_molecular_df

    # Combine the two proteomics dataframes
    prot_tumor = data.get("proteomics_tumor")
    prot_normal = data.get("proteomics_normal") # Normal entries are marked with 'N' on the end of the ID
    prot_combined = prot_tumor.append(prot_normal)
    prot_combined.name = "proteomics"
    data[prot_combined.name] = prot_combined
    del data["proteomics_tumor"]
    del data["proteomics_normal"]

    # Get phosphoproteomics dataframes
    phos_tumor = data.get("phosphoproteomics_tumor")
    phos_normal = data.get("phosphoproteomics_normal") # Normal entries are not marked

    # Mark entries in phosphoproteomics_normal dataframe with an N at the end of the ID, to match proteomics_normal
    phos_normal_indicies = phos_normal.index.values.tolist()
    for i in range(len(phos_normal_indicies)):
        index = phos_normal_indicies[i]
        index_marked = index + 'N'
        phos_normal_indicies[i] = index_marked
    new_phos_index = pd.Index(phos_normal_indicies)
    phos_normal = phos_normal.set_index(new_phos_index)

    # Combine the two phosphoproteomics dataframes into one dataframe. 
    phos_combined = phos_tumor.append(phos_normal)
    phos_combined.name = 'phosphoproteomics'

    # Add the new phosphoproteomics dataframe to the dict, and delete the old ones.
    data[phos_combined.name] = phos_combined
    del data["phosphoproteomics_tumor"]
    del data["phosphoproteomics_normal"]

    # Rename the somamtic_mutation dataframe's "SampleID" column to "PatientID", then set that as the index, to match the other dataframes
    new_somatic = data["somatic_mutation"]
    new_somatic = new_somatic.rename(columns={"SampleID":"Patient_ID"})
    new_somatic = new_somatic.set_index("Patient_ID")
    new_somatic.name = "somatic_mutation"
    data["somatic_mutation"] = new_somatic

    # Get a union of all dataframes' indicies, with duplicates removed
    indicies = [df.index for df in data.values()]
    master_index = pd.Index([])
    for index in indicies:
        master_index = master_index.union(index)
        master_index = master_index.drop_duplicates()

    # Sort this master_index so all the samples with an N suffix are last. Because the N is a suffix, not a prefix, this is kind of messy.
    status_df = pd.DataFrame(master_index, columns=['Patient_ID']) # Create a new dataframe with the master_index as a column called "Patient_ID"
    status_col = []
    for index in master_index:
        if index[-1] == 'N':
            status_col.append("Normal")
        else:
            status_col.append("Tumor")
    status_df = status_df.assign(Status=status_col)
    status_df = status_df.sort_values(by=['Status', 'Patient_ID'], ascending=[False, True])
    master_index = status_df["Patient_ID"].tolist()

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

    # Add a column, Sample_Tumor_Normal, indicating whether each sample is a tumor or normal sample. Samples with a Patient_ID ending in N are normal.
    clinical_status_col = []
    for sample in master_clinical["Patient_ID"]:
        if sample[-1] == 'N':
            clinical_status_col.append("Normal")
        else:
            clinical_status_col.append("Tumor")
    master_clinical.insert(1, "Sample_Tumor_Normal", clinical_status_col)

    data['clinical'] = master_clinical # Replace the clinical dataframe in the data dictionary with our new and improved version!

    # Give the other dataframes Sample_ID indicies, but don't keep the old index, since we have a mapping in the clinical dataframe of all sample ids to their patient ids.
    for name in data.keys(): # Only loop over keys, to avoid changing the structure of the object we're looping over
        if name != 'clinical':
            data[name] = set_sample_id_index(data[name], sample_id_dict, drop_patient_ids=True)

    return data
