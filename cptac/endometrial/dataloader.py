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
import glob
import textwrap

def create_dataframe(path):
    """Create a dataframe for the data a file, parsed based on the data type.

    Parameters:
    path (str): The path from the current directory to the file to load the data from.

    Returns:
    pandas DataFrame: Dataframe of data in file, parsed depending on the data type.
    """
    path_elements = path.split(os.sep) # Get a list of the levels of the path
    file_name = path_elements[-1] # The last element will be the name of the file
    file_name_split = file_name.split(".")
    df_name = file_name_split[0] # Dataframe name will be the first section of file name; i.e. proteomics.txt.gz becomes proteomics
    file_extension = file_name_split[1]

    df = None
    print("Loading {} data...".format(df_name))
    if file_extension == "txt":
        #temp fix for reading error on clinical_v2:
        file = open(path, "r", errors="ignore")
        df = pd.read_csv(file, sep="\t", index_col=0)
        df = df.sort_index()

        df.name = df_name
    elif file_extension == "cct" or file_extension == "cbt":
        df = pd.read_csv(path, sep="\t", index_col=0)
        df = df.transpose()
        df = df.sort_index()

        df.name = df_name
    elif file_extension == "maf":
        df = pd.read_csv(path, sep = "\t")
        if "Tumor_Sample_Barcode" in df.columns:
            split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n = 1, expand = True)
            df["Tumor_Sample_Barcode"] = split_barcode[0]
        df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
        df = df.rename({"Tumor_Sample_Barcode":"Patient_Id","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')

        df.name = df_name + "_maf"
    else:
        error_message = "Error reading " + path
        raise IOError(error_message)
    return df

def warning():
    print("\n","******PLEASE READ******")
    warning = "WARNING: This data is under a publication embargo until July 1, 2019. CPTAC is a community resource project and data are made available rapidly after generation for community research use. The embargo allows exploring and utilizing the data, but the data may not be in a publication until July 1, 2019. Please see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter embargo() to open the webpage for more details."
    wrapped_list = textwrap.wrap(warning)
    for line in wrapped_list:
        print(line)

def get_dictionary():
    """Get a dictionary of terms for the dataset.

    Parameters: None

    Returns:
    dict: keys are terms, and values are definitions.
    """
    print("Loading Dictionary...")
    dictionary = {}
    dir_path = os.path.dirname(os.path.realpath(__file__))
    data_directory = dir_path + os.sep + "data" + os.sep
    with open(data_directory + "definitions.txt", "r") as dict_file:
        for line in dict_file.readlines():
            line = line.strip()
            line = line.split("\t")
            dictionary[line[0]] = line[1]
    return dictionary

def rename_df(data, old, new):
    """Rename a dataframe in the data dictionary, and set that name as the key, deleting the old key.

    Parameters:
    data (dict of str keys, pandas DataFrame values): The data dictionary containing the dataframe to be renamed.
    old (str): The old name.
    new (str): The new name.

    Returns: None
    """
    df = data[old]
    df.name = new
    data[new] = df
    del data[old]

def get_dataframes():
    """Load all of the endometrial dataframes, and format them properly.

    Parameters: None

    Returns:
    dict of pandas DataFrame: A dictionary containing all the endometrial dataframes as values, with their names as keys.
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    data_directory = dir_path + os.sep + "data" + os.sep
    path = data_directory + "*.*"
    files = glob.glob(path) # Puts all files into iterable variable
    files = [one_file for one_file in files if "definitions.txt" not in one_file] # Take out the definition.txt file, since we load it through get_dictionary()
    data = {}
    print("Loading cptac endometrial data:")
    for file in files: # Loops through files variable
        try:
            df = create_dataframe(file)
            data[df.name] = df # Maps dataframe name to dataframe
        except IOError:
            print("Error reading ", file)
            print("Check that all file names coincide with dataframe.py specs.")

    # Separate out clinical, derived_molecular, and experimental_setup dataframes
    all_clinical = data["clinical"]
    clinical = all_clinical[[
        'Proteomics_Participant_ID', 'Case_excluded',  'Proteomics_Tumor_Normal',  'Country',
        'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity',
        'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN', 'Clin_Stage_Dist_Mets-cM', 'Path_Stage_Dist_Mets-pM',
        'tumor_Stage-Pathological', 'FIGO_stage', 'LVSI', 'BMI', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site',
        'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm',   'Num_full_term_pregnancies']] 
    clinical = clinical.rename(columns={"Proteomics_Participant_ID":"Patient_ID"})
    clinical.name = "clinical"
    data["clinical"] = clinical

    derived_molecular = all_clinical.drop(['Proteomics_Participant_ID', 'Case_excluded',  'Proteomics_Tumor_Normal',  'Country',
        'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity',
        'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN', 'Clin_Stage_Dist_Mets-cM', 'Path_Stage_Dist_Mets-pM',
        'tumor_Stage-Pathological', 'FIGO_stage', 'LVSI', 'BMI', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site',
        'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm',   'Num_full_term_pregnancies',
        'Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs',
        'Proteomics_Aliquot_ID', 'Proteomics_OCT', 'WXS_normal_sample_type', 'WXS_normal_filename', 'WXS_normal_UUID', 'WXS_tumor_sample_type', 'WXS_tumor_filename',
        'WXS_tumor_UUID', 'WGS_normal_sample_type', 'WGS_normal_UUID', 'WGS_tumor_sample_type', 'WGS_tumor_UUID', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID',
        'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality'], axis=1)
    derived_molecular.name = "derived_molecular"
    data["derived_molecular"] = derived_molecular

    experimental_setup = all_clinical[['Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs',
        'Proteomics_Aliquot_ID', 'Proteomics_OCT', 'WXS_normal_sample_type', 'WXS_normal_filename', 'WXS_normal_UUID', 'WXS_tumor_sample_type', 'WXS_tumor_filename',
        'WXS_tumor_UUID', 'WGS_normal_sample_type', 'WGS_normal_UUID', 'WGS_tumor_sample_type', 'WGS_tumor_UUID', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID',
        'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality']]
    experimental_setup.name = "experimental_setup"
    data["experimental_setup"] = experimental_setup

    # Add Sample_ID column to somatic_mutations dataframe and make it the index
    clinical = data["clinical"] # We need the Patient_ID column from clinical, to map sample ids to patient ids. The sample ids are the clinical index, and the patient ids are in the Patient_ID column.
    patient_id_col = clinical.loc[clinical["Proteomics_Tumor_Normal"] == "Tumor", "Patient_ID"] # We only need to generate a sample id map for tumor samples, because normal samples are from the same person as tumor samples, and thus map to the same patient id.
    patient_id_col.index.name = "Sample_ID" # Label the index
    patient_id_df = patient_id_col.reset_index() # Make the index accessible as a column
    patient_id_df = patient_id_df.set_index("Patient_ID") # Set Patient_ID as the index, so we can look up a Sample_ID given a Patient_ID
    patient_id_map = patient_id_df["Sample_ID"] # Get the mapping as a series. Patient_ID will be the index.

    mutations = data["somatic_maf"] 
    mutations_patient_renamed = mutations.rename(columns={"Patient_Id":"Patient_ID"})
    mutations_patient_indexed = mutations_patient_renamed.set_index("Patient_ID") # Set the index as the Patient_ID column
    sample_id_col = [] # We're going to create a Sample_ID column for the mutations dataframe
    for patient_id in mutations_patient_indexed.index: 
        if patient_id in patient_id_map.index:
            sample_id_col.append(patient_id_map[patient_id]) # Get the sample id corresponding to the patient id
        else: # If there's not a corresponding sample ID for a patient ID, print an error message and return None
            print("Error mapping sample ids in somatic mutations dataframe. Patient_ID {} did not have corresponding Sample_ID mapped in clinical dataframe. Data loading aborted.".format(patient_id))
            return 
    mutations_with_sample = mutations_patient_indexed.assign(Sample_ID=sample_id_col) # Add in the Sample_ID column
    mutations_sample_indexed = mutations_with_sample.set_index("Sample_ID") # Make the Sample_ID column the index
    mutations_sample_indexed.name = mutations.name
    data["somatic_maf"] = mutations_sample_indexed

    # Drop unfiltered samples
    clinical = data["clinical"]
    cases_to_drop = clinical[clinical["Case_excluded"] == "Yes"].index

    for name in data.keys(): # Loop over keys instead of dictionary directly, so we're not altering the structure we're looping over
        df = data[name]
        df_filtered = df.drop(index=cases_to_drop, errors="ignore")
        df_filtered.name = df.name
        data[name] = df_filtered

    # Drop Case_excluded column from clinical
    clinical = data["clinical"]
    clinical_no_case_excluded = clinical.drop(columns=["Case_excluded"])
    clinical_no_case_excluded.name = clinical.name
    data["clinical"] = clinical_no_case_excluded

    # Fix names
    rename_df(data, old="transcriptomics_linear", new="transcriptomics")
    rename_df(data, old="phosphoproteomics_site", new="phosphoproteomics")
    rename_df(data, old="transcriptomics_circular", new="circular_RNA")
    rename_df(data, old="somatic_maf", new="somatic_mutation")
    rename_df(data, old="somatic_binary", new="somatic_mutation_binary")

    # Rename indicies to "Sample_ID", since that's what they all are.
    for name in data.keys():
        df = data[name]
        df.index.name = "Sample_ID"
        data[name] = df

    warning()
    return data

