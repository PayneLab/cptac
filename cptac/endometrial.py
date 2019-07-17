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
import os
from .dataset import DataSet
from .file_download import update_index
from .file_tools import validate_version, get_version_files_paths
from .dataframe_tools import *

class Endometrial(DataSet):

    def __init__(self, version="latest"):
        """Load all of the endometrial dataframes as values in the self._data dict variable, with names as keys, and format them properly."""

        # Call the parent DataSet __init__ function, which initializes self._data and other variables we need
        super().__init__("endometrial")

        # Update the index, if possible. If there's no internet, update_index will return False, but we don't care in this context.
        update_index(self._cancer_type)

        # Validate the index
        self._version = validate_version(version, self._cancer_type, use_context="init")
        if self._version is None: # Validation error. validate_version already printed an error message.
            return

        # Get the paths to all the data files
        data_files = [
            "acetylproteomics.cct.gz", 
            "clinical.txt", 
            "CNA.cct.gz", 
            "definitions.txt",
            "miRNA.cct.gz", 
            "phosphoproteomics_gene.cct.gz", 
            "phosphoproteomics_site.cct.gz", 
            "proteomics.cct.gz", 
            "somatic_binary.cbt.gz", 
            "somatic.maf.gz", 
            "transcriptomics_circular.cct.gz", 
            "transcriptomics_linear.cct.gz"]

        data_files_paths = get_version_files_paths(self._cancer_type, self._version, data_files)
        if data_files_paths is None: # Data error. get_version_files_paths already printed an error message.
            return 

        # Load the data files into dataframes in the self._data dict
        loading_msg = "Loading dataframes"
        for file_path in data_files_paths: 

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            df_name = file_name.split(".")[0] # Dataframe name will be the first section of file name; i.e. proteomics.txt.gz becomes proteomics

            # Load the file, based on what it is
            if file_name == "clinical.txt":
                # Fix for reading error on clinical.txt:
                with open(file_path, "r", errors="ignore") as clinical_file:
                    df = pd.read_csv(clinical_file, sep="\t", index_col=0)
                df = df.sort_index()
                self._data[df_name] = df # Maps dataframe name to dataframe

            elif file_name == "definitions.txt":
                with open(file_path, "r") as definitions_file:
                    for line in definitions_file.readlines():
                        line = line.strip()
                        line = line.split("\t")
                        term = line[0]
                        definition = line[1]
                        self._definitions[term] = definition

            elif file_name == "somatic.maf.gz":
                df = pd.read_csv(file_path, sep = "\t")
                split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n = 1, expand = True) # The first part of the barcode is the patient id, which we need want to make the index
                df["Tumor_Sample_Barcode"] = split_barcode[0]
                df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
                df = df.rename({"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
                df = df.sort_values(by=["Patient_ID", "Gene"])
                df = df.set_index("Patient_ID")
                self._data["somatic_mutation"] = df # Maps dataframe name to dataframe

            else:
                df = pd.read_csv(file_path, sep="\t", index_col=0)
                df = df.transpose()
                df = df.sort_index()
                self._data[df_name] = df # Maps dataframe name to dataframe

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Separate out clinical, derived_molecular, and experimental_setup dataframes
        all_clinical = self._data["clinical"]
        clinical = all_clinical[[
            'Proteomics_Participant_ID', 'Case_excluded',  'Proteomics_Tumor_Normal',  'Country',
            'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity',
            'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN', 'Clin_Stage_Dist_Mets-cM', 'Path_Stage_Dist_Mets-pM',
            'tumor_Stage-Pathological', 'FIGO_stage', 'LVSI', 'BMI', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site',
            'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm',   'Num_full_term_pregnancies']]
        clinical = clinical.rename(columns={"Proteomics_Participant_ID":"Patient_ID"})
        self._data["clinical"] = clinical

        derived_molecular = all_clinical.drop(['Proteomics_Participant_ID', 'Case_excluded',  'Proteomics_Tumor_Normal',  'Country',
            'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity',
            'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN', 'Clin_Stage_Dist_Mets-cM', 'Path_Stage_Dist_Mets-pM',
            'tumor_Stage-Pathological', 'FIGO_stage', 'LVSI', 'BMI', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site',
            'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm',   'Num_full_term_pregnancies',
            'Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs',
            'Proteomics_Aliquot_ID', 'Proteomics_OCT', 'WXS_normal_sample_type', 'WXS_normal_filename', 'WXS_normal_UUID', 'WXS_tumor_sample_type', 'WXS_tumor_filename',
            'WXS_tumor_UUID', 'WGS_normal_sample_type', 'WGS_normal_UUID', 'WGS_tumor_sample_type', 'WGS_tumor_UUID', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID',
            'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality'], axis=1)
        self._data["derived_molecular"] = derived_molecular

        experimental_setup = all_clinical[['Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs',
            'Proteomics_Aliquot_ID', 'Proteomics_OCT', 'WXS_normal_sample_type', 'WXS_normal_filename', 'WXS_normal_UUID', 'WXS_tumor_sample_type', 'WXS_tumor_filename',
            'WXS_tumor_UUID', 'WGS_normal_sample_type', 'WGS_normal_UUID', 'WGS_tumor_sample_type', 'WGS_tumor_UUID', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID',
            'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality']]
        self._data["experimental_setup"] = experimental_setup

        # Add Sample_ID column to somatic_mutations dataframe and make it the index
        clinical = self._data["clinical"] # We need the Patient_ID column from clinical, to map sample ids to patient ids. The sample ids are the clinical index, and the patient ids are in the Patient_ID column.
        patient_id_col = clinical.loc[clinical["Proteomics_Tumor_Normal"] == "Tumor", "Patient_ID"] # We only want to generate a map for tumor samples, because all the normal samples are from the same patients as the tumor samples, so they have duplicate patient ids.
        patient_id_col.index.name = "Sample_ID" # Label the sample id column (it's currently the index)
        patient_id_map = get_reindex_map(patient_id_col)

        mutations = self._data["somatic_mutation"]
        mutations_reindexed = reindex_dataframe(mutations, patient_id_map, "Sample_ID", keep_old=False)
        if mutations_reindexed is None:
            del self._data["somatic_mutation"]
            print("Error mapping sample ids in somatic_mutation dataframe. At least one Patient_ID did not have corresponding Sample_ID mapped in clinical dataframe. somatic_mutation dataframe not loaded.")
        else:
            self._data["somatic_mutation"] = mutations_reindexed

        # Drop all excluded samples from the dataset. They were excluded due to poor sample quality, etc.
        clinical = self._data["clinical"]
        cases_to_drop = clinical[clinical["Case_excluded"] == "Yes"].index

        for name in self._data.keys(): # Loop over the keys so we can alter the values without any issues
            df = self._data[name]
            df_filtered = df.drop(index=cases_to_drop, errors="ignore")
            self._data[name] = df_filtered

        # Drop Case_excluded column from clinical, now that we've dropped all excluded cases in the dataset.
        clinical = self._data["clinical"]
        clinical_no_case_excluded = clinical.drop(columns=["Case_excluded"])
        self._data["clinical"] = clinical_no_case_excluded

        # Sort CNV dataframe columns alphabetically
        cna = self._data["CNA"]
        cna_sorted = cna.sort_index(axis=1)
        self._data["CNA"] = cna_sorted

        # Fix dataframe names
        rename_dict = { # Keys are old names, values are new names
            "CNA":"CNV",
            "transcriptomics_linear":"transcriptomics",
            "transcriptomics_circular":"circular_RNA",
            "phosphoproteomics_site":"phosphoproteomics",
            "somatic_binary":"somatic_mutation_binary",}
        for old, new in rename_dict.items():
            self._data[new] = self._data[old]
            del self._data[old]

        # Rename indices to "Sample_ID", since that's what they all are.
        for name in self._data.keys(): # Loop over the keys so we can alter the values without any issues
            df_rename_index = self._data[name]
            df_rename_index.index.name = "Sample_ID"
            self._data[name] = df_rename_index

        # Drop name of column axis for all dataframes
        for name in self._data.keys(): # Loop over the keys so we can alter the values without any issues
            df = self._data[name]
            df.columns.name = None
            self._data[name] = df

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message

    # Overload the self._get_sample_status_map function to work with "Proteomics_Tumor_Normal" column instead of default "Sample_Tumor_Normal" column
    def _get_sample_status_map(self):
        """Get a pandas Series from the clinical dataframe, with sample ids as the index, and each sample's status (tumor or normal) as the values."""
        clinical = self.get_clinical()
        raw_map = clinical["Proteomics_Tumor_Normal"] 
        parsed_map = raw_map.where(raw_map == "Tumor", other="Normal") # Replace various types of normal (Adjacent_normal, Myometrium_normal, etc.) with just "Normal"
        parsed_map.name = "Sample_Status"
        return parsed_map