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
import warnings
from .dataset import Dataset
from .dataframe_tools import *
from .exceptions import FailedReindexWarning, ReindexMapError

class Endometrial(Dataset):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the endometrial dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        valid_versions = ["2.1", "2.1.1"] # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.

        data_files = {
            "2.1": [
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
                "transcriptomics_linear.cct.gz"],
            "2.1.1": [
                "acetylproteomics.cct.gz",
                "clinical.txt",
                "CNA.cct.gz",
                "definitions.txt",
                "miRNA.cct.gz",
                "phosphoproteomics_site.cct.gz",
                "proteomics.cct.gz",
                "somatic_binary.cbt.gz",
                "somatic.maf.gz",
                "transcriptomics_circular.cct.gz",
                "transcriptomics_linear.cct.gz",
                "UCEC_followup_9_12.xlsx"],
        }

        super().__init__(cancer_type="endometrial", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data files into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: 

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
                split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n=1, expand=True) # The first part of the barcode is the patient id, which we need want to make the index
                df["Tumor_Sample_Barcode"] = split_barcode[0]
                df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
                df = df.rename({"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
                df = df.sort_values(by=["Patient_ID", "Gene"])
                df = df.set_index("Patient_ID")
                self._data["somatic_mutation"] = df # Maps dataframe name to dataframe

            elif file_name == "acetylproteomics.cct.gz" or file_name == "phosphoproteomics_site.cct.gz":
                df = pd.read_csv(file_path, sep = "\t", index_col=0)
                df.index = df.index.str.rsplit('-', n=1, expand=True) # Separate the index into a multiindex where the 1st level is the gene, and 2nd is the site
                df.index = df.index.set_names(["Name", "Site"]) # Properly name the levels
                df = df.sort_index()
                df = df.transpose()
                self._data[df_name] = df # Maps dataframe name to dataframe

            elif file_name == 'UCEC_followup_9_12.xlsx' and self._version == "2.1.1":
                df = pd.read_excel(file_path)

                # Replace redundant values for 'not reported' with NaN
                nan_equivalents = ['Not Reported/ Unknown', 'Reported/ Unknown', 'Not Applicable', 'na', 'unknown',
                    'Not Performed', 'Unknown tumor status', 'Unknown', 'Unknown Tumor Status', 'Not specified']
                    
                df = df.replace(nan_equivalents, np.nan)

                # Rename, set, and sort index
                df = df.rename(columns={'Case ID': 'Patient_ID'})
                df = df.set_index("Patient_ID")
                df = df.sort_index()

                self._data["followup"] = df

            else:
                df = pd.read_csv(file_path, sep="\t", index_col=0)
                df = df.transpose()
                df = df.sort_index()
                self._data[df_name] = df # Maps dataframe name to dataframe

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Separate out clinical, derived_molecular, and experimental_design dataframes
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
        derived_molecular = derived_molecular.rename(columns={"JAK1_Mutation":"JAK1_Mutation_status"})
        self._data["derived_molecular"] = derived_molecular

        experimental_design = all_clinical[['Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs',
            'Proteomics_Aliquot_ID', 'Proteomics_OCT', 'WXS_normal_sample_type', 'WXS_normal_filename', 'WXS_normal_UUID', 'WXS_tumor_sample_type', 'WXS_tumor_filename',
            'WXS_tumor_UUID', 'WGS_normal_sample_type', 'WGS_normal_UUID', 'WGS_tumor_sample_type', 'WGS_tumor_UUID', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID',
            'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality']]
        self._data["experimental_design"] = experimental_design

        # Drop all excluded samples from the dataset. They were excluded due to poor sample quality, etc.
        clinical = self._data["clinical"]
        cases_to_drop = clinical[clinical["Case_excluded"] == "Yes"].index.union(clinical[clinical["Case_excluded"] == "Yes"]["Patient_ID"])

        for name in self._data.keys(): # Loop over the keys so we can alter the values without any issues
            df = self._data[name]
            df_filtered = df.drop(index=cases_to_drop, errors="ignore")
            self._data[name] = df_filtered

        # Drop Case_excluded column from clinical, now that we've dropped all excluded cases in the dataset.
        clinical = self._data["clinical"]
        clinical = clinical.drop(columns=["Case_excluded"])

        # Add a Sample_Tumor_Normal column to the clinical dataframe, with just "Tumor" or "Normal" values (unlike the Proteomics_Tumor_Normal column, which gives the different types of normal samples)
        raw_map = clinical["Proteomics_Tumor_Normal"] 
        parsed_map = raw_map.where(raw_map == "Tumor", other="Normal") # Replace various types of normal (Adjacent_normal, Myometrium_normal, etc.) with just "Normal"
        clinical.insert(1, "Sample_Tumor_Normal", parsed_map)

        # Mark the Patient_IDs of the normal samples by appending a ".N" to them
        clinical["Patient_ID"] = clinical["Patient_ID"].where(clinical["Sample_Tumor_Normal"] == "Tumor", other=clinical["Patient_ID"] + ".N")

        # Save our new and improved clinical dataframe!
        self._data["clinical"] = clinical

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
            "somatic_binary":"somatic_mutation_binary",
            }
        for old, new in rename_dict.items():
            self._data[new] = self._data[old]
            del self._data[old]

        # Call a function from dataframe_tools.py to reindex all the dataframes with sample IDs instead of patient IDs
        # Skip the followup and somatic_mutation dataframes because they're already indexed with Patient_IDs
        sample_id_to_patient_id_map = self._data["clinical"]["Patient_ID"]
        self._data = reindex_all_sample_id_to_patient_id(self._data, sample_id_to_patient_id_map, skip=["followup", "somatic_mutation"])

        # We no longer need the Patient_ID column in the clinical dataframe, because it's in the index. So we'll remove it.
        clinical = self._data["clinical"]
        clinical = clinical.drop(columns="Patient_ID")
        self._data["clinical"] = clinical

        # Get a union of all dataframes' indices, with duplicates removed
        # Exclude the followup dataframe because it has samples from a different cohort that aren't included anywhere else in the dataset
        master_index = unionize_indices(self._data, exclude="followup")

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset.
        clinical = self._data["clinical"]
        clinical = clinical.reindex(master_index)
        self._data['clinical'] = clinical

        # Call function from dataframe_tools.py to sort all tables first by sample status, and then by the index
        self._data = sort_all_rows(self._data)

        # Call function from dataframe_tools.py to standardize the names of the index and column axes
        self._data = standardize_axes_names(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message

    def how_to_cite(self):
        return super().how_to_cite(cancer_type='endometrial carcinoma (uterine)', pmid=32059776)
