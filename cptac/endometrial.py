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
import glob
import textwrap
import datetime
from .dataset import DataSet
from .fileloader import check_data

class Endometrial(DataSet):

    def __init__(self):
        """Load all of the endometrial dataframes as values in the self._data dict variable, with names as keys, and format them properly."""

        # Call the parent DataSet __init__ function, which initializes self._data and other variables we need
        super().__init__()

        # Check the data files. If they're not downloaded, download them. If they're out of date, update them.
        path_here = os.path.abspath(os.path.dirname(__file__))
        data_directory = os.path.join(path_here, "data_endometrial")
        check_data(data_directory)

        # Print the data version
        data_version = "2.1"
        print("endometrial data version: {}\n".format(data_version))

        # Get the path to the data files
        all_data_path = os.path.join(data_directory, "*.*")
        files = glob.glob(all_data_path) # Put all files into a list
        files = [file for file in files if not os.path.isdir(file)] # Take out the urls directory
        files = sorted(files, key=str.lower)

        # Load the data files into dataframes in the self._data dict
        for file in files: 
            path_elements = file.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            file_name_split = file_name.split(".")
            df_name = file_name_split[0] # Dataframe name will be the first section of file name; i.e. proteomics.txt.gz becomes proteomics

            # Load the file, based on what it is
            print("Loading {} data...".format(df_name))
            if file_name == "clinical.txt":
                # Fix for reading error on clinical.txt:
                with open(file, "r", errors="ignore") as clinical_file:
                    df = pd.read_csv(clinical_file, sep="\t", index_col=0)
                df = df.sort_index()
                df.name = df_name
                self._data[df.name] = df # Maps dataframe name to dataframe
            elif file_name == "somatic.maf.gz":
                df = pd.read_csv(file, sep = "\t")
                split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n = 1, expand = True) # The first part of the barcode is the patient id, which we need want to make the index
                df["Tumor_Sample_Barcode"] = split_barcode[0]
                df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
                df = df.rename({"Tumor_Sample_Barcode":"Patient_Id","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
                df.name = "somatic_mutation"
                self._data[df.name] = df # Maps dataframe name to dataframe
            elif file_name in ("acetylproteomics.cct.gz", "CNA.cct.gz", "miRNA.cct.gz", "phosphoproteomics_gene.cct.gz", "phosphoproteomics_site.cct.gz", "proteomics.cct.gz", "somatic_binary.cbt.gz", "transcriptomics_circular.cct.gz", "transcriptomics_linear.cct.gz"):
                df = pd.read_csv(file, sep="\t", index_col=0)
                df = df.transpose()
                df = df.sort_index()
                df.name = df_name
                self._data[df.name] = df # Maps dataframe name to dataframe
            elif file_name == "definitions.txt":
                pass # We'll load the defintions separately
            else:
                print("Unrecognized file: {}.\nFile not loaded.".format(file))

        # Separate out clinical, derived_molecular, and experimental_setup dataframes
        all_clinical = self._data["clinical"]
        clinical = all_clinical[[
            'Proteomics_Participant_ID', 'Case_excluded',  'Proteomics_Tumor_Normal',  'Country',
            'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity',
            'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN', 'Clin_Stage_Dist_Mets-cM', 'Path_Stage_Dist_Mets-pM',
            'tumor_Stage-Pathological', 'FIGO_stage', 'LVSI', 'BMI', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site',
            'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm',   'Num_full_term_pregnancies']]
        clinical = clinical.rename(columns={"Proteomics_Participant_ID":"Patient_ID"})
        clinical.name = "clinical"
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
        derived_molecular.name = "derived_molecular"
        self._data["derived_molecular"] = derived_molecular

        experimental_setup = all_clinical[['Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs',
            'Proteomics_Aliquot_ID', 'Proteomics_OCT', 'WXS_normal_sample_type', 'WXS_normal_filename', 'WXS_normal_UUID', 'WXS_tumor_sample_type', 'WXS_tumor_filename',
            'WXS_tumor_UUID', 'WGS_normal_sample_type', 'WGS_normal_UUID', 'WGS_tumor_sample_type', 'WGS_tumor_UUID', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID',
            'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality']]
        experimental_setup.name = "experimental_setup"
        self._data["experimental_setup"] = experimental_setup

        # Add Sample_ID column to somatic_mutations dataframe and make it the index
        clinical = self._data["clinical"] # We need the Patient_ID column from clinical, to map sample ids to patient ids. The sample ids are the clinical index, and the patient ids are in the Patient_ID column.
        patient_id_col = clinical.loc[clinical["Proteomics_Tumor_Normal"] == "Tumor", "Patient_ID"] # We only want to generate a map for tumor samples, because all the normal samples are from the same patients as the tumor samples, so they have duplicate patient ids.
        patient_id_col.index.name = "Sample_ID" # Label the sample id column (it's currently the index)
        patient_id_df = patient_id_col.reset_index() # Make the sample id index accessible as a column
        patient_id_df = patient_id_df.set_index("Patient_ID") # Set Patient_ID as the index, so we can look up a Sample_ID given a Patient_ID
        patient_id_map = patient_id_df["Sample_ID"] # Make the mapping a series. Patient_ID will be the index.

        mutations = self._data["somatic_mutation"]
        mutations_patient_renamed = mutations.rename(columns={"Patient_Id":"Patient_ID"})
        mutations_patient_indexed = mutations_patient_renamed.set_index("Patient_ID") # Set the index as the Patient_ID column, dropping the default numerical index
        sample_id_col = [] # We're going to create a Sample_ID column for the mutations dataframe
        map_success = True
        for patient_id in mutations_patient_indexed.index:
            if patient_id in patient_id_map.index:
                sample_id_col.append(patient_id_map[patient_id]) # Get the sample id corresponding to the patient id
            else: # If there's not a corresponding sample ID for a patient ID, print an error message and return None
                print("Error mapping sample ids in somatic_mutation dataframe. Patient_ID {} did not have corresponding Sample_ID mapped in clinical dataframe. somatic_mutation dataframe not loaded.".format(patient_id))
                map_success = False
        if map_success:
            mutations_with_sample = mutations_patient_indexed.assign(Sample_ID=sample_id_col) # Add in the Sample_ID column
            mutations_sample_indexed = mutations_with_sample.set_index("Sample_ID") # Make the Sample_ID column the index
            mutations_sample_indexed.name = mutations.name
            self._data["somatic_mutation"] = mutations_sample_indexed
        else:
            del self._data["somatic_mutation"]

        # Drop all excluded samples from the dataset. They were excluded due to poor sample quality, etc.
        clinical = self._data["clinical"]
        cases_to_drop = clinical[clinical["Case_excluded"] == "Yes"].index

        for name in self._data.keys(): # Loop over the keys instead of directly over the dict, so we're not altering the structure we're looping over
            df = self._data[name]
            df_filtered = df.drop(index=cases_to_drop, errors="ignore")
            df_filtered.name = df.name
            self._data[name] = df_filtered

        # Drop Case_excluded column from clinical, now that we've dropped all excluded cases in the dataset.
        clinical = self._data["clinical"]
        clinical_no_case_excluded = clinical.drop(columns=["Case_excluded"])
        clinical_no_case_excluded.name = clinical.name
        self._data["clinical"] = clinical_no_case_excluded

        # Sort CNA dataframe columns alphabetically
        cna = self._data["CNA"]
        cna_sorted = cna.sort_index(axis=1)
        cna_sorted.name = cna.name
        self._data["CNA"] = cna_sorted

        # Fix dataframe names
        rename_dict = { # Keys are old names, values are new names
            "transcriptomics_linear":"transcriptomics",
            "transcriptomics_circular":"circular_RNA",
            "phosphoproteomics_site":"phosphoproteomics",
            "somatic_binary":"somatic_mutation_binary",}
        for old, new in rename_dict.items():
            rename_df = self._data[old]
            rename_df.name = new
            self._data[new] = rename_df
            del self._data[old]

        # Rename indicies to "Sample_ID", since that's what they all are.
        for name in self._data.keys():
            df_rename_index = self._data[name]
            df_rename_index.index.name = "Sample_ID"
            self._data[name] = df_rename_index

        # Drop name of column axis for all dataframes
        for name in self._data.keys():
            df_rename_col_axis = self._data[name]
            df_rename_col_axis.columns.name = None
            self._data[name] = df_rename_col_axis

        # Load definitions
        definitions_path = os.path.join(data_directory, "definitions.txt")
        with open(definitions_path, "r") as definitions_file:
            for line in definitions_file.readlines():
                line = line.strip()
                line = line.split("\t")
                self._definitions[line[0]] = line[1]

        # Print data embargo warning, if the date hasn't passed yet.
        today = datetime.date.today()
        embargo_date = datetime.date(2019, 7, 1)
        if today < embargo_date:
            print("\n","******PLEASE READ******")
            warning = "WARNING: This data is under a publication embargo until July 1, 2019. CPTAC is a community resource project and data are made available rapidly after generation for community research use. The embargo allows exploring and utilizing the data, but analysis may not be published until July 1, 2019. Please see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter cptac.embargo() to open the webpage for more details."
            wrapped_list = textwrap.wrap(warning)
            for line in wrapped_list:
                print(line)

    # Overload the self._get_sample_status_map function to work with "Proteomics_Tumor_Normal" column instead of default "Sample_Tumor_Normal" column
    def _get_sample_status_map(self):
        """Get a pandas Series from the clinical dataframe, with sample ids as the index, and each sample's status (tumor or normal) as the values."""
        clinical = self.get_clinical()
        raw_map = clinical["Proteomics_Tumor_Normal"] 
        parsed_map = raw_map.where(raw_map == "Tumor", other="Normal") # Replace various types of normal (Adjacent_normal, Myometrium_normal, etc.) with just "Normal"
        parsed_map.name = "Sample_Status"
        return parsed_map
