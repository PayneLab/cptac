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

from cptac.cancers.source import Source
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, ReindexMapError

class AwgUcec(Source):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the endometrial dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        self.valid_versions = ["2.1", "2.1.1"] # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.

        self.data_files = {
            "2.1": {
                "acetylproteomics"          : "acetylproteomics.cct.gz", 
                "clinical"                  : "clinical.txt", 
                "CNV"                       : "CNA.cct.gz", 
                "definitions"               : "definitions.txt",
                "miRNA"                     : "miRNA.cct.gz", 
                "phosphoproteomics_gene"    : "phosphoproteomics_gene.cct.gz", 
                "phosphoproteomics"         : "phosphoproteomics_site.cct.gz", 
                "proteomics"                : "proteomics.cct.gz", 
                "somatic_mutation_binary"   : "somatic_binary.cbt.gz", 
                "somatic_mutation"          : "somatic.maf.gz", 
                "circular_RNA"              : "transcriptomics_circular.cct.gz", 
                "transcriptomics"           : "transcriptomics_linear.cct.gz"},
            "2.1.1": {
                "acetylproteomics"          : "acetylproteomics.cct.gz",
                "clinical"                  : "clinical.txt",
                "CNV"                       : "CNA.cct.gz",
                "definitions"               : "definitions.txt",
                "miRNA"                     : "miRNA.cct.gz",
                "phosphoproteomics"         : "phosphoproteomics_site.cct.gz",
                "proteomics"                : "proteomics.cct.gz",
                "somatic_mutation_binary"   : "somatic_binary.cbt.gz",
                "somatic_mutation"          : "somatic.maf.gz",
                "circular_RNA"              : "transcriptomics_circular.cct.gz",
                "transcriptomics"           : "transcriptomics_linear.cct.gz",
                "followup"                  : "UCEC_followup_9_12.xlsx"},
        }

        self.load_functions = {
            'acetylproteomics'        : self.load_acetylproteomics,
            'clinical'                : self.load_clinical,
            'CNV'                     : self.load_CNV,
            'derived_molecular'       : self.load_derived_molecular,
            'experimental_design'     : self.load_experimental_design,
            'miRNA'                   : self.load_miRNA,
            'phosphoproteomics'       : self.load_phosphoproteomics,
            'proteomics'              : self.load_proteomics,
            'somatic_mutation_binary' : self.load_somatic_mutation_binary,
            'somatic_mutation'        : self.load_somatic_mutation,
            'circular_RNA'            : self.load_circular_RNA,
            'transcriptomics'         : self.load_transcriptomics,
            'followup'                : self.load_followup,
        }

        if version == "2.1":
            self.load_functions['phosphoproteomics_gene'] = self.load_phosphoproteomics_gene

        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        super().__init__(cancer_type="ucec", source='awg', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)


    def load_acetylproteomics(self):
        df_type = 'acetylproteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df.index = df.index.str.rsplit('-', n=1, expand=True) # Separate the index into a multiindex where the 1st level is the gene, and 2nd is the site
            df.index = df.index.set_names(["Name", "Site"]) # Properly name the levels
            df = df.sort_index()
            df = df.transpose()

            # Clinical contains information on which cases need to be excluded
            self.load_clinical()
            excluded_cases = self._helper_tables["excluded_cases"]
            df = df.drop(index=excluded_cases, errors="ignore")
            # Change index from sample ids to patient ids
            sample_id_to_patient_id_map = self._helper_tables["sample_id_to_patient_id_map"]
            df = reindex_dataframe(df, sample_id_to_patient_id_map, "Patient_ID", False)

            # save df in self._data
            self.save_df(df_type, df)


    def load_clinical(self):
        df_type = 'clinical'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # Fix for reading error on clinical.txt:
            with open(file_path, "r", errors="ignore") as clinical_file:
                df = pd.read_csv(clinical_file, sep='\t', index_col=0)
            df = df.sort_index()

            # Separate out clinical, derived_molecular, and experimental_design dataframes
            all_clinical = df
            clinical = all_clinical[[
                'Proteomics_Participant_ID', 'Case_excluded',  'Proteomics_Tumor_Normal',  'Country',
                'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity',
                'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN', 'Clin_Stage_Dist_Mets-cM', 'Path_Stage_Dist_Mets-pM',
                'tumor_Stage-Pathological', 'FIGO_stage', 'LVSI', 'BMI', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site',
                'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm',   'Num_full_term_pregnancies']]
            clinical = clinical.rename(columns={"Proteomics_Participant_ID":"Patient_ID"})

            # Create list of samples that were excluded due to poor sample quality, etc.
            # Other load functions will need to drop these cases before being saved
            cases_to_drop = clinical[clinical["Case_excluded"] == "Yes"].index.union(clinical[clinical["Case_excluded"] == "Yes"]["Patient_ID"])
            self._helper_tables["excluded_cases"] = cases_to_drop

            # Drop Case_excluded column from clinical, now that we've saved excluded cases for the other datatypes
            clinical = clinical.drop(columns=["Case_excluded"])

            # Add a Sample_Tumor_Normal column to the clinical dataframe, with just "Tumor" or "Normal" values (unlike the Proteomics_Tumor_Normal column, which gives the different types of normal samples)
            raw_map = clinical["Proteomics_Tumor_Normal"] 
            parsed_map = raw_map.where(raw_map == "Tumor", other="Normal") # Replace various types of normal (Adjacent_normal, Myometrium_normal, etc.) with just "Normal"
            clinical.insert(1, "Sample_Tumor_Normal", parsed_map)

            # Mark the Patient_IDs of the normal samples by appending a ".N" to them
            clinical["Patient_ID"] = clinical["Patient_ID"].where(clinical["Sample_Tumor_Normal"] == "Tumor", other=clinical["Patient_ID"] + ".N")

            # Reindex the dataframe with sample IDs instead of patient IDs
            # Also do this for other datatypes in their load functions, but
            # skip the followup and somatic_mutation dataframes because they're already indexed with Patient_IDs
            sample_id_to_patient_id_map = clinical["Patient_ID"]
            self._helper_tables["sample_id_to_patient_id_map"] = sample_id_to_patient_id_map

            clinical.index.name = "Sample_ID" # So that it's labeled properly when we keep it as a column in the clinical dataframe.
            # reindex_dataframe is a function from dataframe_tools.py
            clinical = reindex_dataframe(clinical, sample_id_to_patient_id_map, "Patient_ID", True)

            # We no longer need the Patient_ID column in the clinical dataframe, because it's in the index. So we'll remove it.
            clinical = clinical.drop(columns="Patient_ID")

            # Save other datatypes in clinical file while we are at it
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
            derived_molecular = derived_molecular.drop(index=cases_to_drop, errors="ignore")
            derived_molecular = reindex_dataframe(derived_molecular, sample_id_to_patient_id_map, "Patient_ID", False)
            self.save_df("derived_molecular", derived_molecular)

            experimental_design = all_clinical[['Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs',
                'Proteomics_Aliquot_ID', 'Proteomics_OCT', 'WXS_normal_sample_type', 'WXS_normal_filename', 'WXS_normal_UUID', 'WXS_tumor_sample_type', 'WXS_tumor_filename',
                'WXS_tumor_UUID', 'WGS_normal_sample_type', 'WGS_normal_UUID', 'WGS_tumor_sample_type', 'WGS_tumor_UUID', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID',
                'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality']]

            experimental_design = experimental_design.drop(index=cases_to_drop, errors="ignore")
            experimental_design = reindex_dataframe(experimental_design, sample_id_to_patient_id_map, "Patient_ID", False)
            self.save_df("experimental_design", experimental_design)

            # finally, save clinical in self._data
            self.save_df(df_type, clinical)


    def load_derived_molecular(self):
        if "derived_molecular" not in self._data:
            # This information is contained in the clinical table
            self.load_clinical()


    def load_experimental_design(self):
        if "experimental_design" not in self._data:
            # This information is contained in the clinical table
            self.load_clinical()


    def load_CNV(self):
        df_type = 'CNV'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df = df.sort_index()
            # Sort CNV dataframe columns alphabetically
            df = df.sort_index(axis=1)

            # Clinical contains information on which cases need to be excluded
            self.load_clinical()
            excluded_cases = self._helper_tables["excluded_cases"]
            df = df.drop(index=excluded_cases, errors="ignore")
            # Change index from sample ids to patient ids
            sample_id_to_patient_id_map = self._helper_tables["sample_id_to_patient_id_map"]
            df = reindex_dataframe(df, sample_id_to_patient_id_map, "Patient_ID", False)

            # save df in self._data
            self.save_df(df_type, df)


    def load_definitions(self):
        # This is almost certainly not how we want to handle this
        # Figure out what to do with definitions in the refactor, probably something like with helper tables
        df_type = 'definitions'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            with open(file_path, "r") as definitions_file:
                for line in definitions_file.readlines():
                    line = line.strip()
                    line = line.split("\t")
                    term = line[0]
                    definition = line[1]
                    self._definitions[term] = definition

            # save df in self._data
            self.save_df(df_type, df)


    def load_miRNA(self):
        df_type = 'miRNA'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df = df.sort_index()

            # Clinical contains information on which cases need to be excluded
            self.load_clinical()
            excluded_cases = self._helper_tables["excluded_cases"]
            df = df.drop(index=excluded_cases, errors="ignore")
            # Change index from sample ids to patient ids
            sample_id_to_patient_id_map = self._helper_tables["sample_id_to_patient_id_map"]
            df = reindex_dataframe(df, sample_id_to_patient_id_map, "Patient_ID", False)

            # save df in self._data
            self.save_df(df_type, df)


    def load_phosphoproteomics_gene(self):
        df_type = 'phosphoproteomics_gene'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df = df.sort_index()

            # Clinical contains information on which cases need to be excluded
            self.load_clinical()
            excluded_cases = self._helper_tables["excluded_cases"]
            df = df.drop(index=excluded_cases, errors="ignore")
            # Change index from sample ids to patient ids
            sample_id_to_patient_id_map = self._helper_tables["sample_id_to_patient_id_map"]
            df = reindex_dataframe(df, sample_id_to_patient_id_map, "Patient_ID", False)

            # save df in self._data
            self.save_df(df_type, df)


    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df.index = df.index.str.rsplit('-', n=1, expand=True) # Separate the index into a multiindex where the 1st level is the gene, and 2nd is the site
            df.index = df.index.set_names(["Name", "Site"]) # Properly name the levels
            df = df.sort_index()
            df = df.transpose()

            # Clinical contains information on which cases need to be excluded
            self.load_clinical()
            excluded_cases = self._helper_tables["excluded_cases"]
            df = df.drop(index=excluded_cases, errors="ignore")
            # Change index from sample ids to patient ids
            sample_id_to_patient_id_map = self._helper_tables["sample_id_to_patient_id_map"]
            df = reindex_dataframe(df, sample_id_to_patient_id_map, "Patient_ID", False)

            # save df in self._data
            self.save_df(df_type, df)


    def load_proteomics(self):
        df_type = 'proteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df = df.sort_index()

            # Clinical contains information on which cases need to be excluded
            self.load_clinical()
            excluded_cases = self._helper_tables["excluded_cases"]
            df = df.drop(index=excluded_cases, errors="ignore")
            # Change index from sample ids to patient ids
            sample_id_to_patient_id_map = self._helper_tables["sample_id_to_patient_id_map"]
            df = reindex_dataframe(df, sample_id_to_patient_id_map, "Patient_ID", False)

            # save df in self._data
            self.save_df(df_type, df)


    def load_somatic_mutation_binary(self):
        df_type = 'somatic_mutation_binary'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df = df.sort_index()

            # Clinical contains information on which cases need to be excluded
            self.load_clinical()
            excluded_cases = self._helper_tables["excluded_cases"]
            df = df.drop(index=excluded_cases, errors="ignore")
            # Change index from sample ids to patient ids
            sample_id_to_patient_id_map = self._helper_tables["sample_id_to_patient_id_map"]
            df = reindex_dataframe(df, sample_id_to_patient_id_map, "Patient_ID", False)

            # save df in self._data
            self.save_df(df_type, df)


    def load_somatic_mutation(self):
        df_type = 'somatic_mutation'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n=1, expand=True) # The first part of the barcode is the patient id, which we need want to make the index
            df["Tumor_Sample_Barcode"] = split_barcode[0]
            df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
            df = df.rename({"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
            df = df.sort_values(by=["Patient_ID", "Gene"])
            df = df.set_index("Patient_ID")

            # Clinical contains information on which cases need to be excluded
            self.load_clinical()
            excluded_cases = self._helper_tables["excluded_cases"]
            df = df.drop(index=excluded_cases, errors="ignore")

            # save df in self._data
            self.save_df(df_type, df)


    def load_circular_RNA(self):
        df_type = 'circular_RNA'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df = df.sort_index()

            # Clinical contains information on which cases need to be excluded
            self.load_clinical()
            excluded_cases = self._helper_tables["excluded_cases"]
            df = df.drop(index=excluded_cases, errors="ignore")
            # Change index from sample ids to patient ids
            sample_id_to_patient_id_map = self._helper_tables["sample_id_to_patient_id_map"]
            df = reindex_dataframe(df, sample_id_to_patient_id_map, "Patient_ID", False)

            # save df in self._data
            self.save_df(df_type, df)


    def load_transcriptomics(self):
        df_type = 'transcriptomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()
            df = df.sort_index()

            # Clinical contains information on which cases need to be excluded
            self.load_clinical()
            excluded_cases = self._helper_tables["excluded_cases"]
            df = df.drop(index=excluded_cases, errors="ignore")
            # Change index from sample ids to patient ids
            sample_id_to_patient_id_map = self._helper_tables["sample_id_to_patient_id_map"]
            df = reindex_dataframe(df, sample_id_to_patient_id_map, "Patient_ID", False)

            # save df in self._data
            self.save_df(df_type, df)


    def load_followup(self):
        df_type = 'followup'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_excel(file_path)

            # Replace redundant values for 'not reported' with NaN
            nan_equivalents = ['Not Reported/ Unknown', 'Reported/ Unknown', 'Not Applicable', 'na', 'unknown',
                'Not Performed', 'Unknown tumor status', 'Unknown', 'Unknown Tumor Status', 'Not specified']

            df = df.replace(nan_equivalents, np.nan)

            # Rename, set, and sort index
            df = df.rename(columns={'Case ID': 'Patient_ID'})
            df = df.set_index("Patient_ID")
            df = df.sort_index()

            # Clinical contains information on which cases need to be excluded
            self.load_clinical()
            excluded_cases = self._helper_tables["excluded_cases"]
            df = df.drop(index=excluded_cases, errors="ignore")

            # save df in self._data
            self.save_df(df_type, df)



# This should technically not be necessary in the future
#             # Get a union of all dataframes' indices, with duplicates removed
#             # Exclude the followup dataframe because it has samples from a different cohort that aren't included anywhere else in the dataset
#             master_index = unionize_indices(self._data, exclude="followup")
#             # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset.
#             clinical = self._data["clinical"]
#             clinical = clinical.reindex(master_index)
#             self._data['clinical'] = clinical

# This all should happen in save_df
#             # Call function from dataframe_tools.py to sort all tables first by sample status, and then by the index
#             self._data = sort_all_rows(self._data)
#             # Call function from dataframe_tools.py to standardize the names of the index and column axes
#             self._data = standardize_axes_names(self._data)


    def how_to_cite(self):
        return super().how_to_cite(cancer_type='endometrial carcinoma (uterine)', pmid=32059776)
