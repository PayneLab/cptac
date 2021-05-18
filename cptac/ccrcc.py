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

class Ccrcc(Dataset):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the ccrcc dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        valid_versions = ["0.0", "0.1", "0.1.1"] # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.

        data_files = {
            "0.0": [
                "6_CPTAC3_CCRCC_Phospho_abundance_gene_protNorm=2_CB_imputed.tsv.gz",
                "6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv.gz",
                "6_CPTAC3_CCRCC_Whole_abundance_protein_pep=unique_protNorm=2_CB.tsv.gz",
                "Clinical Table S1.xlsx",
                "ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf.gz",
                "ccrccMethylGeneLevelByMean.txt.gz",
                "cptac-metadata.xls.gz",
                "kirc_wgs_cnv_gene.csv.gz",
                "RNA_Normal_Tumor_185_samples.tsv.gz",
                "S044_CPTAC_CCRCC_Discovery_Cohort_Clinical_Data_r3_Mar2019.xlsx"],
            "0.1": [
                "6_CPTAC3_CCRCC_Phospho_abundance_gene_protNorm=2_CB_imputed.tsv.gz",
                "6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv.gz",
                "6_CPTAC3_CCRCC_Whole_abundance_protein_pep=unique_protNorm=2_CB.tsv.gz",
                "ccrccMethylGeneLevelByMean.txt.gz",
                "ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf.gz",
                "Clinical Table S1.xlsx",
                "cptac-metadata.xls.gz",
                "kirc_wgs_cnv_gene.csv.gz",
                "RNA_Normal_Tumor_185_samples.tsv.gz",
                "S044_CPTAC_CCRCC_Discovery_Cohort_Clinical_Data_r3_Mar2019.xlsx",
                "Table S7.xlsx"],
            "0.1.1": [
                "6_CPTAC3_CCRCC_Phospho_abundance_gene_protNorm=2_CB_imputed.tsv.gz",
                "6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv.gz",
                "6_CPTAC3_CCRCC_Whole_abundance_protein_pep=unique_protNorm=2_CB.tsv.gz",
                "CCRCC_followup_9_12.xlsx",
                "ccrccMethylGeneLevelByMean.txt.gz",
                "ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf.gz",
                "Clinical Table S1.xlsx",
                "cptac-metadata.xls.gz",
                "kirc_wgs_cnv_gene.csv.gz",
                "RNA_Normal_Tumor_185_samples.tsv.gz",
                "S044_CPTAC_CCRCC_Discovery_Cohort_Clinical_Data_r3_Mar2019.xlsx",
                "Table S7.xlsx"],
        }

        super().__init__(cancer_type="ccrcc", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # We're going to need to drop the samples below from a couple dataframes
        nci_labels = ["NCI7-1", "NCI7-2", "NCI7-3", "NCI7-4", "NCI7-5"]
        nci_dotted_labels = [label.replace("-", ".") for label in nci_labels]
        qc_labels = ["QC1", "QC2", "QC3", "QC4", "QC5", "QC6", "QC7", "QC8"]

        # We have multiple clinical files, so we'll load those into a separate dict, and combine them later
        clinical_dfs = {}

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            if file_name == "6_CPTAC3_CCRCC_Phospho_abundance_gene_protNorm=2_CB_imputed.tsv.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                ref_intensities = df["ReferenceIntensity"] # Copy this out, so we can subtract the reference intensities later
                df = df.drop(columns=["NumberPSM", "Proteins", "ReferenceIntensity"] + nci_dotted_labels + qc_labels)
                df = df.subtract(ref_intensities, axis="index") # Subtract reference intensities from all the values, to get ratios
                df = df.transpose()
                self._data["phosphoproteomics_gene"] = df

            elif file_name == "6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')

                # Drop unlocalized sites
                unlocalized_sites = (df["Index"].str.rsplit("_", n=1, expand=True)[1] == '0') 
                df = df[~unlocalized_sites]

                # Drop unwanted samples
                df = df.drop(columns=nci_labels + qc_labels)

                # Subtract reference intensities from numerical data columns, to get ratios
                df = df.rename(columns={"Gene": "Name"})
                metadata_cols = ["Index", "Name", "Peptide", "ReferenceIntensity"]
                metadata = df[metadata_cols] # Extract these for later
                df = df.drop(columns=metadata_cols) # Get the df to contain just the numerical data columns
                ref_intensities = metadata["ReferenceIntensity"]
                df = df.subtract(ref_intensities, axis="index") # Subtract reference intensities
                df = metadata.join(df, how="outer") # Put the metadata columns back in
                df = df.drop(columns="ReferenceIntensity") # Don't need this anymore

                # Parse a few columns out of the "Index" column that we'll need for our multiindex
                split_ids = df["Index"].str.split('_', expand=True)
                df = df.drop(columns="Index")
                sites = split_ids.iloc[:, -1]
                database_ids = split_ids[0].str.cat(split_ids[1], sep='_')
                df = df.assign(**{"Site": sites, "Database_ID": database_ids})

                # Some rows have at least one localized phosphorylation site, but also have other phosphorylations that aren't localized. We'll drop those rows, if their localized sites are duplicated in another row, to avoid creating duplicates, because we only preserve information about the localized sites in a given row. However, if the localized sites aren't duplicated in another row, we'll keep the row.
                unlocalized_to_drop = df.index[~split_ids[4].eq(split_ids[5]) & df.duplicated(["Name", "Site", "Peptide", "Database_ID"], keep=False)] # Column 4 of the split "Index" column is number of phosphorylations detected, and column 5 is number of phosphorylations localized, so if the two values aren't equal, the row has at least one unlocalized site
                df = df.drop(index=unlocalized_to_drop)

                # Give it a multiindex
                df = df.set_index(["Name", "Site", "Peptide", "Database_ID"]) # This will create a multiindex from these columns, in this order.
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                self._data["phosphoproteomics"] = df
            
            elif file_name == "6_CPTAC3_CCRCC_Whole_abundance_protein_pep=unique_protNorm=2_CB.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')
                df = df.rename(columns={"Proteins": "Name", "Index": "Database_ID"})
                df = df.set_index(["Name", "Database_ID"])
                ref_intensities = df["ReferenceIntensity"] # Copy this out, so we can subtract the reference intensities later
                df = df.drop(columns=["NumberPSM", "ReferenceIntensity"] + nci_labels + qc_labels)
                df = df.subtract(ref_intensities, axis="index") # Subtract reference intensities from all the values, to get ratios
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                self._data["proteomics"] = df
            
            elif file_name == "Clinical Table S1.xlsx":
                df = pd.\
                    read_excel(file_path, sheet_name="ccrcc_clinical_characteristics", index_col=0).\
                    dropna(how="all", axis=0).\
                    dropna(how="all", axis=1)
                df.index.name = "Patient_ID" # The index is currently "case_id", but we call that "Patient_ID"
                clinical_dfs["authoritative_clinical"] = df
            
            elif file_name == "ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf.gz":
                df = pd.read_csv(file_path, sep='\t', dtype={"PUBMED":object}) # "PUBMED" column has mixed types, so we specify object as the dtype to prevent a warning from printing. We don't actually use the column, so that's all we need to do.
                split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n=1, expand=True) # The first part of the barcode is the patient id, which we need want to make the index
                df["Tumor_Sample_Barcode"] = split_barcode[0]
                df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
                df = df.rename(columns={"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"})                
                df = df.sort_values(by=["Patient_ID", "Gene"])
                df = df.set_index("Patient_ID")
                self._data["somatic_mutation"] = df

            elif file_name == "ccrccMethylGeneLevelByMean.txt.gz":
                df = pd.read_csv(file_path, sep='\t')
                df = df.sort_index()
                df = df.transpose()
                df.index.name = "Patient_ID"
                self._data["methylation"] = df

            elif file_name == "cptac-metadata.xls.gz":
                df = pd.read_csv(file_path, index_col=0)
                df = df.drop(index=["pooled sample"] + nci_labels + qc_labels) # Drop the pooled samples in addition to other samples to exclude
                clinical_dfs["metadata_and_keys"] = df
            
            elif file_name == "kirc_wgs_cnv_gene.csv.gz":
                df = pd.read_csv(file_path)
                df = df.rename(columns={"gene_name": "Name", "gene_id": "Database_ID"})
                df = df.set_index(["Name", "Database_ID"])
                df = df.sort_index()
                df = df.transpose()

                # The dataframe contains 4 rows for each sample: lr.loc_***,  lr.seg_***, mzd.loc_***, or  mzd.seg_*** where *** is the 
                # Patient_ID. "lr" stands for log ratio, and "mzd" stands for Mean Zygosity Deviation. “lr.loc” is based on average of lr 
                # of probes belonging to each gene. “lr.seg” is based on segmented CNV result (i.e. first performed segmentation based on 
                # probe level data, and then use the lr of the representative segment of each gene as the gene-level lr). They used the
                # lr.seg values in the paper, so we'll use those.
                df = df.drop(index=df[~df.index.str.startswith("lr.seg")].index)

                # Parse a Patient_ID index out of the current index
                barcode_col = df.index.to_series()
                split_barcode = barcode_col.str.split("_", n=1, expand=True) # The second part of the barcode is the patient id, which we want to make the index
                df.index = pd.Index(split_barcode[1])
                df.index.name = "Patient_ID"

                df = df.sort_index()
                self._data["CNV"] = df

            elif file_name == "RNA_Normal_Tumor_185_samples.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')
                df = df.sort_index()
                df = df.transpose()

                # There are a couple duplicate column headers, but they're full of just zeros. We'll drop them.
                # You can do this with a one liner: df =  df.loc[:, ~df.columns.duplicated(keep=False) | (df != 0).any(axis=0)]
                # But the one liner is about 100 times slower.
                dups = df.loc[:, df.columns.duplicated(keep=False)] # Select all the columns with duplicated headers
                dups = dups.loc[:, (dups != 0).any(axis=0)] # Get only the columns that aren't all zeros
                df = df.loc[:, ~df.columns.duplicated(keep=False)] # Get rid of the duplicate columns from the original dataframe
                df = df.join(dups, how="outer") # Sub in our un-duplicated selections
                df = df.sort_index(axis="columns") # Get all the column names in order again

                self._data["transcriptomics"] = df
            
            elif file_name == "S044_CPTAC_CCRCC_Discovery_Cohort_Clinical_Data_r3_Mar2019.xlsx":

                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", "Unknown extension is not supported and will be removed") # This warning is just due to some formatting on the spreadsheet. We don't need to worry about it.
                    clinical_sheets = pd.read_excel(file_path,
                        sheet_name=['Patient_Clinical_Attributes', 'Other_Medical_Information', 'Specimen_Attributes'],
                        index_col=0,
                        usecols=lambda x: x != "tumor_code") # Don't load the tumor_code column in any of them--it's just "CCRCC" for every row.

                for sheet, df in clinical_sheets.items(): # Put them in the clinical_dfs dictionary for processing later
                    df.index.name = "Patient_ID" # The indices are currently "case_id", but we call that "Patient_ID"
                    clinical_dfs[sheet] = df

            elif file_name == "Table S7.xlsx":
                immune_groups = pd.\
                    read_excel(
                        file_path, 
                        sheet_name="xCell Signatures", 
                        index_col = 0, 
                        header=None,
                        skiprows=range(0, 3)
                    ).\
                    transpose()
                immune_groups = immune_groups[["Samples", "Immune Group"]] # We only need these columns
                immune_groups = immune_groups.set_index("Samples")

            elif file_name == 'CCRCC_followup_9_12.xlsx' and self._version == "0.1.1":
                df = pd.read_excel(file_path)

                # Replace redundant values for "not reported" with NaN
                nan_equivalents = ['Not Reported/ Unknown', 'Reported/ Unknown', 'Not Applicable', 'na',
                    'unknown', 'Not Performed', 'Unknown tumor status', 'Unknown', ' Unknown', 'Unknown ',
                    'Unknown Tumor Status', 'Not specified']

                df = df.replace(nan_equivalents, np.nan)

                # Replace redundanct values in "Cause of Death" column
                disease_prog_equivalents = ['Metastatic Renal Cell Carcinoma', 'Tumor progression',
                    'Progression of disease', 'Progression of disease ', 'Tumor', 'Disease progression',
                    'Progressive Disease', 'Disease progression', 'disease progression ', 'main disease ']

                df['Cause of Death'] = df['Cause of Death'].replace(disease_prog_equivalents, 'Disease progression')

                # Rename, set, and sort by index
                df = df.rename(columns={"Case ID": "Patient_ID"})
                df = df.set_index("Patient_ID")
                df = df.sort_index()

                self._data["followup"] = df

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Process and combine the multiple clinical dataframes
        clinical = clinical_dfs["metadata_and_keys"] # We'll start with this dataframe, and add the others to it.
        clinical.index = clinical.index.where(~(clinical["Type"] == "Normal"), 'N' + clinical.index) # Prepend an 'N' to the patient IDs of normal samples
        clinical.index.name = "Patient_ID" # Name the index Patient_ID, since the index values are in fact the patient IDs (also called case IDs)
        clinical = clinical.sort_index()
        clinical = clinical.rename(columns={"Type": "Sample_Tumor_Normal"}) # This matches the other datasets.
        clinical = clinical.join(clinical_dfs["Patient_Clinical_Attributes"], how="outer")

        auth_clin = clinical_dfs["authoritative_clinical"] # This table has more authoritative values for histologic type, so we'll replace that column
        clinical["histologic_type"] = auth_clin["Histologic_Type"]

        specimen_attributes = clinical_dfs["Specimen_Attributes"] # Before we can join in this dataframe, we need to prepend "N" to the normal samples' index values
        specimen_attributes.index = np.where(specimen_attributes["tissue_type"] == "normal", 'N' + specimen_attributes.index, specimen_attributes.index)
        specimen_attributes.index.name = "Patient_ID" # Name the index Patient_ID, since the index values are in fact the patient IDs (also called case IDs)
        clinical = clinical.join(specimen_attributes, how="outer") # Join in our nicely reindexed dataframe

        medical_info = clinical_dfs["Other_Medical_Information"] # We're going to join one column from this to the clinical dataframe, and later take some other columns to make a medical_history dataframe
        medication_col = medical_info["medication_name"]
        medication_col = medication_col.str.replace("|", ",", regex=False)
        clinical = clinical.assign(patient_medications=medication_col) # Add it in as a new "patient_medications" column
        clinical = clinical.drop(columns=["MS.Directory.Name", "Batch", "Data.Set", "TMT.Channel", "Mass.Spectrometer", "Mass.Spectrometer.Operator", "Set.A", "Set.B", "tissue_type"]) # tissue_type column is a duplicate of Sample_Tumor_Normal

        # If this if version 0.1 or greater, join the immune group data into the clinical dataframe
        if self._version == "0.1":
            clinical = clinical.join(immune_groups, on="Specimen.Label", how="outer")

        # Save this final amalgamation of all the clinical dataframes
        self._data["clinical"] = clinical 

        # Create the medical_history dataframe
        medical_info = clinical_dfs["Other_Medical_Information"]
        medical_history_unparsed = medical_info[["medical_condition", "medical_condition_other", "condition_year_of_onset", "history_of_medical_treatment", "medical_history_source"]]
        medical_history_unparsed = medical_history_unparsed.dropna(axis="index", how="all")

        first_time = True # We'll use this to get special behavior on our first iteration
        for col_name in medical_history_unparsed.columns:
            med_hist_col = medical_history_unparsed[col_name].str.split("|", expand=True) # Cells have mutliple values separated by "|", so we'll splits those into individual columns
            med_hist_col = med_hist_col.reset_index() # Make the sample ids a column so we can use them in melt
            med_hist_col = med_hist_col.melt(id_vars="Patient_ID") # Use each row's patient id as its index, then combine everything else into one column
            if first_time:
                med_hist_col = med_hist_col.drop(columns="variable") # Had the old column headers
            else:
                med_hist_col = med_hist_col.drop(columns=["Patient_ID", "variable"]) # Also dropping patient ids so we're not trying to join dfs with duplicate columns at the end of the loop
            med_hist_col = med_hist_col.rename(columns={"value": col_name})
            med_hist_col = med_hist_col.dropna()
            if first_time:
                medical_history_parsed = med_hist_col 
            else:
                medical_history_parsed = medical_history_parsed.join(med_hist_col, how="outer")
            first_time = False

        med_cond = medical_history_parsed["medical_condition"] # Where medical_condition column is "Other, specify", we're going to fill in value from medical_condition_other column
        med_cond = med_cond.where(med_cond != "Other, specify", medical_history_parsed["medical_condition_other"])
        medical_history_parsed["medical_condition"] = med_cond
        medical_history_parsed = medical_history_parsed.drop(columns="medical_condition_other") # Since we copied all of the useful information into the medical_condition column

        medical_history_parsed = medical_history_parsed.set_index("Patient_ID")
        medical_history_parsed = medical_history_parsed.sort_index()
        self._data["medical_history"] = medical_history_parsed

        # Use the RNA.ID column from clinical dataframe to reindex transcriptomics dataframe with patient ids
        tran = self._data["transcriptomics"]
        try:
            tran_map = get_reindex_map(clinical["RNA.ID"])
            tran_reindexed = reindex_dataframe(tran, tran_map, new_index_name="Patient_ID", keep_old=False)
        except ReindexMapError:
            del self._data["transcriptomics"]
            warnings.warn("Error mapping sample ids in transcriptomics dataframe. At least one RNA.ID did not have a corresponding Patient_ID mapped in the clinical dataframe. transcriptomics dataframe not loaded.", FailedReindexWarning, stacklevel=2)
        else:
            self._data["transcriptomics"] = tran_reindexed

        # Use the Specimen.Label columns from clinical dataframe to reindex the proteomics, phosphoproteomics, and phosphoproteomics_gene dataframes with patient ids
        specimen_indexed_dfs = ["proteomics", "phosphoproteomics", "phosphoproteomics_gene"]
        try:
            specimen_label_map = get_reindex_map(clinical["Specimen.Label"])
        except ReindexMapError:
            for df_name in specimen_indexed_dfs:
                del self._data[df_name]
            warnings.warn(f"Error mapping sample ids in these dataframes: {' '.join(df for df in specimen_indexed_dfs)}. Specimen.Label mapping in clinical dataframe was not one-to-one. Dataframes not loaded.", FailedReindexWarning, stacklevel=2)
        else:
            for df_name in specimen_indexed_dfs:
                df = self._data[df_name]
                try:
                    df_reindexed = reindex_dataframe(df, specimen_label_map, new_index_name="Patient_ID", keep_old=False)
                except ReindexMapError as error:
                    del self._data[df_name]
                    warnings.warn(f"Error mapping sample ids in {df_name} dataframe. RNA.ID {str(error)} did not have a corresponding Patient_ID mapped in the clinical dataframe. {df_name} dataframe not loaded.", FailedReindexWarning, stacklevel=2)
                else:
                    self._data[df_name] = df_reindexed

        # Now that we've used the RNA.ID and Specimen.Label columns to reindex the dataframes that needed it, we can drop them from the clinical dataframe
        clinical = clinical.drop(columns=["Specimen.Label", "RNA.ID"])

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = clinical

        # There are 28 samples in the CNV table, 26 of which are also in the methylation table, that aren't in any other tables, and were actually excluded from analysis. We're going to drop them.
        # All samples used for analysis should be in the proteomics dataframe, so we'll use that test to determine which samples to drop from CNV and methylation.
        prot = self._data["proteomics"]

        cnv = self._data["CNV"]
        to_drop = cnv.index.difference(prot.index)
        cnv = cnv.drop(index=to_drop)
        self._data["CNV"] = cnv

        methylation = self._data["methylation"]
        to_drop = methylation.index.difference(prot.index)
        methylation = methylation.drop(index=to_drop)
        self._data["methylation"] = methylation

        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data, exclude="followup")

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
        clinical = self._data["clinical"]
        clinical = clinical.reindex(master_index)
        self._data['clinical'] = clinical

        # Edit the format of the Patient_IDs to have normal samples marked the same way as in other datasets. Currently, normal patient IDs have an "N" prepended. We're going to erase that and append a ".N"
        self._data = reformat_normal_patient_ids(self._data, existing_identifier="N", existing_identifier_location="start")

        # Call function from dataframe_tools.py to sort all tables first by sample status, and then by the index
        self._data = sort_all_rows(self._data)

        # Call function from dataframe_tools.py to standardize the names of the index and column axes
        self._data = standardize_axes_names(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message

    def how_to_cite(self):
        return super().how_to_cite(cancer_type='clear cell renal cell carcinoma (kidney)', pmid=31675502)