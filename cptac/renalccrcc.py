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

class RenalCcrcc(DataSet):

    def __init__(self, version="latest"):
        """Load all of the renalccrcc dataframes as values in the self._data dict variable, with names as keys, and format them properly."""

        # Call the parent DataSet __init__ function, which initializes self._data and other variables we need
        super().__init__("renalccrcc")

        # Update the index, if possible. If there's no internet, update_index will return False, but we don't care in this context.
        update_index(self._cancer_type)

        # Validate the index
        self._version = validate_version(version, self._cancer_type, use_context="init")
        if self._version is None: # Validation error. validate_version already printed an error message.
            return

        # Get the paths to all the data files
        data_files = [
            "6_CPTAC3_CCRCC_Phospho_abundance_gene_protNorm=2_CB_imputed.tsv.gz",
            "6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv.gz",
            "6_CPTAC3_CCRCC_Whole_abundance_protein_pep=unique_protNorm=2_CB.tsv.gz",
            "Clinical Table S1.xlsx",
            "ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf.gz",
            "ccrccMethylGeneLevelByMean.txt.gz",
            "cptac-metadata.xls.gz",
            "kirc_wgs_cnv_gene.csv.gz",
            "RNA_Normal_Tumor_185_samples.tsv.gz",
            "S044_CPTAC_CCRCC_Discovery_Cohort_Clinical_Data_r3_Mar2019.xlsx",]
        data_files_paths = get_version_files_paths(self._cancer_type, self._version, data_files)
        if data_files_paths is None: # Data error. get_version_files_paths already printed an error message.
            return

        # We're going to need to drop the samples below from a couple dataframes
        nci_labels = ["NCI7-1", "NCI7-2", "NCI7-3", "NCI7-4", "NCI7-5"]
        nci_dotted_labels = [label.replace("-", ".") for label in nci_labels]
        qc_labels = ["QC1", "QC2", "QC3", "QC4", "QC5", "QC6", "QC7", "QC8"]

        # We have multiple clinical files, so we'll load those into a separate dict, and combine them later
        clinical_dfs = {}

        # Load the data into dataframes in the self._data dict
        loading_msg = "Loading dataframes"
        for file_path in data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            if file_name == "6_CPTAC3_CCRCC_Phospho_abundance_gene_protNorm=2_CB_imputed.tsv.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df.index.name = None # It's going to become our column headers, and we don't want the column axes to have names
                ref_intensities = df["ReferenceIntensity"] # Copy this out, so we can subtract the reference intensities later
                df = df.drop(columns=["NumberPSM", "Proteins", "ReferenceIntensity"] + nci_dotted_labels + qc_labels)
                df = df.subtract(ref_intensities, axis="index") # Subtract reference intensities from all the values, to get ratios
                df = df.transpose()
                self._data["phosphoproteomics_gene"] = df

            elif file_name == "6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')

                # Parse out the sites and append to gene names
                old_index = df["Index"]
                split_index = old_index.str.rsplit("_", n=1, expand=True)
                sites = split_index[1]
                genes = df["Gene"]
                genes_with_sites = genes.str.cat(sites, sep='-')
                df["Gene"] = genes_with_sites

                # Sites where the site location is 0 are unlocalized, and we're going to drop them.
                unlocalized_sites = (sites == '0') 
                df = df[~unlocalized_sites]

                df = df.set_index("Gene")
                ref_intensities = df["ReferenceIntensity"] # Copy this out, so we can subtract the reference intensities later
                df = df.drop(columns=["Index", "Peptide", "ReferenceIntensity"] + nci_labels + qc_labels)
                df = df.subtract(ref_intensities, axis="index") # Subtract the reference intensities from all the values, to get ratios
                df = df.sort_index()
                df = df.transpose()
                self._data["phosphoproteomics"] = df
            
            elif file_name == "6_CPTAC3_CCRCC_Whole_abundance_protein_pep=unique_protNorm=2_CB.tsv.gz":
                df = pd.read_csv(file_path, sep='\t')
                df = df.set_index("Proteins")
                df.index.name = None # It's going to become our column headers, and we don't want the column axes to have names
                ref_intensities = df["ReferenceIntensity"] # Copy this out, so we can subtract the reference intensities later
                df = df.drop(columns=["NumberPSM", "Index", "ReferenceIntensity"] + nci_labels + qc_labels)
                df = df.subtract(ref_intensities, axis="index") # Subtract reference intensities from all the values, to get ratios
                df = df.sort_index()
                df = df.transpose()
                self._data["proteomics"] = df
            
            elif file_name == "Clinical Table S1.xlsx":
                df = pd.read_excel(file_path, sheet_name="ccrcc_clinical_characteristics", index_col=0) # This file has multiple sheets, but we only want the one
                df.index.name = "Patient_ID" # The index is currently "case_id", but we call that "Patient_ID"
                clinical_dfs["authoritative_clinical"] = df.copy()
            
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
                df = df.drop(columns="gene_id")
                df = df.set_index("gene_name")
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
                self._data["transcriptomics"] = df
            
            elif file_name == "S044_CPTAC_CCRCC_Discovery_Cohort_Clinical_Data_r3_Mar2019.xlsx":
                clinical_sheets = pd.read_excel(file_path, # This file has multiple sheets, but we only need the ones specified on the next line.
                    sheet_name=['Patient_Clinical_Attributes', 'Other_Medical_Information', 'Specimen_Attributes'],
                    index_col=0,
                    usecols=lambda x: x != "tumor_code") # Don't load the tumor_code column in any of them--it's just "CCRCC" for every row.

                for sheet, df in clinical_sheets.items(): # Put them in the clinical_dfs dictionary for processing later
                    df.index.name = "Patient_ID" # The indices are currently "case_id", but we call that "Patient_ID"
                    clinical_dfs[sheet] = df.copy()

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Process and combine the multiple clinical dataframes
        clinical = clinical_dfs["metadata_and_keys"] # We'll start with this dataframe, and add the others to it.
        new_clinical_index = [] # We gonna reindex this here dataframe
        for index, row in clinical.iterrows(): # Prepend "N" to the index values of normal samples
            if row["Type"] == "Normal":
                index = "N" + index
            new_clinical_index.append(index)

        clinical.index = new_clinical_index
        clinical.index.name = "Patient_ID" # Name the index Patient_ID, since the index values are in fact the patient IDs (also called case IDs)
        clinical = clinical.sort_index()
        clinical = clinical.rename(columns={"Type": "Sample_Tumor_Normal"}) # This matches the other datasets.
        clinical = clinical.join(clinical_dfs["Patient_Clinical_Attributes"], how="outer")

        auth_clin = clinical_dfs["authoritative_clinical"] # This table has more authoritative values for histologic type, so we'll replace that column
        clinical["histologic_type"] = auth_clin["Histologic_Type"]

        specimen_attributes = clinical_dfs["Specimen_Attributes"] # Before we can join in this dataframe, we need to prepend "N" to the normal samples' index values
        new_specimen_attributes_index = []
        for index, row in specimen_attributes.iterrows(): # Prepend "N" to the index values of normal samples
            if row["tissue_type"] == "normal":
                index = "N" + index
            new_specimen_attributes_index.append(index)
        specimen_attributes.index = new_specimen_attributes_index
        specimen_attributes.index.name = "Patient_ID" # Name the index Patient_ID, since the index values are in fact the patient IDs (also called case IDs)
        clinical = clinical.join(specimen_attributes, how="outer") # Join in our nicely reindexed dataframe

        medical_info = clinical_dfs["Other_Medical_Information"] # We're going to join one column from this to the clinical dataframe, and later take some other columns to make a medical_history dataframe
        medication_col = medical_info["medication_name"]
        medication_col = medication_col.str.replace("|", ",")
        clinical = clinical.assign(patient_medications=medication_col) # Add it in as a new "patient_medications" column
        clinical = clinical.drop(columns=["MS.Directory.Name", "Batch", "Data.Set", "TMT.Channel", "Mass.Spectrometer", "Mass.Spectrometer.Operator", "Set.A", "Set.B", "tissue_type"]) # tissue_type column is a duplicate of Sample_Tumor_Normal
        self._data["clinical"] = clinical # Save this final amalgamation of all the clinical dataframes

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
        tran_map = get_reindex_map(clinical["RNA.ID"])
        tran = self._data["transcriptomics"]
        tran_reindexed = reindex_dataframe(tran, tran_map, new_index_name="Patient_ID", keep_old=False)
        self._data["transcriptomics"] = tran_reindexed

        # Use the Specimen.Label columns from clinical dataframe to reindex the proteomics, phosphoproteomics, and phosphoproteomics_gene dataframes with patient ids
        specimen_label_map = get_reindex_map(clinical["Specimen.Label"])

        prot = self._data["proteomics"]
        phos = self._data["phosphoproteomics"]
        phos_gene = self._data["phosphoproteomics_gene"]

        prot_reindexed = reindex_dataframe(prot, specimen_label_map, new_index_name="Patient_ID", keep_old=False)
        phos_reindexed = reindex_dataframe(phos, specimen_label_map, new_index_name="Patient_ID", keep_old=False)
        phos_gene_reindexed = reindex_dataframe(phos_gene, specimen_label_map, new_index_name="Patient_ID", keep_old=False)

        self._data["proteomics"] = prot_reindexed
        self._data["phosphoproteomics"] = phos_reindexed
        self._data["phosphoproteomics_gene"] = phos_gene_reindexed

        # Now that we've used the RNA.ID and Specimen.Label columns to reindex the dataframes that needed it, we can drop them from the clinical dataframe
        clinical = clinical.drop(columns=["Specimen.Label", "RNA.ID"])

        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data)

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
        clinical = clinical.reindex(master_index)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = clinical

        # Generate a sample ID for each patient ID
        sample_id_dict = generate_sample_id_map(master_index)

        # Give all the dataframes Sample_ID indices
        dfs_to_delete = [] #
        for name in self._data.keys(): # Only loop over keys, to avoid changing the structure of the object we're looping over
            df = self._data[name]
            df.index.name = "Patient_ID"
            keep_old = name == "clinical" # Keep the old Patient_ID index as a column in the clinical dataframe, so we have a record of it.

            df = reindex_dataframe(df, sample_id_dict, "Sample_ID", keep_old)
            if df is None:
                print(f"Error mapping sample ids in {name} dataframe. At least one Patient_ID did not have corresponding Sample_ID mapped in clinical dataframe. {name} dataframe not loaded.")
                dfs_to_delete.append(name)
                continue

            self._data[name] = df

        for name in dfs_to_delete: # Delete any dataframes that had issues reindexing
            del self._data[name]

        # Drop name of column axis for all dataframes
        for name in self._data.keys():
            df = self._data[name]
            df.columns.name = None
            self._data[name] = df

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
