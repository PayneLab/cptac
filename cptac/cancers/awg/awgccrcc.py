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

class AwgCcrcc(Source):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the ccrcc dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        self.valid_versions = ["0.0", "0.1", "0.1.1"] # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.

        self.data_files = {
            "0.0": {
                "phosphoproteomics_gene" : "6_CPTAC3_CCRCC_Phospho_abundance_gene_protNorm=2_CB_imputed.tsv.gz",
                "phosphoproteomics"      : "6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv.gz",
                "proteomics"             : "6_CPTAC3_CCRCC_Whole_abundance_protein_pep=unique_protNorm=2_CB.tsv.gz",
                "annotation"             : ["Clinical Table S1.xlsx", "cptac-metadata.xls.gz", "S044_CPTAC_CCRCC_Discovery_Cohort_Clinical_Data_r3_Mar2019.xlsx"],
                "somatic_mutation"       : "ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf.gz",
                "methylation"            : "ccrccMethylGeneLevelByMean.txt.gz",
                "CNV"                    : "kirc_wgs_cnv_gene.csv.gz",
                "transcriptomics"        : "RNA_Normal_Tumor_185_samples.tsv.gz"},
            "0.1": {
                "phosphoproteomics_gene" : "6_CPTAC3_CCRCC_Phospho_abundance_gene_protNorm=2_CB_imputed.tsv.gz",
                "phosphoproteomics"      : "6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv.gz",
                "proteomics"             : "6_CPTAC3_CCRCC_Whole_abundance_protein_pep=unique_protNorm=2_CB.tsv.gz",
                "methylation"            : "ccrccMethylGeneLevelByMean.txt.gz",
                "somatic_mutation"       : "ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf.gz",
                "annotation"             : ["Clinical Table S1.xlsx", "cptac-metadata.xls.gz", "S044_CPTAC_CCRCC_Discovery_Cohort_Clinical_Data_r3_Mar2019.xlsx", "Table S7.xlsx"],
                "CNV"                    : "kirc_wgs_cnv_gene.csv.gz",
                "transcriptomics"        : "RNA_Normal_Tumor_185_samples.tsv.gz"},
            "0.1.1": {
                "phosphoproteomics_gene" : "6_CPTAC3_CCRCC_Phospho_abundance_gene_protNorm=2_CB_imputed.tsv.gz",
                "phosphoproteomics"      : "6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv.gz",
                "proteomics"             : "6_CPTAC3_CCRCC_Whole_abundance_protein_pep=unique_protNorm=2_CB.tsv.gz",
                "followup"               : "CCRCC_followup_9_12.xlsx",
                "methylation"            : "ccrccMethylGeneLevelByMean.txt.gz",
                "somatic_mutation"       : "ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf.gz",
                "annotation"             : ["Clinical Table S1.xlsx", "cptac-metadata.xls.gz", "S044_CPTAC_CCRCC_Discovery_Cohort_Clinical_Data_r3_Mar2019.xlsx", "Table S7.xlsx"],
                "CNV"                    : "kirc_wgs_cnv_gene.csv.gz",
                "transcriptomics"        : "RNA_Normal_Tumor_185_samples.tsv.gz"},
        }

        self.load_functions = {
            "clinical"               : self.load_annotation,
            "medical_history"        : self.load_annotation,
            "CNV"                    : self.load_CNV,
            "followup"               : self.load_followup,
            "methylation"            : self.load_methylation,
            "phosphoproteomics"      : self.load_phosphoproteomics,
            "phosphoproteomics_gene" : self.load_phosphoproteomics_gene,
            "proteomics"             : self.load_proteomics,
            "somatic_mutation"       : self.load_somatic_mutation,
            "transcriptomics"        : self.load_transcriptomics,
        }

        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        super().__init__(cancer_type="ccrcc", source='awg', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

        # We're going to need to drop the samples below from a couple dataframes
        # TODO: WHY?
        self.nci_labels = ["NCI7-1", "NCI7-2", "NCI7-3", "NCI7-4", "NCI7-5"]
        self.nci_dotted_labels = [label.replace("-", ".") for label in self.nci_labels]
        self.qc_labels = ["QC1", "QC2", "QC3", "QC4", "QC5", "QC6", "QC7", "QC8"]

        # specimen labels and RNA ids for reindexing other dfs
        # loaded during self.load_annotation
        self.specimen_labels = None
        self.rna_ids = None

        # TODO: Why is proteomics used to drop unique columns from CNV and methylation and not clinical?

    def load_annotation(self):
        df_type = 'annotation'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_paths = self.locate_files(df_type)

            # process the files
            clinical_dfs = {}
            for file_path in file_paths:
                file_name = file_path.split(os.sep)[-1]
            
                if file_name == "Clinical Table S1.xlsx":
                    df = pd.\
                        read_excel(file_path, sheet_name="ccrcc_clinical_characteristics", index_col=0).\
                        dropna(how="all", axis=0).\
                        dropna(how="all", axis=1)
                    df.index.name = "Patient_ID" # The index is currently "case_id", but we call that "Patient_ID"
                    clinical_dfs["authoritative_clinical"] = df

                elif file_name == "cptac-metadata.xls.gz":
                    df = pd.read_csv(file_path, index_col=0)
                    df = df.drop(index=["pooled sample"] + self.nci_labels + self.qc_labels) # Drop the pooled samples in addition to other samples to exclude
                    clinical_dfs["metadata_and_keys"] = df
                
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
                    immune_groups = pd.read_excel(file_path, sheet_name="xCell Signatures", index_col = 0, header=None, skiprows=range(0, 3)).transpose()
                    immune_groups = immune_groups[["Samples", "Immune Group"]] # We only need these columns
                    immune_groups = immune_groups.set_index("Samples")

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
            if self.version == "0.1":
                clinical = clinical.join(immune_groups, on="Specimen.Label", how="outer")

            # Save the RNA.ID and Specimen.Label columns to reindex the dataframes that need it, and drop from clinical
            self.specimen_labels = clinical["Specimen.Label"].copy(deep=True)
            self.rna_ids = clinical["RNA.ID"].copy(deep=True)
            clinical = clinical.drop(columns=["Specimen.Label", "RNA.ID"])
            
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

            # Move the prepended N to a .N at the end to match other normal sample labeling in cptac
            clinical.index = clinical.index.where(~clinical.index.str.startswith('N'), clinical.index.str[1:] + ".N")

            # save dfs in self._data
            self.save_df('clinical', clinical)
            self.save_df("medical_history", medical_history_parsed)
 

    def load_phosphoproteomics_gene(self):
        # annotation dfs required to fully load proteomics
        self.load_annotation()

        df_type = 'phosphoproteomics_gene'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # process the file
            df = pd.read_csv(file_path, sep='\t', index_col=0)
            ref_intensities = df["ReferenceIntensity"] # Copy this out, so we can subtract the reference intensities later
            df = df.drop(columns=["NumberPSM", "Proteins", "ReferenceIntensity"] + self.nci_dotted_labels + self.qc_labels)
            df = df.subtract(ref_intensities, axis="index") # Subtract reference intensities from all the values, to get ratios
            df = df.transpose()

             # reindex with clinical data
            df = self._specimen_reindex(df_type, df)

            # Move the prepended N to a .N at the end to match other normal sample labeling in cptac
            df.index = df.index.where(~df.index.str.startswith('N'), df.index.str[1:] + ".N")

            # save df in self._data
            self.save_df(df_type, df)

    def load_phosphoproteomics(self):
        # annotation dfs required to fully load proteomics
        self.load_annotation()

        df_type = 'phosphoproteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # process the file
            df = pd.read_csv(file_path, sep='\t')

            # Drop unlocalized sites
            unlocalized_sites = (df["Index"].str.rsplit("_", n=1, expand=True)[1] == '0') 
            df = df[~unlocalized_sites]

            # Drop unwanted samples
            df = df.drop(columns=self.nci_labels + self.qc_labels)

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

             # reindex with clinical data
            df = self._specimen_reindex(df_type, df)

            # save df in self._data
            self.save_df(df_type, df)

    def load_proteomics(self):
        # annotation dfs required to fully load proteomics
        self.load_annotation()

        df_type = 'proteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # process the file
            df = pd.read_csv(file_path, sep='\t')
            df = df.rename(columns={"Proteins": "Name", "Index": "Database_ID"})
            df = df.set_index(["Name", "Database_ID"])
            ref_intensities = df["ReferenceIntensity"] # Copy this out, so we can subtract the reference intensities later
            df = df.drop(columns=["NumberPSM", "ReferenceIntensity"] + self.nci_labels + self.qc_labels)
            df = df.subtract(ref_intensities, axis="index") # Subtract reference intensities from all the values, to get ratios
            df = df.sort_index()
            df = df.transpose()
            df = df.sort_index()
            
            # reindex with clinical data
            df = self._specimen_reindex(df_type, df)

            # Move the prepended N to a .N at the end to match other normal sample labeling in cptac
            df.index = df.index.where(~df.index.str.startswith('N'), df.index.str[1:] + ".N")

            # save df in self._data
            self.save_df(df_type, df)

    def load_somatic_mutation(self):
        df_type = 'somatic_mutation'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # process the file
            df = pd.read_csv(file_path, sep='\t', dtype={"PUBMED":object}) # "PUBMED" column has mixed types, so we specify object as the dtype to prevent a warning from printing. We don't actually use the column, so that's all we need to do.
            split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n=1, expand=True) # The first part of the barcode is the patient id, which we need want to make the index
            df["Tumor_Sample_Barcode"] = split_barcode[0]
            df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
            df = df.rename(columns={"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"})                
            df = df.sort_values(by=["Patient_ID", "Gene"])
            df = df.set_index("Patient_ID")
            
            # save df in self._data
            self.save_df(df_type, df)

    def load_methylation(self):
        # proteomics dfs required to fully load methylation
        self.load_proteomics()

        df_type = 'methylation'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # process the file
            df = pd.read_csv(file_path, sep='\t')
            df = df.sort_index()
            df = df.transpose()
            df.index.name = "Patient_ID"

            # drop 26 rows that aren't in any of the other tables (except CNV) and were excluded from the original analysis
            to_drop = df.index.difference(self._data['proteomics'].index)
            df = df.drop(index=to_drop)
            
            # save df in self._data
            self.save_df(df_type, df)

    def load_CNV(self):
        # proteomics df required to complete the CNV index
        self.load_proteomics()

        df_type = 'CNV'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # process the file
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

            # drop 28 rows that don't show up in the other tables (except 26 of the 28 which are in methylation) 
            # and were excluded from the original analysis
            to_drop = df.index.difference(self._data['proteomics'].index)
            df = df.drop(index=to_drop)

            # save df in self._data
            self.save_df(df_type, df)

    def load_transcriptomics(self):
        # annotation dfs required to fully load transcriptomics df
        self.load_annotation()

        df_type = 'transcriptomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # process the file
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

            # Use the RNA.ID column from clinical dataframe to reindex transcriptomics dataframe with patient ids
            try:
                tran_map = get_reindex_map(self.rna_ids)
                df = reindex_dataframe(df, tran_map, new_index_name="Patient_ID", keep_old=False)
            except ReindexMapError:
                warnings.warn("Error mapping sample ids in transcriptomics dataframe. At least one RNA.ID did not have a corresponding Patient_ID mapped in the clinical dataframe. transcriptomics dataframe not loaded.", FailedReindexWarning, stacklevel=2)

            # Move the prepended N to a .N at the end to match other normal sample labeling in cptac
            df.index = df.index.where(~df.index.str.startswith('N'), df.index.str[1:] + ".N")

            # save df in self._data
            self.save_df(df_type, df)

    def load_followup(self): 
        df_type = 'followup'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # process the file elif file_name == 'CCRCC_followup_9_12.xlsx' and self.version == "0.1.1":
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

            # save df in self._data
            self.save_df(df_type, df)

    def how_to_cite(self):
        return self.super().how_to_cite(cancer_type='clear cell renal cell carcinoma (kidney)', pmid=31675502)

    def _specimen_reindex(self, df_type, df):
        """Use the Specimen.Label columns from clinical dataframe to reindex the 
        proteomics, phosphoproteomics, and phosphoproteomics_gene dataframes with patient ids

        Returns: reindexed df
        """
        try:
            specimen_label_map = get_reindex_map(self.specimen_labels)
        except ReindexMapError:
            warnings.warn(f"Error mapping sample ids in {df_type} dataframe. Specimen.Label mapping in clinical dataframe was not one-to-one. Dataframe not loaded.", FailedReindexWarning, stacklevel=2)
        
        try:
            df = reindex_dataframe(df, specimen_label_map, new_index_name="Patient_ID", keep_old=False)
        except ReindexMapError as error:
            warnings.warn(f"Error mapping sample ids in {df_type} dataframe. RNA.ID {str(error)} did not have a corresponding Patient_ID mapped in the clinical dataframe. {df_type} dataframe not loaded.", FailedReindexWarning, stacklevel=2)
        
        return df