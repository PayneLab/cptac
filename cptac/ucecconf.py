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
import datetime

from .dataset import Dataset
from .dataframe_tools import *
from .exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError, InvalidParameterError

class UcecConf(Dataset):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        valid_versions = ["1.0"]

        data_files = {
            "1.0": [
            #"UCEC_confirmatory_acetyl_gene_ratio_median_polishing_log2_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
            #"UCEC_confirmatory_Direct_SRM_tumor_v1.0.cct.gz", #SRM not to be included in 1.0
            #"UCEC_confirmatory_IMAC_SRM_tumor_v1.0.cct.gz",
            "UCEC_confirmatory_meta_table_v1.0.xlsx",
            "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v1.0.cct.gz",
            "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v1.0.cct.gz",
            #"UCEC_confirmatory_nglycoform-site_ratio_median_polishing_log2_tumor_normal_v1.0.cct.gz",
            #"UCEC_confirmatory_phospho_gene_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_phospho_site_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_proteomics_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v1.0.txt.gz",
            "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v1.0.cct.gz",
            #"UCEC_confirmatory_RNAseq_isoform_FPKM_removed_circRNA_log2(x+1)_tumor_normal_v1.0.cct.gz",
            "UCEC_confirmatory_WES_cnv_gistic_thresholded_tumor_v1.0.cct.gz",
            "UCEC_confirmatory_WES_cnv_log2_ratio_tumor_v1.0.cct.gz",
            "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.0.cbt.gz",
            "UCEC_confirmatory_WES_somatic_mutation_v1.0.maf.gz",
            #"UCEC_confirmatory_WGS_SV_tumor_v1.0.txt.gz" #structural_variation - not to be included in 1.0
            ],
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="ucecconf", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

            if file_name == "UCEC_confirmatory_acetyl_site_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.reset_index()
                df[['Name','Database_ID','Site']] = df.idx.str.split("@", expand=True)
                df['Site'] = df['Site'].str.split('-',expand=True)[1]
                df = df.set_index(["Name", "Site", "Database_ID"])
                df = df.drop(columns=["idx"])
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                self._data["acetylproteomics"] = df
            
            elif file_name == "UCEC_confirmatory_meta_table_v1.0.xlsx":
                df = pd.read_excel(file_path)
                df.insert(6, "Proteomics_Tumor_Normal", df["Group"])
                df.loc[df['Group'] == 'Enriched_Normal', 'Idx'] = df['Idx'] + '.N'
                df.loc[df['Group'] == 'Adjacent_normal', 'Idx'] = df['Idx'].str[:-2] + '.N'
                df = df.set_index("Idx")
                df.loc[df['Group'] != 'Tumor', 'Group'] = 'Normal'
                df = df.rename({'Group': 'Sample_Tumor_Normal'}, axis=1)
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["clinical"] = df
                
            elif file_name == "UCEC_confirmatory_methylation_gene_level_beta_value_tumor_v1.0.cct.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0, dtype=object)
                #Not sure what the dtype thing actually does, but it was giving me a big warning
                #This made it go away. I need to find out exactly what is happening
                df = df.transpose()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                #I also am not 100% sure this is methylation, but I don't know what else it would be
                self._data["methylation"] = df
            
            elif file_name == "UCEC_confirmatory_miRNAseq_miRNA_TPM_log2(x+1)_tumor_normal_v1.0.cct.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                #If we wanted to turn every -A into a .N on the Patient_IDs
                #df = df.reset_index()
                #df.loc[df['Patient_ID'].str[-2:] == '-A', 'Patient_ID'] = df['Patient_ID'].str[:-2] + '.N'
                #df = df.set_index("Patient_ID")
                df.columns.name = "Name"
                self._data["miRNA"] = df
                
            elif file_name == "UCEC_confirmatory_phospho_site_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.reset_index()
                df[['Name','Database_ID','Site']] = df.idx.str.split("@", expand=True)
                df['Site'] = df['Site'].str.split('-',expand=True)[1]
                df = df.set_index(["Name", "Site", "Database_ID"])
                df = df.drop(columns=["idx"])
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                self._data["phosphoproteomics"] = df
                
            elif file_name == "UCEC_confirmatory_proteomics_ratio_median_polishing_log22_tumor_normal_v1.0.cct.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["proteomics"] = df
                
            elif file_name == "UCEC_confirmatory_RNAseq_circRNA_RSEM_UQ_log2(x+1)_tumor_normal_v1.0.cct.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["circular_RNA"] = df
                
            elif file_name == "UCEC_confirmatory_RNAseq_gene_fusion_tumor_v1.0.txt.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.reset_index()
                df = df.set_index("Sample")
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["gene_fusion"] = df

            elif file_name == "UCEC_confirmatory_RNAseq_gene_RSEM_removed_circRNA_UQ_log2(x+1)_tumor_normal_v1.0.cct.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["transcriptomics"] = df
                
            elif file_name == "UCEC_confirmatory_WES_cnv_gistic_thresholded_tumor_v1.0.cct.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["CNV_gistic"] = df
                
            elif file_name == "UCEC_confirmatory_WES_cnv_log2_ratio_tumor_v1.0.cct.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["CNV_log2ratio"] = df
                
            elif file_name == "UCEC_confirmatory_WES_somatic_mutation_gene_level_V1.0.cbt.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df.index.name = "Patient_ID"
                df.columns.name = "Name"
                self._data["somatic_mutation_binary"] = df

            elif file_name == "UCEC_confirmatory_WES_somatic_mutation_v1.0.maf.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.reset_index()
                df['Tumor_Sample_Barcode'].apply(lambda s: s[:-2])
                df = df[['Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification','HGVSp_Short']]
                df = df.rename(columns={
                    "Tumor_Sample_Barcode": "Patient_ID",
                    "Hugo_Symbol":"Gene",
                    "Variant_Classification":"Mutation",
                    "HGVSp_Short":"Location"})
                df = df.sort_values(by=["Patient_ID", "Gene"])
                df.columns.name = "Name"
                df = df.set_index("Patient_ID")
                self._data["somatic_mutation"] = df

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')


        # NOTE: The code below will not work properly until you have all the
        # dataframes formatted properly and loaded into the self._data
        # dictionary. That's why they're commented out for now. Go ahead and
        # uncomment them when all the data tables are ready. Note that some of
        # the lines are marked as just examples, though, and you'll still need
        # to adapt them to your specific situation.

        # Get a union of all dataframes' indices, with duplicates removed
        ###FILL: If there are any tables whose index values you don't want
        ### included in the master index, pass them to the optional 'exclude'
        ### parameter of the unionize_indices function. This was useful, for
        ### example, when some datasets' followup data files included samples
        ### from cohorts that weren't in any data tables besides the followup
        ### table, so we excluded the followup table from the master index since
        ### there wasn't any point in creating empty representative rows for
        ### those samples just because they existed in the followup table.

        # master_index = unionize_indices(self._data)

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.

        # new_clinical = self._data["clinical"]
        # new_clinical = new_clinical.reindex(master_index)

        # Add a column called Sample_Tumor_Normal to the clinical dataframe indicating whether each sample was a tumor or normal sample. Use a function from dataframe_tools to generate it.

        ###FILL: Your dataset should have some way that it marks the Patient IDs
        ### of normal samples. The example code below is for a dataset that
        ### marks them by putting an 'N' at the beginning of each one. You will
        ### need to write a lambda function that takes a given Patient_ID string
        ### and returns a bool indicating whether it corresponds to a normal
        ### sample. Pass that lambda function to the 'normal_test' parameter of
        ### the  generate_sample_status_col function when you call it. See
        ### cptac/dataframe_tools.py for further function documentation.
        ###START EXAMPLE CODE###################################################

        # sample_status_col = generate_sample_status_col(new_clinical, normal_test=lambda sample: sample[0] == 'N')

        ###END EXAMPLE CODE#####################################################

        # new_clinical.insert(0, "Sample_Tumor_Normal", sample_status_col)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        # self._data['clinical'] = new_clinical

        # Edit the format of the Patient_IDs to have normal samples marked the same way as in other datasets.

        ###FILL: You may need to use the code below to reformat the patient IDs
        ### in your dataset. This applies if all of the normal samples are
        ### already marked in the original data files in some way, but just not
        ### in the way we want (e.g. they have an "N" at the beginning of the
        ### sample ID, instead of a ".N" at the end). Be aware that the case
        ### with some datasets such as PDAC is different; instead of the normal
        ### samples already being marked, just not in the way we want, they're
        ### actually contained in a separate table, with no special marking on
        ### the sample ids. In those cases you wouldn't use the
        ### reformat_normal_patient_ids function, and would instead just mark
        ### the samples in the normal tables with the ".N" before appending them
        ### to the tumor tables.
        ### If you do use this function: the standard normal ID format is to
        ### have the string '.N' appended to the end of the normal patient IDs,
        ### e.g. the  normal patient ID corresponding to C3L-00378 would be
        ### C3L-00378.N (this way we can easily match two samples from the same
        ### patient). The example code below is for a dataset where all the
        ### normal samples have  an "N" prepended to the patient IDs. The
        ### reformat_normal_patient_ids function erases that and puts a ".N" at
        ### the end. See cptac/dataframe_tools.py for further function
        ### documentation.
        ###START EXAMPLE CODE###################################################
        # self._data = reformat_normal_patient_ids(self._data, existing_identifier="N", existing_identifier_location="start")
        ###END EXAMPLE CODE#####################################################

        # Call function from dataframe_tools.py to sort all tables first by sample status, and then by the index
        # self._data = sort_all_rows(self._data)

        # Call function from dataframe_tools.py to standardize the names of the index and column axes
        # self._data = standardize_axes_names(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message

        # Print password access only warning
        warnings.warn("The UcecConf data is currently strictly reserved for CPTAC investigators. Otherwise, you are not authorized to access these data. Additionally, even after these data become publicly available, they will be subject to a publication embargo (see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter cptac.embargo() to open the webpage for more details).", PublicationEmbargoWarning, stacklevel=2)
        
    def get_CNV(self, algorithm=None):
        if (not algorithm):
            message = ("Please specify which type of UcecConf CNV data you want: "
            "'log2ratio' or 'gistic'. i.e. get_CNV('gistic')")
            raise InvalidParameterError(message)
        elif (algorithm == "log2ratio"):
            return super().get_dataframe("CNV_log2ratio")
        elif (algorithm == "gistic"):
            return super().get_dataframe("CNV_gistic")
        else: 
            message = ("Please specify a valid algorithm type for UcecConf CNV data: "
            "'log2ratio' or 'gistic'. i.e. get_CNV('gistic')")
            raise InvalidParameterError(message)

    def how_to_cite(self):
        super().how_to_cite(cancer_type='endometrial confirmatory carcinoma', pmid='', unpublished=True)
        

