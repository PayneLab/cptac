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
import logging
from gtfparse import read_gtf

from cptac.dataset import Dataset
from cptac.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError
from .mssmclinical import MssmClinical


class WashuLscc(Dataset):

    def __init__(self, no_internet, version):
        """Load all of the washulscc dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """
        #ignore logging messages
        logger = logging.getLogger()
        logger.setLevel(logging.CRITICAL)
        
        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        valid_versions = ["1.0"]

        data_files = {
            "1.0": [
                "LSCC_discovery.dnp.annotated.exonic.maf.gz",
                "LSCC_tumor_RNA-Seq_Expr_WashU_FPKM.tsv.gz",
                "LSCC_NAT_RNA-Seq_Expr_WashU_FPKM.tsv.gz",
                "LSCC_mature_miRNA_combined.tsv",
                "LSCC_precursor_miRNA_combined.tsv",
                "LSCC_total_miRNA_combined.tsv",
                "LSCC_xCell.txt",
                "CIBERSORT.Output_Abs_LSCC.txt",
                "LSCC.gene_level.from_seg.filtered.tsv",
                "gencode.v22.annotation.gtf.gz",
                "CPTAC_pancan_RNA_tumor_purity_ESTIMATE_WashU.tsv.gz"
            ]
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="washulscc", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)
        
        # get clinical df (used to slice out cancer specific patient_IDs in tumor_purity file)
        mssmclin = MssmClinical(no_internet=no_internet, version="latest", filter_type='pancanlscc') #_get_version - pancandataset
        clinical_df = mssmclin.get_clinical()
        

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

            

            if file_name == "LSCC_discovery.dnp.annotated.exonic.maf.gz": # Note that we use the "file_name" variable to identify files. That way we don't have to use the whole path.
                df = pd.read_csv(file_path, sep='\t')    
                # Rename the columns we want to keep to the appropriate names
                df = pd.read_csv(file_path, sep='\t')    
                df['Patient_ID'] = df.loc[:, 'Tumor_Sample_Barcode']
                df = df.rename(columns={
                         "Hugo_Symbol":"Gene",
                         "Gene":"Gene_Database_ID",
                         "Variant_Classification":"Mutation",
                         "HGVSp_Short":"Location"})

                df = df.set_index("Patient_ID")
                df = df[ ['Gene'] + ["Mutation"] + ["Location"] + [ col for col in df.columns if col not in ["Gene","Mutation","Location"] ] ]
                df.index = df.index.str.replace(r"_T", "", regex=True)     
                self._data["somatic_mutation"] = df
             
           
            if file_name == "LSCC_tumor_RNA-Seq_Expr_WashU_FPKM.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.rename(columns={"gene_name": "Name","gene_id": "Database_ID"})
                df = df.set_index(["Name", "Database_ID"])
                df = df.sort_index()
                df = df.T
                df.index.name = "Patient_ID"
                df.index = df.index.str.replace(r"-T", "", regex=True) #remove label for tumor samples
                self._helper_tables["transcriptomics_tumor"] = df
                
            if file_name == "LSCC_NAT_RNA-Seq_Expr_WashU_FPKM.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.rename(columns={"gene_name": "Name","gene_id": "Database_ID"})
                df = df.set_index(["Name", "Database_ID"])
                df = df.sort_index()
                df = df.T
                df.index.name = "Patient_ID"
                df.index = df.index.str.replace(r"-A", ".N", regex=True) #remove label for tumor samples
                self._helper_tables["transcriptomics_norm"] = df
                
                
            elif 'miRNA_combined' in file_name:
                miRNA_type = file_name.split('_')[1] # get type of miRNA data (precursor, mature, or total)
                if miRNA_type == 'mature':
                    df = pd.read_csv(file_path, delimiter = '\t', index_col = ['Name', 'ID','Alias', 'Derives_from'])
                else:
                    df = pd.read_csv(file_path, delimiter = '\t', index_col = ['Name', 'ID','Alias'])
                df = df.transpose()
                df.index = df.index.str.replace('\.T$','', regex = True)
                df.index = df.index.str.replace('\.A$','.N', regex = True)
                df.index.name = 'Patient_ID'                
                # Sort
                normal = df.loc[df.index.str.contains('\.N$', regex =True)]
                normal = normal.sort_values(by=["Patient_ID"])
                tumor = df.loc[~ df.index.str.contains('\.N$', regex =True)]
                tumor = tumor.sort_values(by=["Patient_ID"])
                all_df = tumor.append(normal)
                self._data[miRNA_type+'_miRNA'] = all_df
                
            elif file_name == "LSCC_xCell.txt":
                df = pd.read_csv(file_path, sep = '\t', index_col = 0) 
                df = df.transpose()
                df.columns.name = 'Name'
                df.index.name = 'Patient_ID'
                df.index = df.index.str.replace(r'-T$', '', regex=True) # remove label for tumor samples
                df.index = df.index.str.replace(r'-A$', '.N', regex=True) # change label for normal samples
                self._data["xcell"] = df
                
            elif file_name == "CIBERSORT.Output_Abs_LSCC.txt":
                df = pd.read_csv(file_path, sep = '\t', index_col = 0) 
                df.index.name = 'Patient_ID'
                df.columns.name = 'Name'
                df.index = df.index.str.replace(r'-T$', '', regex=True) 
                df.index = df.index.str.replace(r'-A$', '.N', regex=True)
                self._data["cibersort"] = df
            #CNV    
            elif file_name == "LSCC.gene_level.from_seg.filtered.tsv":
                df = pd.read_csv(file_path, sep="\t")
                df = df.rename(columns={"Gene": "Name"})
                df = df.set_index("Name")
                self._data["CNV"] = df
                
            elif file_name == "gencode.v22.annotation.gtf.gz":
                df = read_gtf(file_path)
                df = df[["gene_name","gene_id"]]
                df = df.drop_duplicates()
                df = df.rename(columns={"gene_name": "Name","gene_id": "Database_ID"})
                df = df.set_index("Name")
                self._helper_tables["CNV_gene_ids"] = df
                
            # tumor_purity    
            elif file_name == "CPTAC_pancan_RNA_tumor_purity_ESTIMATE_WashU.tsv.gz":
                df = pd.read_csv(file_path, sep = "\t", na_values = 'NA')
                df.Sample_ID = df.Sample_ID.str.replace(r'-T', '', regex=True) # only tumor samples in file
                df = df.set_index('Sample_ID') 
                df.index.name = 'Patient_ID' 
                # Use list of patient_ids to slice out cancers                
                patient_ids = clinical_df.index.to_list()
                df = df.loc[df.index.isin(patient_ids)]                
                self._data["tumor_purity"] = df
        
# combine and create transcriptomic dataframe            
        rna_tumor = self._helper_tables.get("transcriptomics_tumor")
        rna_normal = self._helper_tables.get("transcriptomics_norm") # Normal entries are already marked with 'N' on the end of the ID
        rna_combined = rna_tumor.append(rna_normal)
        self._data["transcriptomics"] = rna_combined
      
        
        # CNV
        cnv = self._data["CNV"]
        gene_ids = self._helper_tables["CNV_gene_ids"]
        df = cnv.join(gene_ids,how = "left") #merge in gene_ids 
        df = df.reset_index()
        df = df.set_index(["Name", "Database_ID"]) #create multi-index
        df = df.T
        df.index.name = 'Patient_ID'
        self._data["CNV"] = df
        
        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        
        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
#        new_clinical = self._data["clinical"]
#        new_clinical = new_clinical.reindex(master_index)

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
#        sample_status_col = generate_sample_status_col(new_clinical, normal_test=lambda sample: sample[0] == 'N')
        ###END EXAMPLE CODE#####################################################

#        new_clinical.insert(0, "Sample_Tumor_Normal", sample_status_col)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
#        self._data['clinical'] = new_clinical

        # Edit the format of the Patient_IDs to have normal samples marked the same way as in other datasets. 
        
        ###FILL: You will need to pass the proper parameters to correctly
        ### reformat the patient IDs in your dataset. The standard format is to
        ### have the string '.N' appended to the end of the normal patient IDs,
        ### e.g. the  normal patient ID corresponding to C3L-00378 would be
        ### C3L-00378.N (this way we can easily match two samples from the same
        ### patient). The example code below is for a dataset where all the
        ### normal samples have  an "N" prepended to the patient IDs. The
        ### reformat_normal_patient_ids function erases that and puts a ".N" at
        ### the end. See cptac/dataframe_tools.py for further function
        ### documentation.
        ###START EXAMPLE CODE###################################################
#        self._data = reformat_normal_patient_ids(self._data, existing_identifier="N", existing_identifier_location="start")
        ###END EXAMPLE CODE#####################################################

        # Call function from dataframe_tools.py to sort all tables first by sample status, and then by the index
#        self._data = sort_all_rows(self._data)

        # Call function from dataframe_tools.py to standardize the names of the index and column axes
#        self._data = standardize_axes_names(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
