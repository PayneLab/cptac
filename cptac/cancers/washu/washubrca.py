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

from cptac.cancer import Cancer
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError
from cptac.utils import get_boxnote_text
from cptac.cancers.mssm.mssmclinical import MssmClinical


class WashuBrca(Source):

    def __init__(self, no_internet, version):
        """Load all of the bcmbrca dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """
        #ignore logging messages
        logger = logging.getLogger()
        logger.setLevel(logging.CRITICAL)

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        self.valid_versions = ["1.0"]

        self.data_files = {
            "1.0": {
                #"not_used" : #"BR_total_miRNA_combined.tsv",  no file on box yet
                "cibersort"         : "CIBERSORT.Output_Abs_BR.txt",
                "CNV"               : "BR.gene_level.from_seg.filtered.tsv",
                "mapping"           : "gencode.v22.annotation.gtf.gz",
                "readme"            : ["README_miRNA", "README_CIBERSORT", "README_xCell", "README_somatic_mutation_WXS", "README_gene_expression", "README.boxnote", "README_ESTIMATE_WashU"],
                "somatic_mutation"  : "BR_prospective.dnp.annotated.exonic.addrecovercases.maf.gz",
                "transcriptomics"   : "BR_tumor_RNA-Seq_Expr_WashU_FPKM.tsv.gz",
                "tumor_purity"      : "CPTAC_pancan_RNA_tumor_purity_ESTIMATE_WashU.tsv.gz",
                "xcell"             : "BR_xCell.txt",
         
            }
        }

        self.load_functions = {
            'transcriptomics'   : self.load_transcriptomics,
            'somatic_mutation'  : self.load_somatic_mutation,
            'xcell'             : self.load_xcell,
            'cibersort'         : self.load_cibersort,
            'mapping'           : self.load_mapping,
            'CNV'               : self.load_CNV,
            'tumor_purity'      : self.load_tumor_purity,
            'readme'            : self.load_readme,
        }

        # Add version check?

        # Call the parent class __init__ function
        super().__init__(cancer_type="washubrca", source="washu", version=version, valid_versions=self.valid_versions, data_files=self.data_files, no_internet=no_internet)

        # TODO: figure out how to do refactor this part

        # get clinical df (used to slice out cancer specific patient_IDs in tumor_purity file)
        mssmclin = MssmClinical(no_internet=no_internet, version=version, filter_type='pancanbrca') #_get_version - pancandataset
        self._clinical_df = mssmclin.get_clinical()

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

            if file_name == "BR_tumor_RNA-Seq_Expr_WashU_FPKM.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.rename(columns={"gene_name": "Name","gene_id": "Database_ID"})
                df = df.set_index(["Name", "Database_ID"])
                df = df.sort_index()
                df = df.T
                df.index.name = "Patient_ID"
                #remove label for tumor samples. All samples are tumors 
                df.index = df.index.str.replace(r"-T", "", regex=True) 
                self._data["transcriptomics"] = df

            elif file_name == "BR_prospective.dnp.annotated.exonic.addrecovercases.maf.gz": # Note that we use the "file_name" variable to identify files. That way we don't have to use the whole path.
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
        
            elif file_name == "BR_xCell.txt":
                df = pd.read_csv(file_path, sep = '\t', index_col = 0) 
                df = df.transpose()
                df.columns.name = 'Name'
                df.index.name = 'Patient_ID'
                df.index = df.index.str.replace(r'-T$', '', regex=True) # remove label for tumor samples
                df.index = df.index.str.replace(r'-A$', '.N', regex=True) # change label for normal samples
                self._data["xcell"] = df
                
            elif file_name == "CIBERSORT.Output_Abs_BR.txt":
                df = pd.read_csv(file_path, sep = '\t', index_col = 0) 
                df.index.name = 'Patient_ID'
                df.columns.name = 'Name'
                df.index = df.index.str.replace(r'-T$', '', regex=True) 
                df.index = df.index.str.replace(r'-A$', '.N', regex=True)
                self._data["cibersort"] = df
                
            elif file_name == "BR.gene_level.from_seg.filtered.tsv":
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
                
            elif file_name == "CPTAC_pancan_RNA_tumor_purity_ESTIMATE_WashU.tsv.gz":
                df = pd.read_csv(file_path, sep = "\t", na_values = 'NA')
                df.Sample_ID = df.Sample_ID.str.replace(r'-T', '', regex=True) # only tumor samples in file
                df = df.set_index('Sample_ID') 
                df.index.name = 'Patient_ID' 
                # Use list of patient_ids to slice out cancers                
                patient_ids = clinical_df.index.to_list()
                df = df.loc[df.index.isin(patient_ids)]                
                self._data["tumor_purity"] = df
                
            elif file_name == "README_miRNA":
                with open(file_path, 'r') as reader:
                    self._readme_files["readme_miRNA"] = reader.read()
                    
            elif file_name == "README_CIBERSORT":
                with open(file_path, 'r') as reader:
                    self._readme_files["readme_cibersort"] = reader.read()
                    
            elif file_name == "README_xCell":
                with open(file_path, 'r') as reader:
                    self._readme_files["readme_xcell"] = reader.read()
            
            elif file_name == "README_somatic_mutation_WXS":
                with open(file_path, 'r') as reader:
                    self._readme_files["readme_somatic_mutation"] = reader.read()
                    
            elif file_name == "README_gene_expression":
                with open(file_path, 'r') as reader:
                    self._readme_files["readme_transcriptomics"] = reader.read()
               
            elif file_name == "README.boxnote":
                self._readme_files["readme_cnv"] = get_boxnote_text(file_path)
                
            elif file_name == "README_ESTIMATE_WashU":
                with open(file_path, 'r') as reader:
                    self._readme_files["readme_tumor_purity"] = reader.read()


        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = f"Formatting {self.get_cancer_type()} dataframes..."
        print(formatting_msg, end='\r')

        # CNV
        cnv = self._data["CNV"]
        gene_ids = self._helper_tables["CNV_gene_ids"]
        df = cnv.join(gene_ids,how = "left") #merge in gene_ids 
        df = df.reset_index()
        df = df.set_index(["Name", "Database_ID"]) #create multi-index
        df = df.T
        df.index.name = 'Patient_ID'
        self._data["CNV"] = df
        
        self._data = sort_all_rows_pancan(self._data)  # Sort IDs (tumor first then normal)
        
        print(" " * len(formatting_msg), end='\r') # Erase the formatting message

### NEW STUFF

    def load_transcriptomics(self):
        df_type = 'transcriptomics'
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.perform_initial_checks(df_type)

            # process the file and add it to self._data
            df = pd.read_csv(file_path, sep="\t")
            df = df.rename(columns={"gene_name": "Name","gene_id": "Database_ID"})
            df = df.set_index(["Name", "Database_ID"])
            df = df.sort_index()
            df = df.T
            df.index.name = "Patient_ID"
            #remove label for tumor samples. All samples are tumors 
            df.index = df.index.str.replace(r"-T", "", regex=True) 
            self._data["transcriptomics"] = df

    def load_somatic_mutation(self):
        df_type = 'load_somatic_mutation'
        if df_type not in self._data:
            file_path = self.perform_initial_checks(df_type)

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

    def load_xcell(self):
        df_type = 'xcell'
        if df_type not in self._data:
            file_path = self.perform_initial_checks(df_type)

            df = pd.read_csv(file_path, sep = '\t', index_col = 0) 
            df = df.transpose()
            df.columns.name = 'Name'
            df.index.name = 'Patient_ID'
            df.index = df.index.str.replace(r'-T$', '', regex=True) # remove label for tumor samples
            df.index = df.index.str.replace(r'-A$', '.N', regex=True) # change label for normal samples
            self._data["xcell"] = df

    def load_cibersort(self):
        df_type = 'cibersort'
        if df_type not in self._data:
            file_path = self.perform_initial_checks(df_type)
            
            df = pd.read_csv(file_path, sep = '\t', index_col = 0) 
            df.index.name = 'Patient_ID'
            df.columns.name = 'Name'
            df.index = df.index.str.replace(r'-T$', '', regex=True) 
            df.index = df.index.str.replace(r'-A$', '.N', regex=True)
            self._data["cibersort"] = df

    def load_mapping(self):
        df_type = 'mapping'
        if df_type not in self._data:
            file_path = self.perform_initial_checks(df_type)

            df = read_gtf(file_path)
            df = df[["gene_name","gene_id"]]
            df = df.drop_duplicates()
            df = df.rename(columns={"gene_name": "Name","gene_id": "Database_ID"})
            df = df.set_index("Name")
            self._helper_tables["CNV_gene_ids"] = df

    def load_CNV(self):
        df_type = 'CNV'
        if df_type not in self._data:
            file_path = self.perform_initial_checks(df_type)

            df = pd.read_csv(file_path, sep="\t")
            df = df.rename(columns={"Gene": "Name"})
            df = df.set_index("Name")
            cnv = df

            self.load_mapping()
            gene_ids = self._helper_tables["CNV_gene_ids"]
            df = cnv.join(gene_ids,how = "left") #merge in gene_ids 
            df = df.reset_index()
            df = df.set_index(["Name", "Database_ID"]) #create multi-index
            df = df.T
            df.index.name = 'Patient_ID'
            self._data["CNV"] = df

    def load_tumor_purity(self):
        df_type = 'tumor_purity'
        if df_type not in self._data:
            file_path = self.perform_initial_checks(df_type)

            df = pd.read_csv(file_path, sep = "\t", na_values = 'NA')
            df.Sample_ID = df.Sample_ID.str.replace(r'-T', '', regex=True) # only tumor samples in file
            df = df.set_index('Sample_ID') 
            df.index.name = 'Patient_ID' 
            # Use list of patient_ids to slice out cancers                
            patient_ids = self._clinical_df.index.to_list()
            df = df.loc[df.index.isin(patient_ids)]                
            self._data["tumor_purity"] = df
    
    def load_readme(self):
        df_type = 'readme'
        if df_type not in self._data:
            # TODO perform initial checks will 
            file_path = self.perform_initial_checks(df_type)
            # need file path to be list in order to
            elif file_name == "README_miRNA":
            with open(file_path, 'r') as reader:
                self._readme_files["readme_miRNA"] = reader.read()
                    
            elif file_name == "README_CIBERSORT":
                with open(file_path, 'r') as reader:
                    self._readme_files["readme_cibersort"] = reader.read()
                    
            elif file_name == "README_xCell":
                with open(file_path, 'r') as reader:
                    self._readme_files["readme_xcell"] = reader.read()
            
            elif file_name == "README_somatic_mutation_WXS":
                with open(file_path, 'r') as reader:
                    self._readme_files["readme_somatic_mutation"] = reader.read()
                    
            elif file_name == "README_gene_expression":
                with open(file_path, 'r') as reader:
                    self._readme_files["readme_transcriptomics"] = reader.read()
               
            elif file_name == "README.boxnote":
                self._readme_files["readme_cnv"] = get_boxnote_text(file_path)
                
            elif file_name == "README_ESTIMATE_WashU":
                with open(file_path, 'r') as reader:
                    self._readme_files["readme_tumor_purity"] = reader.read()


    



    
