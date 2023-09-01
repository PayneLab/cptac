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
from pyranges import read_gtf

from cptac.cancers.source import Source
import cptac.tools.dataframe_tools as df_tools
from cptac.cancers.mssm.mssm import Mssm

class WashuHnscc(Source):
    def __init__(self, no_internet=False):
        """Define which dataframes as are available in the self.load_functions dictionary variable, with names as keys.

        Parameters:
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """
        
        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        self.data_files = {
            "cibersort"         : "CIBERSORT.Output_Abs_HNSCC.txt.gz",
            "CNV"               : "HNSCC.gene_level.from_seg.filtered.tsv.gz",
            "mapping"           : "gencode.v22.annotation.gtf.gz",
            "mature_miRNA"      : "HNSCC_mature_miRNA_combined.tsv.gz",
            "precursor_miRNA"   : "HNSCC_precursor_miRNA_combined.tsv.gz",
            "total_miRNA"       : "HNSCC_total_miRNA_combined.tsv.gz",
            "miRNA"             : ["HNSCC_mature_miRNA_combined.tsv.gz","HNSCC_precursor_miRNA_combined.tsv.gz","HNSCC_total_miRNA_combined.tsv.gz"],
            # "readme"            : ["README_miRNA","README_CIBERSORT","README_xCell","README_somatic_mutation_WXS","README_gene_expression","README.boxnote","README_ESTIMATE_WashU"],
            "somatic_mutation"  : "HNSCC_discovery.dnp.annotated.exonic.maf.gz",
            "transcriptomics"   : ["HNSCC_NAT_RNA-Seq_Expr_WashU_FPKM.tsv.gz","HNSCC_tumor_RNA-Seq_Expr_WashU_FPKM.tsv.gz"],
            "tumor_purity"      : "CPTAC_pancan_RNA_tumor_purity_ESTIMATE_WashU.tsv.gz",
            "xcell"             : "HNSCC_xCell.txt.gz",
            "hla_typing": "hla.sample.ct.10152021.sort.tsv.gz"
        }

        #self._readme_files = {}

        self.load_functions = {
            'transcriptomics'   : self.load_transcriptomics,
            'somatic_mutation'  : self.load_somatic_mutation,
            'miRNA'             : self.load_miRNA,
            'xcell'             : self.load_xcell,
            'cibersort'         : self.load_cibersort,
            'CNV'               : self.load_CNV,
            'tumor_purity'      : self.load_tumor_purity,
            #'readme'            : self.load_readme,
            "hla_typing": self.load_hla_typing
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="hnscc", source='washu', data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

    def load_transcriptomics(self):
        """
        Load transcriptomics data.

        This function scans all file paths, loads the required transcriptomics files, performs the necessary transformations, and saves them in self._data. The files loaded include 'HNSCC_tumor_RNA-Seq_Expr_WashU_FPKM.tsv.gz' and 'HNSCC_NAT_RNA-Seq_Expr_WashU_FPKM.tsv.gz'. 

        Raises:
            ValueError: If the data loaded is not in the pandas DataFrame format.
        """

        df_type = 'transcriptomics'

        # If the df_type does not exist in self._data
        if df_type not in self._data:
            file_path_list = self.locate_files(df_type) # Locate the files

            # Loop over the list of file paths
            for file_path in file_path_list:
                file_name = os.path.basename(file_path) # Extract file name

                # Load tumor data
                if file_name == "HNSCC_tumor_RNA-Seq_Expr_WashU_FPKM.tsv.gz":
                    df = pd.read_csv(file_path, sep='\t')

                    # Change column names to match package-wide naming convention
                    df = df.rename(columns={"gene_name": "Name","gene_id": "Database_ID"})

                    # Set multi-index and sort
                    df = df.set_index(["Name", "Database_ID"]).sort_index()

                    # Transpose dataframe for easier comparison and remove label for tumor samples
                    df = df.T
                    df.index.name = "Patient_ID"
                    df.index = df.index.str.replace(r"-T", "", regex=True)

                    # Save the transformed data frame to the helper tables
                    self._helper_tables["transcriptomics_tumor"] = df

                # Load normal tissue data
                if file_name == "HNSCC_NAT_RNA-Seq_Expr_WashU_FPKM.tsv.gz":
                    df_norm = pd.read_csv(file_path, sep='\t')

                    # Change column names to match package-wide naming convention
                    df_norm = df_norm.rename(columns={"gene_name": "Name","gene_id": "Database_ID"})

                    # Set multi-index and sort
                    df_norm = df_norm.set_index(["Name", "Database_ID"]).sort_index()

                    # Transpose dataframe for easier comparison and adjust naming for normal samples
                    df_norm = df_norm.T
                    df_norm.index.name = "Patient_ID"
                    df_norm.index = df_norm.index.str.replace(r"-A", ".N", regex=True)

                    # Save the transformed data frame to the helper tables
                    self._helper_tables["transcriptomics_normal"] = df_norm

            # Combine the two transcriptomics dataframes
            rna_tumor = self._helper_tables.get("transcriptomics_tumor")
            rna_normal = self._helper_tables.get("transcriptomics_normal") 

            # Check for None values or invalid types
            if rna_tumor is None or rna_normal is None:
                raise ValueError("rna_tumor or rna_normal is None")
            if not isinstance(rna_tumor, pd.DataFrame) or not isinstance(rna_normal, pd.DataFrame):
                raise ValueError("rna_tumor or rna_normal is not a DataFrame")

            # Combine the tumor and normal data
            rna_combined = pd.concat([rna_tumor, rna_normal])

            # Save the combined data frame in self._data
            self.save_df(df_type, rna_combined)
    
    def load_somatic_mutation(self):
        df_type = 'somatic_mutation'
        if df_type not in self._data:
            file_path = self.locate_files(df_type)

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
            # save df in self._data
            self.save_df(df_type, df)

    def load_miRNA(self):
        self.load_precursor_miRNA()
        self.load_mature_miRNA()
        self.load_total_mRNA()

    def load_precursor_miRNA(self):
        df_type = 'precursor_miRNA'
        if df_type not in self._data:
            file_path = self.locate_files(df_type)

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
            all_df = pd.concat([tumor, normal])
            # save df in self._data
            self.save_df('miRNA', all_df)

    def load_mature_miRNA(self):
        df_type = 'mature_miRNA'
        if df_type not in self._data:
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, delimiter = '\t', index_col = ['Name', 'ID','Alias', 'Derives_from'])
            df = df.transpose()
            df.index = df.index.str.replace('\.T$','', regex = True)
            df.index = df.index.str.replace('\.A$','.N', regex = True)
            df.index.name = 'Patient_ID'                
            # Sort
            normal = df.loc[df.index.str.contains('\.N$', regex =True)]
            normal = normal.sort_values(by=["Patient_ID"])
            tumor = df.loc[~ df.index.str.contains('\.N$', regex =True)]
            tumor = tumor.sort_values(by=["Patient_ID"])
            all_df = pd.concat([tumor, normal])
            # save df in self._data
            self.save_df('miRNA', all_df)

    def load_total_mRNA(self):
        df_type = 'total_miRNA'
        if df_type not in self._data:
            file_path = self.locate_files(df_type)
            
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
            all_df = pd.concat([tumor, normal])
            # save df in self._data
            self.save_df('miRNA', all_df)

    def load_xcell(self):
        df_type = 'xcell'
        if df_type not in self._data:
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep = '\t', index_col = 0) 
            df = df.transpose()
            df.columns.name = 'Name'
            df.index.name = 'Patient_ID'
            df.index = df.index.str.replace(r'-T$', '', regex=True) # remove label for tumor samples
            df.index = df.index.str.replace(r'-A$', '.N', regex=True) # change label for normal samples
            # save df in self._data
            self.save_df(df_type, df)

    def load_cibersort(self):
        df_type = 'cibersort'
        if df_type not in self._data:
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep = '\t', index_col = 0) 
            df.index.name = 'Patient_ID'
            df.columns.name = 'Name'
            df.index = df.index.str.replace(r'-T$', '', regex=True) 
            df.index = df.index.str.replace(r'-A$', '.N', regex=True)
            # save df in self._data
            self.save_df(df_type, df)

    def load_mapping(self):
        df_type = 'mapping'
        if "CNV_gene_ids" not in self._helper_tables:
            file_path = self.locate_files(df_type)

            df = read_gtf(file_path)
            df = df.as_df()
            df = df[["gene_name","gene_id"]]
            df = df.drop_duplicates()
            df = df.rename(columns={"gene_name": "Name","gene_id": "Database_ID"})
            df = df.set_index("Name")
            self._helper_tables["CNV_gene_ids"] = df

    def load_CNV(self):
        df_type = 'CNV'
        if df_type not in self._data:
            file_path = self.locate_files(df_type)

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
            # save df in self._data
            self.save_df(df_type, df)

    def load_tumor_purity(self):
        df_type = 'tumor_purity'
        if df_type not in self._data:
            file_path = self.locate_files(df_type)
        
            df = pd.read_csv(file_path, sep = "\t", na_values = 'NA')
            df.Sample_ID = df.Sample_ID.str.replace(r'-T', '', regex=True) # only tumor samples in file
            df = df.set_index('Sample_ID') 
            df.index.name = 'Patient_ID'

            # get clinical df (used to slice out cancer specific patient_IDs in tumor_purity file)
            mssmclin = Mssm(filter_type='hnscc', no_internet=self.no_internet)
            clinical_df = mssmclin.get_df('clinical')                
            patient_ids = clinical_df.index.to_list()
            df = df.loc[df.index.isin(patient_ids)]

            # save df in self._data
            self.save_df(df_type, df)

    def load_hla_typing(self):
        df_type = 'hla_typing'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            # which cancer_type goes with which cancer in the mssm table
            tumor_codes = {'brca':'BR', 'ccrcc':'CCRCC',
                           'ucec':'UCEC', 'gbm':'GBM', 'hnscc':'HNSCC',
                           'lscc':'LSCC', 'luad':'LUAD', 'pdac':'PDA',
                           'hcc':'HCC', 'coad':'CO', 'ov':'OV'}

            df = pd.read_csv(file_path, sep='\t')
            df = df.loc[df['Cancer'] == tumor_codes[self.cancer_type]]
            df = df.set_index("Sample")
            df.index.name = 'Patient_ID'
            df = df.sort_values(by=["Patient_ID"])

            self.save_df(df_type, df)

        return self._data[df_type]