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

from cptac.cancers.source import Source
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, ReindexMapError, PublicationEmbargoWarning

class AwgHnscc(Source):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the hnscc dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        self.valid_versions = ["0.1", "2.0"] # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.

        self.data_files = {
            "0.1": {
                "somatic_mutation"  : "HNSCC.strelka.sorted.filtered.annovar.hg19_multianno_filtered.maf.txt.gz",
                "proteomics"        : ["Proteomics_DIA_Gene_level_Normal.cct.gz", "Proteomics_DIA_Gene_level_Tumor.cct.gz"],
                "transcriptomics"   : "RNAseq_RSEM_UQ_log2.cct.gz",
                "circular_RNA"      : "RNAseq_circ_RSEM_UQ_log2.cct.gz",
                "CNV"               : "SCNA_gene_level.cct.gz",
                "annotation"        : "clinic.tsi.gz"},
            "2.0": {
                "circular_RNA"      : "circRNAseq_RSEM_UQ_log2_Combined.cct.gz",
                "followup"          : "HN_followUp_9_24.xlsx",
                "annotation"        : "Meta_table.tsv.gz",
                "miRNA"             : "microRNA_log2_Combined.cct.gz",
                "phosphoproteomics" : "Phosphoproteomics_TMT_site_level_combined_all.cct.gz",
                "proteomics"        : "Proteomics_TMT_gene_level_combined_all.cct.gz",
                "transcriptomics"   : "RNAseq_RSEM_UQ_Combined.cct.gz",
                "CNV"               : "SCNA_log2_gene_level.cct.gz",
                "somatic_mutation"  : "SomaticMutations_maf.tsv.gz"},
        }

        self.load_functions = {
            'clinical'                : self.load_clinical,
            'derived_molecular'       : self.load_derived_molecular,
            'CNV'                     : self.load_CNV,
            'followup'                : self.load_followup,
            'circular_RNA'            : self.load_circular_RNA,
            'miRNA'                   : self.load_miRNA,
            'phosphoproteomics'       : self.load_phosphoproteomics,
            'proteomics'              : self.load_proteomics,
            'somatic_mutation'        : self.load_somatic_mutation,
            'transcriptomics'         : self.load_transcriptomics,
        }

        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        super().__init__(cancer_type="hnscc", source='awg', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)


    def load_CNV(self):
        df_type = 'CNV'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep='\t')

            if self.version == "2.0":
                df = df.set_index('gene_symbol')

            df = df.sort_index()
            df = df.transpose()
            df = df.sort_index()
            df.columns.name=None
            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_miRNA(self):
        df_type = 'miRNA'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.sort_index()
            df = df.transpose()

            # Reformat patient ids
            df.index = df.index.str.replace(r'-T$', '', 1, regex=True)
            df.index = df.index.str.replace(r'-N$', '.N', 1, regex=True)

            # save df in self._data
            self.save_df(df_type, df)


    def load_transcriptomics(self):
        df_type = 'transcriptomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep='\t')

            if self.version == "2.0":
                df = df.set_index('Idx')

            df = df.sort_index()
            df = df.transpose()
            df = df.sort_index()
            df.columns.name=None

            if self.version == "0.1":
                df.index = df.index.str.replace(r'\.', '-', 1, regex=True)
                df.index = df.index.str.replace(r'\.T$', '', 1, regex=True)
            elif self.version == "2.0":
                df.index = df.index.str.replace(r'-T$', '', 1, regex=True)
                df.index = df.index.str.replace(r'-N$', '.N', 1, regex=True)

            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_circular_RNA(self):
        df_type = 'circular_RNA'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df = df.sort_index()
            df = df.transpose()
            df = df.sort_index()
            df.columns.name=None

            if self.version == "0.1":
                # We want all the patientIDs to have the the format C3L-00977, and these have the form C3L.00977.N, so we need to replace the first "." with a "-"
                df.index = df.index.str.replace(r'\.', '-', 1, regex=True)
                df.index = df.index.str.replace(r'\.T$', '', 1, regex=True)

            elif self.version == "2.0":
                df.index = df.index.str.replace(r'-T$', '', 1, regex=True)
                df.index = df.index.str.replace(r'-N$', '.N', 1, regex=True)

            df.index.name = "Patient_ID"

            # save df in self._data
            self.save_df(df_type, df)


    def load_followup(self):
        df_type = 'followup'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_excel(file_path)

            # Rename, set, and sort by index
            df = df.rename(columns={"CASE_ID": "Patient_ID"})
            df = df.set_index("Patient_ID")
            df = df.sort_index()

            # save df in self._data
            self.save_df(df_type, df)


    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep='\t')
            df = df.rename(columns={"Gene": "Name"})

            # Drop unlocalized sites
            unlocalized_sites = (df["Index"].str.rsplit("_", n=1, expand=True)[1] == '0')
            df = df[~unlocalized_sites]

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
            df.index = df.index.str.replace(r'-T$', '', 1, regex=True)
            df.index = df.index.str.replace(r'-N$', '.N', 1, regex=True)
            df.index = df.index.str.replace(r'-C$', '.C', 1, regex=True) #-C is cored NAT samples
            df = df.sort_index()

            # save df in self._data
            self.save_df(df_type, df)


    def load_proteomics(self):
        df_type = 'proteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            if self.version == "0.1":
                for file in file_path:
                    path_elements = file.split(os.sep) # Get a list of the levels of the path
                    file_name = path_elements[-1]# The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below

                    # Combine the two proteomics dataframes
                    if file_name == "Proteomics_DIA_Gene_level_Normal.cct.gz":
                        df_normal = pd.read_csv(file, sep='\t')
                        df_normal = df_normal.transpose()
                        df_normal.columns.name=None
                        df_normal.index.name = "Patient_ID"

                    elif file_name == "Proteomics_DIA_Gene_level_Tumor.cct.gz":
                        df_tumor = pd.read_csv(file, sep='\t')
                        df = df.transpose()
                        df.columns.name=None
                        df.index.name = "Patient_ID"

                df_normal.index = df_normal.index + ".N" #concatenate a ".N" onto the end of the normal data so we can identify it as normal after it's appended to tumor
                prot_combined = pd.concat([df_tumor, df_normal]) #append the normal data onto the end of the tumor data
                prot_combined = prot_combined.sort_index(axis='columns') # Put all the columns in alphabetical order
                prot_combined = prot_combined.sort_index()
                df = prot_combined

            if self.version == "2.0":
                df = pd.read_csv(file_path, sep='\t')
                df = df.set_index('Index')

                df = df.transpose()
                df.columns.name=None
                df.index.name = "Patient_ID"

                df.index = df.index.str.replace(r'-T$', '', 1, regex=True)
                df.index = df.index.str.replace(r'-N$', '.N', 1, regex=True)
                df.index = df.index.str.replace(r'-C$', '.C', 1, regex=True) #-C is cored NAT samples

            # save df in self._data
            self.save_df(df_type, df)


    def load_somatic_mutation(self):
        df_type = 'somatic_mutation'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep='\t')

            if self.version == "0.1":
                #Rename the columns we want to keep to the appropriate names
                df = df.rename(columns={"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol_Annovar":"Gene","Variant_Classification_Annovar":"Mutation"})
                df['Location'] = df['Annovar_Info_protein'].str.extract(r'([^:]+$)') #The location that we care about is stored after the last colon
                df = df[['Patient_ID', 'Gene', 'Mutation', 'Location']]

            elif self.version == "2.0":
                df = df[['Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification','HGVSp_Short']]
                df = df.rename(columns={
                    "Tumor_Sample_Barcode":"Patient_ID",
                    "Hugo_Symbol":"Gene",
                    "Variant_Classification":"Mutation",
                    "HGVSp_Short":"Location"}) #Rename the columns we want to keep to the appropriate names

            df = df.sort_values(by=["Patient_ID", "Gene"])
            df = df.set_index("Patient_ID")
            df = df.sort_index()
            df.columns.name=None

            # save df in self._data
            self.save_df(df_type, df)


    def load_annotation(self):
        # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
        file_path = self.locate_files('annotation')
        df = pd.read_csv(file_path, sep='\t')

        if self.version == "2.0":
            df = df.set_index('case_id')
        elif self.version == "0.1":
            df = df.set_index('CASE_ID')

        df.columns.name=None
        df.index.name="Patient_ID"

        # Split the clinical data in to clinical data and derived molecular data
        if self.version == "0.1":
            derived_molecular_cols = ['P53GENE_ANALYSIS', 'EGFR_AMP_STATUS']

        elif self.version == "2.0":
            derived_molecular_cols = ['NAT_pathology_review', 'tumor_pathology_review',
               'ESTIMATE_stromal_score', 'ESTIMATE_immune_score', 'stemness_score',
               'mutation_count', 'TP53_mutation', 'CDKN2A_mutation', 'FAT1_mutation',
               'NOTCH1_mutation', 'CSMD3_mutation', 'DNAH5_mutation', 'KMT2D_mutation',
               'transcriptomic_subtype', 'chr_instability_idx', 'tumor_proportion',
               'normal_epithelial_proportion', 'immune_proportion',
               'muscle_proportion', 'fibroblast_proportion', 'EGFR_pathway',
               'Hypoxia_pathway', 'JAK.STAT_pathway', 'MAPK_pathway', 'NFkB_pathway',
               'PI3K_pathway', 'TGFb_pathway', 'TNFa_pathway', 'Trail_pathway',
               'VEGF_pathway', 'p53_pathway']

        derived_molecular_df = df[derived_molecular_cols]
        derived_molecular_df = derived_molecular_df.sort_index(axis='columns')
        derived_molecular_df = derived_molecular_df.sort_index()

        df = df.drop(columns=derived_molecular_cols)
        df = df.sort_index()
        df = df.sort_index(axis='columns')

        # save dfs in self._data
        self.save_df("clinical", df)
        self.save_df("derived_molecular", derived_molecular_df)

    # load_annotation takes care of the data for these datatypes
    def load_clinical(self):
        if 'clinical' not in self._data:
            self.load_annotation()

    def load_derived_molecular(self):
        if 'derived_molecular' not in self._data:
            self.load_annotation()


    def how_to_cite(self):
        return super().how_to_cite(cancer_type='head and neck squamous cell carcinoma', pmid=33417831)