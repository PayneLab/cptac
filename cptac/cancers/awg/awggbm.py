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

from asyncore import file_dispatcher
import pandas as pd
import numpy as np
import os
import warnings
import datetime

from cptac.cancers.source import Source
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError

class AwgGbm(Source):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the gbm dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        self.valid_versions = ["1.0", "2.0", "2.1", "3.0"]

        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        self.data_files = {
            "1.0": {
                "clinical"              : "clinical_data_core.v1.0.20190802.tsv.gz",
                "miRNA"                 : "mirnaseq_mirna_mature_tpm.v1.0.20190802.tsv.gz",
                "phosphoproteomics"     : "phosphoproteome_pnnl_d6.v1.0.20190802.tsv.gz",
                "proteomics"            : "proteome_pnnl_per_gene_d4.v1.0.20190802.tsv.gz",
                "experimental_design"   : "proteome_tmt_design.v1.0.20190802.tsv.gz",
                "transcriptomics"       : "rnaseq_gdc_fpkm_uq.v1.0.20190802.tsv.gz",
                "somatic_mutation"      : "tindaisy_all_cases_filtered.v1.0.20190802.maf.gz",
                "CNV"                   : "wgs_somatic_cnv_per_gene.v1.0.20190802.tsv.gz"},
            "2.0": {
                "acetylproteomics"      : "acetylome_pnnl_d6.v2.0.20190905.tsv.gz",
                "clinical"              : "clinical_data_core.v2.0.20190905.tsv.gz",
                "metabolomics"          : "metabolome_pnnl.v2.0.20190905.tsv.gz",
                "experimental_design"   : ["metabolome_sample_info.v2.0.20190905.tsv.gz", "proteome_tmt_design.v2.0.20190905.tsv.gz"],
                "miRNA"                 : "mirnaseq_mirna_mature_tpm.v2.0.20190905.tsv.gz",
                "lipidomics"            : [ "negative_lipidome_pnnl.v2.0.20190905.tsv.gz", "positive_lipidome_pnnl.v2.0.20190905.tsv.gz"],
                "phosphoproteomics"     : "phosphoproteome_pnnl_d6.v2.0.20190905.tsv.gz",
                "proteomics"            : "proteome_pnnl_per_gene_d4.v2.0.20190905.tsv.gz",
                "circular_RNA"          : "rnaseq_bcm_circular_rna_expression_rsem_uq.v2.0.20190905.tsv.gz",
                "gene_fusion"           : "rnaseq_gene_fusion.v2.0.20190905.tsv.gz",
                "transcriptomics"       : "rnaseq_washu_fpkm_uq.v2.0.20190905.tsv.gz",
                "somatic_mutation"      : "tindaisy_all_cases_filtered.v2.0.20190905.maf.gz",
                "CNV"                   : "wgs_somatic_cnv_per_gene.v2.0.20190905.tsv.gz"},
            "2.1": {
                "acetylproteomics"      : "acetylome_mssm_per_gene_clean.v2.1.20190927.tsv.gz",
                "clinical"              : "clinical_data_core.v2.1.20190927.tsv.gz",
                "metabolomics"          : "metabolome_pnnl.v2.1.20190927.tsv.gz",
                "experimental_design"   : ["metabolome_sample_info.v2.1.20190927.tsv.gz", "proteome_tmt_design.v2.1.20190927.tsv.gz"],
                "miRNA"                 : "mirnaseq_mirna_mature_tpm.v2.1.20190927.tsv.gz",
                "lipidomics"            : ["negative_lipidome_pnnl.v2.1.20190927.tsv.gz", "positive_lipidome_pnnl.v2.1.20190927.tsv.gz"],
                "phosphoproteomics"     : "phosphoproteome_mssm_per_gene_clean.v2.1.20190927.tsv.gz",
                "proteomics"            : "proteome_mssm_per_gene_clean.v2.1.20190927.tsv.gz",
                "circular_RNA"          : "rnaseq_bcm_circular_rna_expression_rsem_uq.v2.1.20190927.tsv.gz",
                "gene_fusion"           : "rnaseq_gene_fusion.v2.1.20190927.tsv.gz",
                "transcriptomics"       : "rnaseq_washu_fpkm_uq.v2.1.20190927.tsv.gz",
                "somatic_mutation"      : "tindaisy_all_cases_filtered.v2.1.20190927.maf.gz",
                "CNV"                   : "wgs_somatic_cnv_per_gene.v2.1.20190927.tsv.gz"},
            "3.0": {
                "acetylproteomics"      : "acetylome_mssm_per_gene_clean.v3.0.20191121.tsv.gz",
                "clinical"              : "clinical_data_core.v3.0.20191121.tsv.gz",
                "derived_molecular"     : "gbm_all_subtype_collections.2020-01-13.tsv.gz",
                "metabolomics"          : "metabolome_pnnl.v3.0.20191121.tsv.gz",
                "experimental_design"   : ["metabolome_sample_info.v3.0.20191121.tsv.gz", "proteome_tmt_design.v3.0.20191121.tsv.gz"],
                "miRNA"                 : "mirnaseq_mirna_mature_tpm.v3.0.20191121.tsv.gz",
                "lipidomics"            : ["negative_lipidome_pnnl.v3.0.20191121.tsv.gz", "positive_lipidome_pnnl.v3.0.20191121.tsv.gz"],
                "phosphoproteomics"     : "phosphoproteome_mssm_per_gene_clean.v3.0.20191121.tsv.gz",
                "proteomics"            : "proteome_mssm_per_gene_clean.v3.0.20191121.tsv.gz",
                "circular_RNA"          : "rnaseq_bcm_circular_rna_expression_rsem_uq.v3.0.20191121.tsv.gz",
                "gene_fusion"           : "rnaseq_gene_fusion.v3.0.20191121.tsv.gz",
                "transcriptomics"       : "rnaseq_washu_fpkm_uq.v3.0.20191121.tsv.gz",
                "somatic_mutation"      : "tindaisy_all_cases_filtered.v3.0.20191121.maf.gz",
                "CNV"                   : "wgs_somatic_cnv_per_gene.v3.0.20191121.tsv.gz"},
        }

        self.load_functions = {
            "clinical"               : self.load_clinical,
            "CNV"                    : self.load_CNV,
            "experimental_design"    : self.load_experimental_design,
            "miRNA"                  : self.load_miRNA,
            "phosphoproteomics"      : self.load_phosphoproteomics,
            "proteomics"             : self.load_proteomics,
            "somatic_mutation"       : self.load_somatic_mutation,
            "transcriptomics"        : self.load_transcriptomics,
        }

        if version != "1.0":
            self.load_functions["acetylproteomics"] = self.load_acetylproteomics
            self.load_functions["metabolomics"] = self.load_metabolomics
            self.load_functions["lipidomics"] = self.load_lipidomics
            self.load_functions["circular_RNA"] = self.load_circular_RNA
            self.load_functions["gene_fusion"] = self.load_gene_fusion

        super().__init__(cancer_type="gbm", source='awg', version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)


    def how_to_cite(self):
        return super().how_to_cite(cancer_type='glioblastoma', pmid=33577785)

    def load_acetylproteomics(self):
        df_type = 'acetylproteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            # Split the genes from the sites, splitting from the right since some genes have hyphens in their names, but the genes and sites are also separated by hyphens
            split_genes = df["site"].str.rsplit("-", n=1,expand=True)
            df = df.drop(columns="site")
            df = df.assign(Site=split_genes[1])
            # Get rid of all lowercase k delimeters in the sites
            df["Site"] = df["Site"].str.replace(r"k", r"", regex=True)

            # Create the multiindex
            df = df.rename(columns={
                    "gene": "Name",
                    "peptide": "Peptide",
                    "refseq_id": "Database_ID",
                })
            df = df.set_index(["Name", "Site", "Peptide", "Database_ID"])  # Turn these columns into a multiindex
            df = df.sort_index()
            df = df.transpose()

            # save df in self._data
            self.save_df(df_type, df)

    def load_circular_RNA(self):
        df_type = 'circular_RNA'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df["circRNA_id"] = df["circRNA_id"].str.split('_', n=1, expand=True)[1] # Drop the "circ_" prefix on all the keys
            df = df.set_index("circRNA_id")
            df = df.drop(columns=["gene_id", "gene_name", "gene_type", "alias"])
            df = df.transpose()

            # save df in self._data
            self.save_df(df_type, df)

    def load_clinical(self):
        df_type = 'clinical'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0).assign(Stage="IV") # By definition they're all stage IV, since it's glioblastoma

            # For versions 1.0, 2.0, and 2.1, the gender is mis-entered for two samples in the clinical dataframe.
            # Both C3N-01196 and C3N-01856 are entered as Female, but are actually Male. Let's fix that.
            if self.version in ("1.0", "2.0", "2.1"):
                df.loc[df.index.isin(["C3N-01196", "C3N-01856"]), "gender"] = "Male"

            # save df in self._data
            self.save_df(df_type, df)

    def load_CNV(self):
        df_type = 'CNV'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.drop(columns=["gene_id", "gene_id_version", "original_symbol"])
            df = df.sort_index()
            df = df.transpose()
            df = df.sort_index()

            # save df in self._data
            self.save_df(df_type, df)

    def load_derived_molecular(self):
        df_type = 'derived_molecular'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.drop(columns="sample_type")

            # save df in self._data
            self.save_df(df_type, df)

    def load_experimental_design(self):
        df_type = 'experimental_design'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_paths = self.locate_files(df_type)

            for file_path in file_paths:
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df.index.name = "Patient_ID"
                if "metabolome_sample_info" in file_path:
                    sample_info = df
                else:
                    experimental_design = df

            # Add useful columns from sample_info table to experimental_design
            useful_cols = sample_info[["mass_mg", "is_oct"]]
            useful_cols = useful_cols.add_prefix("sample_")
            experimental_design = experimental_design.join(useful_cols, how="outer")

            # save df in self._data
            self.save_df(df_type, experimental_design)

    def load_gene_fusion(self):
        df_type = 'gene_fusion'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)

            # save df in self._data
            self.save_df(df_type, df)

    def load_lipidomics(self):
        df_type = 'lipidomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_paths = self.locate_files(df_type)

            for file_path in file_paths:
                if "positive_lipidome_pnnl" in file_path:
                    df = pd.read_csv(file_path, sep='\t', index_col=0)
                    df = df.transpose()
                    df = df.add_suffix("_positive")
                    lipidomics_positive = df

                elif "negative_lipidome_pnnl" in file_path:
                    df = pd.read_csv(file_path, sep='\t', index_col=0)
                    df = df.transpose()
                    df = df.add_suffix("_negative")
                    lipidomics_negative= df

            # Combine positive and negative lipidomics tables
            lipidomics = lipidomics_positive.join(lipidomics_negative, how="outer")

            # save df in self._data
            self.save_df(df_type, lipidomics)

    def load_metabolomics(self):
        df_type = 'metabolomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df = df.transpose()

            # save df in self._data
            self.save_df(df_type, df)

    def load_miRNA(self):
        df_type = 'miRNA'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df = df.rename(columns={"name": "Name", "unique_id": "Database_ID"})
            df = df.set_index(["Name", "Database_ID"]) # We use a multiindex with database IDs, not just names, to avoid duplicate column headers
            df = df.drop(columns=["chromosome", "start", "end", "strand", "mirna_type", "mirbase_id", "precursor_id"])
            df = df.sort_index()
            df = df.transpose()

            # save df in self._data
            self.save_df(df_type, df)    

    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')

            # Create our multiindex
            # Split the genes from the sites, splitting from the right since some genes have hyphens in their names, but the genes and sites are also separated by hyphens
            split_genes = df["site"].str.rsplit("-", n=1, expand=True)
            df = df.drop(columns="site")
            df = df.assign(Site=split_genes[1])
            # Get rid of all lowercase s, t, and y delimeters in the sites
            df["Site"] = df["Site"].str.replace(r"[sty]", r"", regex=True)

            if self.version == "1.0":
                df = df.rename(columns={"gene": "Name", "peptide": "Peptide"})
                df = df.set_index(["Name", "Site", "Peptide"]) # Turn these columns into a multiindex

            elif self.version in ("2.0", "2.1", "3.0"):
                df = df.rename(columns={
                        "gene": "Name",
                        "peptide": "Peptide",
                        "refseq_id": "Database_ID",
                    })
                df = df.set_index(["Name", "Site", "Peptide", "Database_ID"]) # Turn these columns into a multiindex

            df = df.sort_index()
            df = df.transpose()

            # save df in self._data
            self.save_df(df_type, df)

    def load_proteomics(self):
        df_type = 'proteomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t', index_col=0)

            if self.version in ("2.0", "2.1", "3.0"):
                df = df.drop(columns="refseq_id") # We don't need this database ID, because the gene name index is already unique

            df = df.sort_index()
            df = df.transpose()

            # save df in self._data
            self.save_df(df_type, df)

    def load_somatic_mutation(self):
        df_type = 'somatic_mutation'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n = 1, expand = True) # The first part of the barcode is the patient id, which we need want to make the index
            df["Tumor_Sample_Barcode"] = split_barcode[0]
            df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
            df = df.rename({"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
            df = df.sort_values(by=["Patient_ID", "Gene"])
            df = df.set_index("Patient_ID")

            # save df in self._data
            self.save_df(df_type, df)

    def load_transcriptomics(self):
        df_type = 'transcriptomics'
        if df_type not in self._data:
            # verify the df_type is valid for the current version and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep='\t')
            df = df.rename(columns={"gene_name": "Name", "gene_id": "Database_ID"})
            df = df.set_index(["Name", "Database_ID"]) # We use a multiindex with Ensembl IDs, not just gene names, to avoid duplicate column headers
            df = df.drop(columns=["gene_type", "gene_status", "havana_gene", "full_length", "exon_length", "exon_num"])
            df = df.sort_index()
            df = df.transpose()
            df = df.sort_index()

            # save df in self._data
            self.save_df(df_type, df)

    # Override the save_df function from source.py so we can mark normal samples
    def save_df(self, datatype, df):
        # Append a ".N" to the Patient_IDs of the normal samples, to match the other datasets
        # In GBM all normal samples begin with "PT"
        df.index = df.index.where(~df.index.str.startswith('PT'), df.index + ".N")

        # Inherit the parent event
        super().save_df(datatype, df)
