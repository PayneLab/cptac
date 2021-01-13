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
from .exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError

class Gbm(Dataset):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the gbm dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        valid_versions = ["1.0", "2.0", "2.1", "3.0"]

        data_files = {
            "1.0": [
                "clinical_data_core.v1.0.20190802.tsv.gz",
                "mirnaseq_mirna_mature_tpm.v1.0.20190802.tsv.gz",
                "phosphoproteome_pnnl_d6.v1.0.20190802.tsv.gz",
                "proteome_pnnl_per_gene_d4.v1.0.20190802.tsv.gz",
                "proteome_tmt_design.v1.0.20190802.tsv.gz",
                "rnaseq_gdc_fpkm_uq.v1.0.20190802.tsv.gz",
                "tindaisy_all_cases_filtered.v1.0.20190802.maf.gz",
                "wgs_somatic_cnv_per_gene.v1.0.20190802.tsv.gz"],
            "2.0": [
                "acetylome_pnnl_d6.v2.0.20190905.tsv.gz",
                "clinical_data_core.v2.0.20190905.tsv.gz",
                "metabolome_pnnl.v2.0.20190905.tsv.gz",
                "metabolome_sample_info.v2.0.20190905.tsv.gz",
                "mirnaseq_mirna_mature_tpm.v2.0.20190905.tsv.gz",
                "negative_lipidome_pnnl.v2.0.20190905.tsv.gz",
                "phosphoproteome_pnnl_d6.v2.0.20190905.tsv.gz",
                "positive_lipidome_pnnl.v2.0.20190905.tsv.gz",
                "proteome_pnnl_per_gene_d4.v2.0.20190905.tsv.gz",
                "proteome_tmt_design.v2.0.20190905.tsv.gz",
                "rnaseq_bcm_circular_rna_expression_rsem_uq.v2.0.20190905.tsv.gz",
                "rnaseq_gene_fusion.v2.0.20190905.tsv.gz",
                "rnaseq_washu_fpkm_uq.v2.0.20190905.tsv.gz",
                "tindaisy_all_cases_filtered.v2.0.20190905.maf.gz",
                "wgs_somatic_cnv_per_gene.v2.0.20190905.tsv.gz"],
            "2.1": [
                "acetylome_mssm_per_gene_clean.v2.1.20190927.tsv.gz",
                "clinical_data_core.v2.1.20190927.tsv.gz",
                "metabolome_pnnl.v2.1.20190927.tsv.gz",
                "metabolome_sample_info.v2.1.20190927.tsv.gz",
                "mirnaseq_mirna_mature_tpm.v2.1.20190927.tsv.gz",
                "negative_lipidome_pnnl.v2.1.20190927.tsv.gz",
                "phosphoproteome_mssm_per_gene_clean.v2.1.20190927.tsv.gz",
                "positive_lipidome_pnnl.v2.1.20190927.tsv.gz",
                "proteome_mssm_per_gene_clean.v2.1.20190927.tsv.gz",
                "proteome_tmt_design.v2.1.20190927.tsv.gz",
                "rnaseq_bcm_circular_rna_expression_rsem_uq.v2.1.20190927.tsv.gz",
                "rnaseq_gene_fusion.v2.1.20190927.tsv.gz",
                "rnaseq_washu_fpkm_uq.v2.1.20190927.tsv.gz",
                "tindaisy_all_cases_filtered.v2.1.20190927.maf.gz",
                "wgs_somatic_cnv_per_gene.v2.1.20190927.tsv.gz"],
            "3.0": [
                "acetylome_mssm_per_gene_clean.v3.0.20191121.tsv.gz",
                "clinical_data_core.v3.0.20191121.tsv.gz",
                "gbm_all_subtype_collections.2020-01-13.tsv.gz",
                "metabolome_pnnl.v3.0.20191121.tsv.gz",
                "metabolome_sample_info.v3.0.20191121.tsv.gz",
                "mirnaseq_mirna_mature_tpm.v3.0.20191121.tsv.gz",
                "negative_lipidome_pnnl.v3.0.20191121.tsv.gz",
                "phosphoproteome_mssm_per_gene_clean.v3.0.20191121.tsv.gz",
                "positive_lipidome_pnnl.v3.0.20191121.tsv.gz",
                "proteome_mssm_per_gene_clean.v3.0.20191121.tsv.gz",
                "proteome_tmt_design.v3.0.20191121.tsv.gz",
                "rnaseq_bcm_circular_rna_expression_rsem_uq.v3.0.20191121.tsv.gz",
                "rnaseq_gene_fusion.v3.0.20191121.tsv.gz",
                "rnaseq_washu_fpkm_uq.v3.0.20191121.tsv.gz",
                "tindaisy_all_cases_filtered.v3.0.20191121.maf.gz",
                "wgs_somatic_cnv_per_gene.v3.0.20191121.tsv.gz"],
        }

        super().__init__(cancer_type="gbm", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name, so we don't include the version

            if df_name in ("acetylome_pnnl_d6", "acetylome_mssm_per_gene_clean"):
                df = pd.read_csv(file_path, sep='\t')
                split_genes = df["site"].str.rsplit("-", n=1,expand=True)  # Split the genes from the sites, splitting from the right since some genes have hyphens in their names, but the genes and sites are also separated by hyphens
                df = df.drop(columns="site")
                df = df.assign(Site=split_genes[1])
                df["Site"] = df["Site"].str.replace(r"k", r"", regex=True)  # Get rid of all lowercase k delimeters in the sites

                # Create the multiindex
                df = df.rename(columns={
                        "gene": "Name",
                        "peptide": "Peptide",
                        "refseq_id": "Database_ID",
                    })
                df = df.set_index(["Name", "Site", "Peptide", "Database_ID"])  # Turn these columns into a multiindex
                df = df.sort_index()

                df = df.transpose()
                self._data["acetylproteomics"] = df

            elif df_name == "clinical_data_core":
                df = pd.read_csv(file_path, sep='\t', index_col=0).\
                    assign(Stage="IV") # By definition they're all stage IV, since it's glioblastoma
                self._data["clinical"] = df

            elif file_name == "gbm_all_subtype_collections.2020-01-13.tsv.gz":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.drop(columns="sample_type")
                self._data["derived_molecular"] = df

            elif df_name == "metabolome_pnnl":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                self._data["metabolomics"] = df

            elif df_name == "metabolome_sample_info":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df.index.name = "Patient_ID"
                self._data["sample_info"] = df

            elif df_name == "mirnaseq_mirna_mature_tpm":
                df = pd.read_csv(file_path, sep='\t')
                df = df.rename(columns={"name": "Name", "unique_id": "Database_ID"})
                df = df.set_index(["Name", "Database_ID"]) # We use a multiindex with database IDs, not just names, to avoid duplicate column headers
                df = df.drop(columns=["chromosome", "start", "end", "strand", "mirna_type", "mirbase_id", "precursor_id"])
                df = df.sort_index()
                df = df.transpose()
                self._data["miRNA"] = df

            elif df_name == "negative_lipidome_pnnl":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df = df.add_suffix("_negative")
                self._data["lipidomics_negative"] = df

            elif df_name in ("phosphoproteome_pnnl_d6", "phosphoproteome_mssm_per_gene_clean"):
                df = pd.read_csv(file_path, sep='\t')

                # Create our multiindex
                split_genes = df["site"].str.rsplit("-", n=1, expand=True) # Split the genes from the sites, splitting from the right since some genes have hyphens in their names, but the genes and sites are also separated by hyphens
                df = df.drop(columns="site")
                df = df.assign(Site=split_genes[1])
                df["Site"] = df["Site"].str.replace(r"[sty]", r"", regex=True) # Get rid of all lowercase s, t, and y delimeters in the sites

                if self._version == "1.0":
                    df = df.rename(columns={"gene": "Name", "peptide": "Peptide"})
                    df = df.set_index(["Name", "Site", "Peptide"]) # Turn these columns into a multiindex

                elif self._version in ("2.0", "2.1", "3.0"):
                    df = df.rename(columns={
                            "gene": "Name",
                            "peptide": "Peptide",
                            "refseq_id": "Database_ID",
                        })
                    df = df.set_index(["Name", "Site", "Peptide", "Database_ID"]) # Turn these columns into a multiindex

                df = df.sort_index()
                df = df.transpose()
                self._data["phosphoproteomics"] = df

            elif df_name == "positive_lipidome_pnnl":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.transpose()
                df = df.add_suffix("_positive")
                self._data["lipidomics_positive"] = df

            elif df_name in ("proteome_pnnl_per_gene_d4", "proteome_mssm_per_gene_clean"):
                df = pd.read_csv(file_path, sep='\t', index_col=0)

                if self._version in ("2.0", "2.1", "3.0"):
                    df = df.drop(columns="refseq_id") # We don't need this database ID, because the gene name index is already unique

                df = df.sort_index()
                df = df.transpose()
                self._data["proteomics"] = df

            elif df_name == "proteome_tmt_design":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df.index.name = "Patient_ID"
                self._data["experimental_design"] = df

            elif df_name == "rnaseq_bcm_circular_rna_expression_rsem_uq":
                df = pd.read_csv(file_path, sep='\t')
                df["circRNA_id"] = df["circRNA_id"].str.split('_', n=1, expand=True)[1] # Drop the "circ_" prefix on all the keys
                df = df.set_index("circRNA_id")
                df = df.drop(columns=["gene_id", "gene_name", "gene_type", "alias"])
                df = df.transpose()
                self._data["circular_RNA"] = df

            elif df_name == "rnaseq_gene_fusion":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                self._data["gene_fusion"] = df

            elif df_name in ("rnaseq_gdc_fpkm_uq", "rnaseq_washu_fpkm_uq"):
                df = pd.read_csv(file_path, sep='\t')
                df = df.rename(columns={"gene_name": "Name", "gene_id": "Database_ID"})
                df = df.set_index(["Name", "Database_ID"]) # We use a multiindex with Ensembl IDs, not just gene names, to avoid duplicate column headers
                df = df.drop(columns=["gene_type", "gene_status", "havana_gene", "full_length", "exon_length", "exon_num"])
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                self._data["transcriptomics"] = df

            elif df_name == "tindaisy_all_cases_filtered":
                df = pd.read_csv(file_path, sep='\t')
                split_barcode = df["Tumor_Sample_Barcode"].str.split("_", n = 1, expand = True) # The first part of the barcode is the patient id, which we need want to make the index
                df["Tumor_Sample_Barcode"] = split_barcode[0]
                df = df[["Tumor_Sample_Barcode","Hugo_Symbol","Variant_Classification","HGVSp_Short"]]
                df = df.rename({"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol":"Gene","Variant_Classification":"Mutation","HGVSp_Short":"Location"}, axis='columns')
                df = df.sort_values(by=["Patient_ID", "Gene"])
                df = df.set_index("Patient_ID")
                self._data["somatic_mutation"] = df

            elif df_name == "wgs_somatic_cnv_per_gene":
                df = pd.read_csv(file_path, sep='\t', index_col=0)
                df = df.drop(columns=["gene_id", "gene_id_version", "original_symbol"])
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                self._data["CNV"] = df

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        if self._version in ("2.0", "2.1", "3.0"):
            # Combine positive and negative lipidomics tables
            lipidomics_positive = self._data["lipidomics_positive"]
            lipidomics_negative = self._data["lipidomics_negative"]

            lipidomics = lipidomics_positive.join(lipidomics_negative, how="outer")
            self._data["lipidomics"] = lipidomics

            del self._data["lipidomics_positive"]
            del self._data["lipidomics_negative"]

            # Add useful columns from sample_info table to experimental_design
            sample_info = self._data["sample_info"]
            useful_cols = sample_info[["mass_mg", "is_oct"]]
            useful_cols = useful_cols.add_prefix("sample_")

            experimental_design = self._data["experimental_design"]
            experimental_design = experimental_design.join(useful_cols, how="outer")

            self._data["experimental_design"] = experimental_design
            del self._data["sample_info"]

        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data)

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
        clinical = self._data["clinical"]
        clinical = clinical.reindex(master_index)

        # Construct the sample status column
        sample_status_col = np.where(clinical.index.str.startswith("PT"), "Normal", "Tumor")
        clinical.insert(0, "Sample_Tumor_Normal", sample_status_col)

        # For versions 1.0, 2.0, and 2.1, the gender is mis-entered for two samples in the clinical dataframe. Both C3N-01196 and C3N-01856 are entered as Female, but are actually Male. Let's fix that, if we're loading one of those versions.
        if self._version in ("1.0", "2.0", "2.1"):
            clinical.loc[clinical.index.isin(["C3N-01196", "C3N-01856"]), "gender"] = "Male"

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = clinical

        # Append a ".N" to the Patient_IDs of the normal samples, to match the other datasets
        self._data = reformat_normal_patient_ids(self._data)

        # Call function from dataframe_tools.py to sort all tables first by sample status, and then by the index
        self._data = sort_all_rows(self._data)

        # Call function from dataframe_tools.py to standardize the names of the index and column axes
        self._data = standardize_axes_names(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message

        # Print data embargo warning, if the date hasn't passed yet.
        today = datetime.date.today()
        embargo_date = datetime.date(year=2021, month=3, day=1)
        if today < embargo_date:
            warnings.warn("The GBM dataset is under publication embargo until March 01, 2021. CPTAC is a community resource project and data are made available rapidly after generation for community research use. The embargo allows exploring and utilizing the data, but analysis may not be published until after the embargo date. Please see https://proteomics.cancer.gov/data-portal/about/data-use-agreement or enter cptac.embargo() to open the webpage for more details.", PublicationEmbargoWarning, stacklevel=2)
