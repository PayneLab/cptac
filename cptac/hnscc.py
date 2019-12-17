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
from .dataset import DataSet
from .dataframe_tools import *
from .exceptions import FailedReindexWarning, ReindexMapError

class Hnscc(DataSet):

    def __init__(self, version="latest"):
        """Load all of the hnscc dataframes as values in the self._data dict variable, with names as keys, and format them properly."""

        # Set some needed variables, and pass them to the parent DataSet class __init__ function

        valid_versions = ["0.1"] # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.

        data_files = {
            "0.1": [
                "HNSCC.strelka.sorted.filtered.annovar.hg19_multianno_filtered.maf.txt.gz",
                "Proteomics_DIA_Gene_level_Normal.cct.gz",
                "Proteomics_DIA_Gene_level_Tumor.cct.gz",
                "RNAseq_RSEM_UQ_log2.cct.gz",
                "RNAseq_circ_RSEM_UQ_log2.cct.gz",
                "SCNA_gene_level.cct.gz",
                "clinic.tsi.gz"],
            "2.0": [
                "circRNAseq_RSEM_UQ_log2_Combined.cct.gz",
                "Meta_table.tsv.gz",
                "microRNA_log2_Combined.cct.gz",
                "Phosphoproteomics_TMT_site_level_Combined_all.cct.gz",
                "Proteomics_TMT_Gene_level_Combined_all.cct.gz",
                "RNAseq_RSEM_UQ_Combined.cct.gz",
                "SCNA_log2_gene_level.cct.gz",
                "SomaticMutations_maf.tsv.gz"]
        }

        super().__init__(cancer_type="hnscc", version=version, valid_versions=valid_versions, data_files=data_files)

        # Load the data into dataframes in the self._data dict
        loading_msg = "Loading dataframes"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            if file_name == "SCNA_gene_level.cct.gz" or file_name == "SCNA_log2_gene_level.cct.gz":
                df = pd.read_csv(file_path, sep="\t")

                if self._version == "2.0":
                    df = df.set_index('gene_symbol')

                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.columns.name=None
                df.index.name = "Patient_ID"
                self._data["CNV"] = df

            elif file_name == "RNAseq_RSEM_UQ_log2.cct.gz" or file_name == "RNAseq_RSEM_UQ_Combined.cct.gz":
                df = pd.read_csv(file_path, sep="\t")

                if self._version == "2.0":
                    df = df.set_index('Idx')

                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.columns.name=None

                if self._version == "0.1":
                    df.index = df.index.str.replace(r'\.', '-', 1)
                    df.index = df.index.str.replace(r'\.T$', '', 1)
                elif self._version == "2.0":
                        df.index = df.index.str.replace(r'-T$', '', 1)
                        df.index = df.index.str.replace(r'-N$', '.N', 1)

                df.index.name = "Patient_ID"
                self._data["transcriptomics"] = df

            elif file_name == "RNAseq_circ_RSEM_UQ_log2.cct.gz" or file_name == "circRNAseq_RSEM_UQ_log2_Combined.cct.gz":
                df = pd.read_csv(file_path, sep='\t')
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.columns.name=None
                if self._version == "0.1":
                    df.index = df.index.str.replace(r'\.', '-', 1) #We want all the patientIDs to have the the format C3L-00977, and these have the form C3L.00977.N, so we need to replace the first "." with a "-"
                    df.index = df.index.str.replace(r'\.T$', '', 1)
                if self._version == "2.0":
                    df.index = df.index.str.replace(r'-T$', '', 1)
                    df.index = df.index.str.replace(r'-N$', '.N', 1)
                df.index.name = "Patient_ID"
                self._data["circular_RNA"] = df

            elif file_name == "HNSCC.strelka.sorted.filtered.annovar.hg19_multianno_filtered.maf.txt.gz" or file_name == "SomaticMutations_maf.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")

                if self._version == "0.1":
                    df = df.rename(columns={"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol_Annovar":"Gene","Variant_Classification_Annovar":"Mutation"}) #Rename the columns we want to keep to the appropriate names
                    df['Location'] = df['Annovar_Info_protein'].str.extract(r'([^:]+$)') #The location that we care about is stored after the last colon
                    df = df[['Patient_ID', 'Gene', 'Mutation', 'Location']]
                elif self._version == "2.0":
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
                self._data["somatic_mutation"] = df

            elif file_name == "clinic.tsi.gz" or file_name == "Meta_table.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")

                if self._version == "2.0":
                    df = df.set_index('case_id')
                elif self._version == "0.1":
                    df = df.set_index('CASE_ID')

                df.columns.name=None
                df.index.name="Patient_ID"
                #Split the clinicl data in to clincial data and derived molecular data

                if self._version == "0.1":
                    derived_molecular_cols = ['P53GENE_ANALYSIS', 'EGFR_AMP_STATUS']
                elif self._version == "2.0":
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
                derived_molecular_df = df[derived_molecular_cols]#resume here on update
                derived_molecular_df = derived_molecular_df.sort_index(axis='columns')
                derived_molecular_df = derived_molecular_df.sort_index()
                df = df.drop(columns=derived_molecular_cols)
                df = df.sort_index()
                df = df.sort_index(axis='columns')
                self._data["clinical"] = df
                self._data["derived_molecular"] = derived_molecular_df

            #the only files left are the proteomics files
            elif file_name in ["Proteomics_DIA_Gene_level_Normal.cct.gz", "Proteomics_DIA_Gene_level_Tumor.cct.gz", "Proteomics_TMT_Gene_level_Combined_all.cct.gz"]:
                df = pd.read_csv(file_path, sep="\t")

                if self._version == "2.0":
                    df = df.set_index('Index')

                df = df.transpose()
                df.columns.name=None
                df.index.name = "Patient_ID"

                if self._version == "2.0":
                    df.index = df.index.str.replace(r'-T$', '', 1)
                    df.index = df.index.str.replace(r'-N$', '.N', 1)
                    df.index = df.index.str.replace(r'-C$', '.C', 1) #-C is cored NAT samples

                #Once the files are formatted correctly load them into self._data
                if file_name == "Proteomics_DIA_Gene_level_Normal.cct.gz":
                    self._data["proteomics_normal"] = df
                elif file_name == "Proteomics_DIA_Gene_level_Tumor.cct.gz":
                    self._data["proteomics_tumor"] = df
                elif file_name == "Proteomics_TMT_Gene_level_Combined_all.cct.gz":
                    self._data["proteomics"] = df

            elif file_name == "Phosphoproteomics_TMT_site_level_Combined_all.cct.gz" and self._version == "2.0":
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
                df.index = df.index.str.replace(r'-T$', '', 1)
                df.index = df.index.str.replace(r'-N$', '.N', 1)
                df.index = df.index.str.replace(r'-C$', '.C', 1) #-C is cored NAT samples
                df = df.sort_index()
                self._data["phosphoproteomics"] = df

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        if self._version == "0.1":
            # Combine the two proteomics dataframes
            df_normal = self._data.get("proteomics_normal")
            df_tumor = self._data.get("proteomics_tumor")

            df_normal.index = df_normal.index + ".N" #concatenate a ".N" onto the end of the normal data so we can identify it as normal after it's appended to tumor
            prot_combined = df_tumor.append(df_normal) #append the normal data onto the end of the tumor data
            prot_combined = prot_combined.sort_index(axis='columns') # Put all the columns in alphabetical order
            prot_combined = prot_combined.sort_index()
            self._data["proteomics"] = prot_combined
            del self._data["proteomics_normal"]
            del self._data["proteomics_tumor"]

        import pdb; pdb.set_trace()
        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data)

        # Sort this master_index so all the samples with an N suffix are last. Because the N is a suffix, not a prefix, this is kind of messy.
        status_col = np.where(master_index.str.endswith("N"), "Normal", "Tumor")
        status_df = pd.DataFrame(data={"Patient_ID": master_index, "Status": status_col}) # Create a new dataframe with the master_index as a column called "Patient_ID"
        status_df = status_df.sort_values(by=["Status", "Patient_ID"], ascending=[False, True]) # Sorts first by status, and in descending order, so "Tumor" samples are first
        master_index = pd.Index(status_df["Patient_ID"])

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
        master_clinical = self._data['clinical'].reindex(master_index)

        # Add a column called Sample_Tumor_Normal to the clinical dataframe indicating whether each sample is a tumor or normal sample. Samples with a Patient_ID ending in N are normal.
        clinical_status_col = generate_sample_status_col(master_clinical, normal_test=lambda sample: sample[-1] == 'N')
        master_clinical.insert(0, "Sample_Tumor_Normal", clinical_status_col)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = master_clinical

        # Call function from dataframe_tools.py to reindex all the dataframes to have Sample_ID indices
        self._data = reindex_all(self._data, master_index)

        # Now that we've reindexed all the dataframes with sample IDs, edit the format of the Patient_IDs in the clinical dataframe to have normal samples marked the same way as in other datasets. Currently, normal patient IDs have a ".N" appended. We're going to erase that and prepend an "N."
        self._data = reformat_normal_patient_ids(self._data, existing_identifier=".N", existing_identifier_location="end")

        # Call function from dataframe_tools.py to standardize the names of the index and column axes
        self._data = standardize_axes_names(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
