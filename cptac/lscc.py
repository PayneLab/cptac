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
class Lscc(DataSet):

    def __init__(self, version="latest"):
        """Load all of the lscc dataframes as values in the self._data dict variable, with names as keys, and format them properly."""

        # Set some needed variables, and pass them to the parent DataSet class __init__ function

        valid_versions = ["1.0"] # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.

        data_files = {
            "1.0": [
                "lscc-v1.0-cnv-gene-level-log2.gct.gz", #done
                "lscc-v1.0-cptac3-lscc-rna-seq-fusion-v2.2-y2.all-20190807.txt.gz", #done
                "lscc-v1.0-cptac3-lscc-wxs-somatic-variant-sw-v1.5-lscc.y2-20191211.maf.gz", #somatic variant, find sample identifier (maybe cl), name, gene, mutation, location
                "lscc-v1.0-mirna-mature-tpm-log2.gct.gz", #micro rna #done
                "lscc-v1.0-phosphoproteome-ratio-norm-NArm.gct.gz", #done
                "lscc-v1.0-proteome-ratio-norm-NArm.gct.gz", #done
                "lscc-v1.0-rnaseq-uq-fpkm-log2-NArm.gct.gz", # done
                "lscc-v1.0-sample-annotation.csv.gz"] #done
        }

        super().__init__(cancer_type="lscc", version=version, valid_versions=valid_versions, data_files=data_files)

        # Load the data into dataframes in the self._data dict
        loading_msg = "Loading dataframes"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            if file_name == "lscc-v1.0-cnv-gene-level-log2.gct.gz": #Done
                df = pd.read_csv(file_path, sep="\t", skiprows=2, dtype=object)
                df = df.set_index("id")
                # df = df.apply(pd.to_numeric)
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name="Patient_ID"
                self._data["CNV"] = df

            elif file_name == "lscc-v1.0-phosphoproteome-ratio-norm-NArm.gct.gz": #Done
                df = pd.read_csv(file_path, sep="\t", skiprows=2, dtype=object)
                gene_filter = df['geneSymbol'] != 'na' #Drop rows of metadata
                df = df[gene_filter]

                # Prepare some columns we'll need later for the multiindex
                df["variableSites"] = df["variableSites"].str.replace(r"[a-z\s]", "") # Get rid of all lowercase delimeters and whitespace in the sites
                df = df.rename(columns={
                "geneSymbol": "Name",
                "variableSites": "Site",
                "sequence": "Peptide", # We take this instead of sequenceVML, to match the other datasets' format
                "accession_numbers": "Database_ID" # We take all accession numbers they have, instead of the singular accession_number column
                })

                # Some rows have at least one localized phosphorylation site, but also have other phosphorylations that aren't localized. We'll drop those rows, if their localized sites are duplicated in another row, to avoid creating duplicates, because we only preserve information about the localized sites in a given row. However, if the localized sites aren't duplicated in another row, we'll keep the row.
                unlocalized_to_drop = df.index[~df['Best_numActualVMSites_sty'].eq(df['Best_numLocalizedVMsites_sty']) & df.duplicated(["Name", "Site", "Peptide", "Database_ID"], keep=False)] # Column 3 of the split "id" column is number of phosphorylations detected, and column 4 is number of phosphorylations localized, so if the two values aren't equal, the row has at least one unlocalized site
                df = df.drop(index=unlocalized_to_drop)

                # Give it a multiindex
                df = df.set_index(["Name", "Site", "Peptide", "Database_ID"])


                cols_to_drop = ['id','id.description', 'numColumnsVMsiteObserved', 'bestScore', 'bestDeltaForwardReverseScore',
                'Best_scoreVML', 'Best_numActualVMSites_sty', 'Best_numLocalizedVMsites_sty', 'sequenceVML',
                'accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA', 'protein_mw', 'species',
                'speciesMulti', 'orfCategory', 'accession_number', 'protein_group_num', 'entry_name', 'GeneSymbol']
                df = df.drop(columns=cols_to_drop)
                df = df.apply(pd.to_numeric)
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name="Patient_ID"
                self._data["phosphoproteomics"] = df



            elif file_name == "lscc-v1.0-proteome-ratio-norm-NArm.gct.gz": #done
                df = pd.read_csv(file_path, skiprows=2, sep='\t', dtype=object)
                gene_filter = df['geneSymbol'] != 'na' #Filter out rows of metadata
                df = df[gene_filter]

                df = df.rename(columns={"GeneSymbol": "Name", 'accession_numbers': "Database_ID"})
                df = df.set_index(["Name", "Database_ID"])
                cols_to_drop = ['id', 'id.description', 'geneSymbol', 'numColumnsProteinObserved',
                'numSpectraProteinObserved', 'protein_mw', 'percentCoverage', 'numPepsUnique',
                'scoreUnique', 'species', 'orfCategory', 'accession_number',
                'subgroupNum', 'entry_name']
                df = df.drop(columns=cols_to_drop)
                df = df.apply(pd.to_numeric)
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name="Patient_ID"
                df.columns.name=None
                self._data["proteomics"] = df


            elif file_name == "lscc-v1.0-cptac3-lscc-rna-seq-fusion-v2.2-y2.all-20190807.txt.gz": #done
                 df = pd.read_csv(file_path, sep="\t", dtype=object)
                 self._data['gene_fusion'] = df


            elif file_name == "lscc-v1.0-sample-annotation.csv.gz": #done
                df = pd.read_csv(file_path, sep=",", dtype=object)
                filter = df['QC.status'] == "QC.pass" #There are some samples that are internal references. IRs are used for scaling purposes, and don't belong to a single patient, so we want to drop them.
                df = df[filter]
                df = df.drop(columns="Participant") #Get rid of the "Participant" column becuase the same information is stored in Sample.ID  which is formatted the way we want.
                df = df.set_index("Sample.ID")
                df.index.name="Patient_ID"
                df = df.rename(columns={"Type":"Sample_Tumor_Normal"})
                #Split the metadata into multiple dataframes
                #Make experiemntal_set up dataframe
                experimental_design_cols = ['Experiment', 'Channel', 'QC.status'] #These are the columns for the experimental_design dataframe
                experimental_design_df = df[experimental_design_cols]
                df = df.drop(columns=experimental_design_cols)
                #Make a derived_molecular dataframe
                derived_molecular_cols = ['TP53.mutation', 'CDKN2A.mutation', 'PTEN.mutation', 'PIK3CA.mutation',
                 'KEAP1.mutation', 'HLA.A.mutation', 'NFE2L2.mutation', 'NOTCH1.mutation', 'RB1.mutation',
                 'HRAS.mutation', 'FBXW7.mutation', 'SMARCA4.mutation', 'NF1.mutation', 'SMAD4.mutation',
                 'EGFR.mutation', 'APC.mutation', 'BRAF.mutation', 'TNFAIP3.mutation', 'CREBBP.mutation',
                 'TP53.mutation.status', 'CDKN2A.mutation.status', 'PTEN.mutation.status', 'PIK3CA.mutation.status',
                 'KEAP1.mutation.status', 'HLA.A.mutation.status', 'NFE2L2.mutation.status', 'NOTCH1.mutation.status',
                 'RB1.mutation.status', 'HRAS.mutation.status', 'FBXW7.mutation.status', 'SMARCA4.mutation.status',
                 'NF1.mutation.status', 'SMAD4.mutation.status', 'EGFR.mutation.status', 'APC.mutation.status',
                 'BRAF.mutation.status', 'TNFAIP3.mutation.status', 'CREBBP.mutation.status']
                derived_molecular_df = df[derived_molecular_cols]
                df = df.drop(columns = derived_molecular_cols)
                self._data["clinical"]= df
                self._data['experimental_design'] = experimental_design_df
                self._data['derived_molecular'] = derived_molecular_df



            elif file_name == "lscc-v1.0-cptac3-lscc-wxs-somatic-variant-sw-v1.5-lscc.y2-20191211.maf.gz":
                df = pd.read_csv(file_path, sep="\t", dtype=object)

                cols_to_drop = ['Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position',
                'Strand', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS',
                'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1',
                'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1',
                'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase',
                'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Sequencer', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID',
                'HGVSc', 'HGVSp', 'Transcript_ID', 'Exon_Number', 't_depth', 't_ref_count', 't_alt_count', 'n_depth',
                'n_ref_count', 'n_alt_count', 'callers', 'all_effects', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence',
                'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM',
                'DISTANCE', 'STRAND_VEP', 'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE', 'CANONICAL', 'CCDS', 'ENSP', 'SWISSPROT',
                'TREMBL', 'UNIPARC', 'RefSeq', 'SIFT', 'PolyPhen', 'EXON', 'INTRON', 'DOMAINS', 'GMAF', 'AFR_MAF', 'AMR_MAF',
                'ASN_MAF', 'EAS_MAF', 'EUR_MAF', 'SAS_MAF', 'AA_MAF', 'EA_MAF', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'MOTIF_NAME',
                'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'IMPACT', 'PICK', 'VARIANT_CLASS', 'TSL', 'HGVS_OFFSET', 'PHENO',
                'MINIMISED', 'ExAC_AF', 'ExAC_AF_AFR', 'ExAC_AF_AMR', 'ExAC_AF_EAS', 'ExAC_AF_FIN', 'ExAC_AF_NFE', 'ExAC_AF_OTH',
                'ExAC_AF_SAS', 'GENE_PHENO', 'FILTER', 'flanking_bps', 'variant_id', 'variant_qual', 'ExAC_AF_Adj', 'ExAC_AC_AN_Adj',
                'ExAC_AC_AN', 'ExAC_AC_AN_AFR', 'ExAC_AC_AN_AMR', 'ExAC_AC_AN_EAS', 'ExAC_AC_AN_FIN', 'ExAC_AC_AN_NFE',
                'ExAC_AC_AN_OTH', 'ExAC_AC_AN_SAS', 'ExAC_FILTER']

                df = df.drop(columns = cols_to_drop)
                df = df.rename(columns={"Sample.ID": "Patient_ID", 'Hugo_Symbol': "Gene", "Variant_Classification": "Mutation", "HGVSp_Short": "Location"})
                df.set_index("Patient_ID")
                df.sort_index()
                self._data['somatic_mutation'] = df

            elif file_name == "lscc-v1.0-mirna-mature-tpm-log2.gct.gz": #Done
                df = pd.read_csv(file_path, skiprows=2, sep='\t', dtype=object)
                gene_filter = df['Name'] != 'na' #Filter out rows of metadata
                df = df[gene_filter]
                df = df.set_index(["Name","ID"])
                #df= df.apply(pd.to_numeric)
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name="Patient_ID"
                df.columns.name=None
                self._data["miRNA"] = df

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data)

        # Sort this master_index so all the samples with an N suffix are last. Because the N is a suffix, not a prefix, this is kind of messy.
        status_col = np.where(master_index.str.endswith("N"), "Normal", "Tumor")
        status_df = pd.DataFrame(data={"Patient_ID": master_index, "Status": status_col}) # Create a new dataframe with the master_index as a column called "Patient_ID"
        status_df = status_df.sort_values(by=["Status", "Patient_ID"], ascending=[False, True]) # Sorts first by status, and in descending order, so "Tumor" samples are first
        master_index = pd.Index(status_df["Patient_ID"])

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
        clinical = self._data["clinical"]
        clinical = clinical.reindex(master_index)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = clinical

        # Replace periods with hyphens in all Patient_IDs
        for name in self._data.keys(): # Loop over just the keys to avoid any issues that would come if we looped over the values while editing them
            df = self._data[name]

            # Make the index a column so we can edit it
            df.index.name = "Patient_ID"
            df = df.reset_index()

            # Replace all '.' with '-'
            df["Patient_ID"] = df["Patient_ID"].str.replace(r"\.", "-")
            df["Patient_ID"] = df["Patient_ID"].str.replace(r"-N$", ".N") # If there's a "-N" at the end, it's part of the normal identifier, which we want to actually be ".N"

            # Set the index back to Patient_ID
            df = df.set_index("Patient_ID")

            self._data[name] = df

        # Call function from dataframe_tools.py to standardize the names of the index and column axes
        self._data = standardize_axes_names(self._data)

        # Call function from dataframe_tools.py to sort all tables first by sample status, and then by the index
        self._data = sort_all_rows(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
