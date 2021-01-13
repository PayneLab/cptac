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

class Luad(Dataset):

    def __init__(self, version="latest", no_internet=False):
        """Load all of the luad dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        valid_versions = ["2.0", "3.1", "3.1.1"] # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.

        data_files = {
            "2.0": [
                "luad-v2.0-cnv-gene-LR.gct.gz",
                "luad-v2.0-phosphoproteome-ratio-norm-NArm.gct.gz",
                "luad-v2.0-proteome-ratio-norm-NArm.gct.gz",
                "luad-v2.0-rnaseq-circ-rna.csv.gz",
                "luad-v2.0-rnaseq-prot-uq-rpkm-log2-NArm-row-norm.gct.gz",
                "luad-v2.0-sample-annotation.csv.gz"],
            "3.1": [
                "LUAD_followup_9_12.xlsx",
                "luad-v3.0-rnaseq-circ-rna.csv.gz",
                "luad-v3.0-rnaseq-gene-fusions.csv.gz",
                "luad-v3.0-wxs-somatic.luad.v1.4.20190517.maf.gz",
                "luad-v3.1-acetylome-ratio-norm-NArm.gct.gz",
                "luad-v3.1-cnv-gene-LR.gct.gz",
                "luad-v3.1-mirna-mature-tpm-log2.gct.gz",
                "luad-v3.1-phosphoproteome-ratio-norm-NArm.gct.gz",
                "luad-v3.1-proteome-ratio-norm-NArm.gct.gz",
                "luad-v3.1-rnaseq-linc-uq-rpkm-log2-NArm.gct.gz",
                "luad-v3.1-rnaseq-prot-uq-rpkm-log2-NArm.gct.gz",
                "luad-v3.1-sample-annotation.csv.gz"],
            "3.1.1": [
                "LUAD_followup_9_12.xlsx",
                "luad-v3.0-rnaseq-circ-rna_parsed.tsv.gz",
                "luad-v3.0-rnaseq-gene-fusions.csv.gz",
                "luad-v3.0-wxs-somatic.luad.v1.4.20190517.maf.gz",
                "luad-v3.1-acetylome-ratio-norm-NArm.gct.gz",
                "luad-v3.1-cnv-gene-LR.gct.gz",
                "luad-v3.1-mirna-mature-tpm-log2.gct.gz",
                "luad-v3.1-phosphoproteome-ratio-norm-NArm.gct.gz",
                "luad-v3.1-proteome-ratio-norm-NArm.gct.gz",
                "luad-v3.1-rnaseq-linc-uq-rpkm-log2-NArm.gct.gz",
                "luad-v3.1-rnaseq-prot-uq-rpkm-log2-NArm.gct.gz",
                "luad-v3.1-sample-annotation.csv.gz"],
        }

        super().__init__(cancer_type="luad", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            if file_name == "luad-v3.0-rnaseq-gene-fusions.csv.gz":
                df = pd.read_csv(file_path)
                df = df.rename(columns={"Sample.ID": "Patient_ID"})
                df = df.set_index("Patient_ID")

                self._data["gene_fusion"] = df

            elif file_name == "luad-v3.0-wxs-somatic.luad.v1.4.20190517.maf.gz":
                df = pd.read_csv(file_path, sep='\t')
                df = df.rename(columns={"Sample.ID": "Patient_ID"})

                df = df[['Patient_ID','Hugo_Symbol','Variant_Classification','HGVSp_Short']]
                df = df.rename(columns={
                    "Hugo_Symbol":"Gene",
                    "Variant_Classification":"Mutation",
                    "HGVSp_Short":"Location"}) # Rename the columns we want to keep to the appropriate names

                df = df.sort_values(by=["Patient_ID", "Gene"])
                df = df.set_index("Patient_ID")

                self._data["somatic_mutation"] = df

            elif file_name == "luad-v3.1-acetylome-ratio-norm-NArm.gct.gz":
                df = pd.read_csv(file_path, sep="\t", skiprows=2, dtype=object)
                gene_filter = df['geneSymbol'] != 'na' #Drop rows of metadata
                df = df[gene_filter]

                # Prepare some columns we'll need later for the multiindex
                df["variableSites"] = df["variableSites"].str.replace(r"[a-z\s]", "", regex=True) # Get rid of all lowercase delimeters and whitespace in the sites
                df = df.rename(columns={
                    "geneSymbol": "Name",
                    "variableSites": "Site",
                    "sequence": "Peptide", # We take this instead of sequenceVML, to match the other datasets' format
                    "accession_numbers": "Database_ID" # We take all accession numbers they have, instead of the singular accession_number column
                    })


                # Some rows have at least one localized acetylation site, but also have other acetylations that aren't localized. We'll drop those rows, if their localized sites are duplicated in another row, to avoid creating duplicates, because we only preserve information about the localized sites in a given row. However, if the localized sites aren't duplicated in another row, we'll keep the row.
                localization = df["accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA"].str.split("_", expand=True)
                unlocalized_to_drop = localization.index[~localization[3].eq(localization[4]) & df.duplicated(["Name", "Site", "Peptide", "Database_ID"], keep=False)] # Column 3 of the split localization data column is number of acetylation sites detected, and column 4 is number of acetylation sites localized, so if the two values aren't equal, the row has at least one unlocalized site
                df = df.drop(index=unlocalized_to_drop)

                # Give it a multiindex
                df = df.set_index(["Name", "Site", "Peptide", "Database_ID"])                

                cols_to_drop = [
                    "id",
                    "id.description",
                    "numColumnsVMsiteObserved",
                    "bestScore",
                    "bestDeltaForwardReverseScore",
                    "Best_scoreVML",
                    "sequenceVML",
                    "accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA",
                    "protein_mw",
                    "species",
                    "speciesMulti",
                    "orfCategory",
                    "accession_number",
                    "protein_group_num",
                    "entry_name",
                    "GeneSymbol",
                    ]
                df = df.drop(columns=cols_to_drop)

                df = df.apply(pd.to_numeric)
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name="Patient_ID"
                
                self._data["acetylproteomics"] = df

            elif file_name == "luad-v3.1-mirna-mature-tpm-log2.gct.gz":
                df = pd.read_csv(file_path, sep="\t", skiprows=2, dtype=object)

                # Filter out metadata rows
                gene_filter = df["Name"] != 'na' 
                df = df[gene_filter]

                df = df.drop(columns=['id', 'Alias', 'Name', 'Derives_from', 'Quantified.in.Percent.Samples'])
                df = df.set_index("ID")
                df = df.apply(pd.to_numeric)
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()

                self._data["miRNA"] = df

            elif file_name == "luad-v3.1-rnaseq-linc-uq-rpkm-log2-NArm.gct.gz":
                df = pd.read_csv(file_path, sep="\t", skiprows=2, dtype=object)

                # Filter out metadata rows
                gene_filter = df["geneSymbol"] != 'na' 
                df = df[gene_filter]

                # Filter out just lincRNA. Current it has a bunch of other RNA types too. We'll worry about them later.
                lincRNA_filter = df["gene_type"] == "lincRNA"
                df = df[lincRNA_filter]

                # Drop unneeded columns, set index, sort, transpose
                df = df.drop(columns=["id", "gene_id", "gene_type", "length"])
                df = df.rename(columns={"geneSymbol": "Name"})
                df = df.set_index("Name")
                df = df.apply(pd.to_numeric)
                df = df.sort_index()
                df = df.transpose()
                df.index.name = "Patient_ID"

                self._data["lincRNA"] = df

            elif file_name == "luad-v2.0-cnv-gene-LR.gct.gz" or file_name == "luad-v3.1-cnv-gene-LR.gct.gz":
                df = pd.read_csv(file_path, sep="\t", skiprows=2, dtype=object)

                # Filter out metadata rows
                if self._version == "2.0":
                    gene_filter = df['Description'] != 'na' 

                elif self._version in ["3.1", "3.1.1"]:
                    gene_filter = df['geneSymbol'] != 'na' 

                df = df[gene_filter]

                # Drop non-quantitative columns
                if self._version == "2.0":
                    cols_to_drop = ["GeneID","Description"]
                    
                elif self._version in ["3.1", "3.1.1"]:
                    cols_to_drop = ["id", "gene_id"]

                df = df.drop(columns=cols_to_drop)

                # Set index
                if self._version == "2.0":
                    df = df.set_index("id")

                elif self._version in ["3.1", "3.1.1"]:
                    df = df.set_index("geneSymbol")

                df = df.apply(pd.to_numeric)
                df = df.sort_index()
                df = df.transpose()

                df = df.sort_index()
                df.index.name="Patient_ID"

                self._data["CNV"] = df

            elif file_name == "luad-v2.0-phosphoproteome-ratio-norm-NArm.gct.gz" or file_name == "luad-v3.1-phosphoproteome-ratio-norm-NArm.gct.gz":
                df = pd.read_csv(file_path, sep="\t", skiprows=2, dtype=object)

                # Drop rows of metadata
                gene_filter = df['geneSymbol'] != 'na' 
                df = df[gene_filter]

                # Prepare some columns we'll need later for the multiindex
                df["variableSites"] = df["variableSites"].str.replace(r"[a-z\s]", "", regex=True) # Get rid of all lowercase delimeters and whitespace in the sites
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

                # Drop non-quantitative columns
                if self._version == "2.0":
                    cols_to_drop = ['id', 'id.1', 'id.description', 'numColumnsVMsiteObserved', 'bestScore', 'bestDeltaForwardReverseScore',
                        'Best_scoreVML', 'Best_numActualVMSites_sty', 'Best_numLocalizedVMsites_sty', 'sequenceVML',
                        'accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA', 'protein_mw', 'species',
                        'speciesMulti', 'orfCategory', 'accession_number', 'protein_group_num', 'entry_name', 'GeneSymbol']
                elif self._version in ["3.1", "3.1.1"]:
                    cols_to_drop = ["id", "id.description", "numColumnsVMsiteObserved", "bestScore", "bestDeltaForwardReverseScore",
                        "Best_scoreVML", "Best_numActualVMSites_sty", "Best_numLocalizedVMsites_sty", "sequenceVML",
                        "accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA", "protein_mw", "species",
                        "speciesMulti", "orfCategory", "accession_number", "protein_group_num", "entry_name", "GeneSymbol"]
                
                df = df.drop(columns=cols_to_drop)

                # Format table
                df = df.apply(pd.to_numeric)
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name="Patient_ID"

                self._data["phosphoproteomics"] = df

            elif file_name == "luad-v2.0-proteome-ratio-norm-NArm.gct.gz" or file_name == "luad-v3.1-proteome-ratio-norm-NArm.gct.gz":
                df = pd.read_csv(file_path, skiprows=2, sep='\t', dtype=object)

                # Filter out rows of metadata
                gene_filter = df['geneSymbol'] != 'na' 
                df = df[gene_filter]

                # Set multiindex
                df = df.rename(columns={"GeneSymbol": "Name", 'accession_numbers': "Database_ID"})
                df = df.set_index(["Name", "Database_ID"])

                # Drop non-quantitative columns
                if self._version == "2.0":
                    cols_to_drop = ['id', 'id.1', 'id.description', 'geneSymbol', 'numColumnsProteinObserved', 'numSpectraProteinObserved',
                        'protein_mw', 'percentCoverage', 'numPepsUnique', 'scoreUnique', 'species', 'orfCategory', 'accession_number',
                        'subgroupNum', 'entry_name']
                elif self._version in ["3.1", "3.1.1"]:
                    cols_to_drop = ["id", "id.description", "geneSymbol", "numColumnsProteinObserved", "numSpectraProteinObserved",
                        "protein_mw", "percentCoverage", "numPepsUnique", "scoreUnique", "species", "orfCategory", "accession_number",
                        "subgroupNum", "entry_name"]

                df = df.drop(columns=cols_to_drop)

                # Format table
                df = df.apply(pd.to_numeric)
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"

                self._data["proteomics"] = df

            elif file_name == "luad-v2.0-rnaseq-prot-uq-rpkm-log2-NArm-row-norm.gct.gz" or file_name == "luad-v3.1-rnaseq-prot-uq-rpkm-log2-NArm.gct.gz":
                df = pd.read_csv(file_path, sep="\t", skiprows=2, dtype=object)

                # Filter out metadata rows
                gene_filter = df['geneSymbol'] != 'na'
                df = df[gene_filter]

                df = df.set_index('geneSymbol')
                cols_to_drop = ['id', 'gene_id', 'gene_type', 'length']
                df = df.drop(columns = cols_to_drop)

                df = df.apply(pd.to_numeric)
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                df = df.sort_index()

                self._data["transcriptomics"] = df

            elif file_name == "luad-v3.0-rnaseq-circ-rna_parsed.tsv.gz" and self._version == "3.1.1":
                df = pd.read_csv(file_path, sep='\t', dtype={"spanning.reads": "int16"}, engine="c")
                df = df.reset_index() # More memory efficient to do this here, rather than save with a range index when we parse the file it originally
                df = df.pivot(index="Sample.ID", columns="geneID", values="spanning.reads")
                df.index.name = "Patient_ID"
                df = df.sort_index()
                self._data['circular_RNA'] = df

            elif file_name in ["luad-v2.0-rnaseq-circ-rna.csv.gz", "luad-v3.0-rnaseq-circ-rna.csv.gz"] and self._version in ["2.0", "3.1"]:
                df = pd.read_csv(file_path, sep=",")

                junct_3_split = df['junction.3'].str.split(':', n=2, expand=True)
                chrm = junct_3_split[0] # Get the chromosome
                three_prime = junct_3_split[1] # Get the nucleotide coordinate of the last base of the acceptor

                junct_5_split = df['junction.5'].str.split(':', n=2, expand=True)
                five_prime = junct_5_split[1] # Get the nucleotide coordinates of the first base of the donor

                # Now we need the gene name
                diff = df['gene.5'] != df['gene.3'] # Create a boolean filter where genes are different
                temp = df['gene.5'].where(diff, other="") # Replace the ones that are the same with an empty string
                gene_name = temp + '_' + df["gene.3"] # Concatentate the temp column(which only has the genes from gene.5 that are different) to gene.3

                # Put all those pieces of information together
                df = df.assign(geneID=chrm + '_' + five_prime + '_' + three_prime + '_' + gene_name)

                # Slice out the columns we want
                df = df[['geneID', 'spanning.reads', 'Sample.ID']]

                #There are about 3,000 duplicates in the file. Duplicate meaning that they have identical Sample IDs and identical geneID, but different spanning reads.
                # Marcin Cieslik said to drop the one with the lowest spanning read.
                df = df.sort_values(by='spanning.reads', ascending=False).drop_duplicates(['Sample.ID','geneID']).sort_index()

                df = df.pivot(index="Sample.ID", columns="geneID")['spanning.reads']
                df.index.name = "Patient_ID"
                df = df.sort_index()
                self._data['circular_RNA'] = df

            elif file_name == "luad-v2.0-sample-annotation.csv.gz" or file_name == "luad-v3.1-sample-annotation.csv.gz":
                df = pd.read_csv(file_path, sep=",")
                
                filter = df['QC.status'] == "QC.pass" # There are some samples that are internal references. IRs are used for scaling purposes, and don't belong to a single patient, so we want to drop them.
                df = df[filter]

                df = df.drop(columns="Participant") # Get rid of the "Participant" column becuase the same information is stored in Sample.ID  which is formatted the way we want.
                df = df.set_index("Sample.ID")
                df.index.name="Patient_ID"
                df = df.rename(columns={"Type":"Sample_Tumor_Normal"})
                df["Sample_Tumor_Normal"] = df["Sample_Tumor_Normal"].replace(to_replace="NAT", value="Normal")

                # Split the metadata into multiple dataframes

                # Make experimental_design dataframe
                if self._version == "2.0":
                    experimental_design_cols = ['Experiment', 'Channel', 'QC.status']
                elif self._version in ["3.1", "3.1.1"]:
                    experimental_design_cols = ["Experiment", "Channel", "Aliquot", "QC.status"]

                experimental_design_df = df[experimental_design_cols]
                experimental_design_df.insert(0, "Sample_Tumor_Normal", df["Sample_Tumor_Normal"].copy()) # This is useful in both tables
                df = df.drop(columns=experimental_design_cols)

                # Make derived_molecular dataframe
                if self._version == "2.0":
                    derived_molecular_cols = ['TP53.mutation', 'KRAS.mutation', 'STK11.mutation', 'EGFR.mutation', 'KEAP1.mutation', 'RB1.mutation',
                        'IL21R.mutation', 'EGFL6.mutation', 'LMO2.mutation', 'C10orf62.mutation', 'DKK3.mutation', 'BIRC6.mutation', 'TP53.mutation.status',
                        'KRAS.mutation.status', 'STK11.mutation.status', 'EGFR.mutation.status', 'KEAP1.mutation.status', 'RB1.mutation.status', 'IL21R.mutation.status',
                        'EGFL6.mutation.status', 'LMO2.mutation.status', 'C10orf62.mutation.status', 'DKK3.mutation.status', 'BIRC6.mutation.status',
                        'Mutation.Signature.Activity.W1.COSMIC5', 'Mutation.Signature.Activity.W2.COSMIC4', 'Mutation.Signature.Activity.W3.COSMIC2', 'fusion.EML4-ALK']

                elif self._version in ["3.1", "3.1.1"]:
                    derived_molecular_cols = ["Smoking.Score.WGS", "Smoking.Signature.Fraction.WGS", "Dominant.Signature.WGS.notSmoking.50perc",
                        "Dominant.Signature.Fraction.WGS.notSmoking", "DNP.GG.to.TT.or.CC.to.AA.Count.WGS", "NMF.consensus", "NMF.cluster.membership",
                        "mRNA.Expression.Subtype.TCGA", "mRNA.stemness.index", "CIMP.status", "Tumor.Purity.byESTIMATE.RNAseq", "TSNet Purity", "ESTIMATEScore",
                        "ESTIMATE ImmuneScore", "ESTIMATE StromalScore", "TP53.mutation", "KRAS.mutation", "STK11.mutation", "EGFR.mutation", "KEAP1.mutation",
                        "RB1.mutation", "IL21R.mutation", "EGFL6.mutation", "LMO2.mutation", "C10orf62.mutation", "DKK3.mutation", "BIRC6.mutation",
                        "BRAF.mutation", "ARAF.mutation", "ERBB2.mutation", "TP53.mutation.status", "KRAS.mutation.status", "STK11.mutation.status",
                        "EGFR.mutation.status", "KEAP1.mutation.status", "RB1.mutation.status", "IL21R.mutation.status", "EGFL6.mutation.status", "LMO2.mutation.status",
                        "C10orf62.mutation.status", "DKK3.mutation.status", "BIRC6.mutation.status", "BRAF.mutation.status", "ARAF.mutation.status",
                        "ERBB2.mutation.status", "Total.Mutation.Count.WGS", "Mutation.Count.ExcludingINDELs.WGS", "Total.DNP.Count.WGS", "Number.somatic.mutations",
                        "Mutation.Signature.Activity.W1.COSMIC5", "Mutation.Signature.Activity.W2.COSMIC4", "Mutation.Signature.Activity.W3.COSMIC2", "ALK.fusion",
                        "ROS1.fusion", "RET.fusion", "Putative.driver.mutation"]

                derived_molecular_df = df[derived_molecular_cols]
                df = df.drop(columns = derived_molecular_cols)

                self._data["clinical"] = df
                self._data['experimental_design'] = experimental_design_df
                self._data['derived_molecular'] = derived_molecular_df

            elif file_name == 'LUAD_followup_9_12.xlsx' and self._version in ["3.1", "3.1.1"]:
                df = pd.read_excel(file_path)

                # Replace redundant values for "not reported" with NaN
                nan_equivalents = ['Not Reported/ Unknown', 'Reported/ Unknown', 'Not Applicable',
                    'na', 'unknown', 'Not Performed', 'Unknown tumor status', 'Unknown',
                    'Unknown Tumor Status', 'Not specified']

                df = df.replace(nan_equivalents, np.nan)

                # Replace redundant values for cause of death
                disease_prog_equivalents = ['Progression of disease', 'Progression of disease ', 'Tumor', 'Disease progression',
                    'Progressive Disease', 'disease progression', 'disease progression ', 'main disease ']

                df['Cause of Death'] = df['Cause of Death'].replace(disease_prog_equivalents, 'Disease progression')

                # Rename, set, and sort by index
                df = df.rename(columns={"Case ID": "Patient_ID"})
                df = df.set_index("Patient_ID")
                df = df.sort_index()

                self._data["followup"] = df

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Get a union of all dataframes' indices, with duplicates removed
        # Exclude the followup dataframe because it has samples from a different cohort that aren't included anywhere else in the dataset
        master_index = unionize_indices(self._data, exclude="followup")

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
        clinical = self._data["clinical"]
        clinical = clinical.reindex(master_index)
        self._data['clinical'] = clinical

        # Drop samples C3N.00545 and C3N.00545.N from the dataset. They were excluded due to poor sample quality (see data freeze README; excluded in data freeze 3.0)
        cases_to_drop = ["C3N.00545", "C3N.00545.N"]
        for name in self._data.keys(): # Loop over the keys so we can alter the values without any issues
            df = self._data[name]
            df = df.drop(index=cases_to_drop, errors="ignore")
            self._data[name] = df

        # Replace periods with hyphens in all Patient_IDs
        for name in self._data.keys(): # Loop over just the keys to avoid any issues that would come if we looped over the values while editing them
            df = self._data[name]

            # Make the index a column so we can edit it
            df.index.name = "Patient_ID"
            df = df.reset_index()
            
            # Replace all '.' with '-'
            df["Patient_ID"] = df["Patient_ID"].str.replace(r"\.", "-", regex=True) 
            df["Patient_ID"] = df["Patient_ID"].str.replace(r"-N$", ".N", regex=True) # If there's a "-N" at the end, it's part of the normal identifier, which we want to actually be ".N"
            
            # Set the index back to Patient_ID
            df = df.set_index("Patient_ID")

            self._data[name] = df

        # Call function from dataframe_tools.py to sort all tables first by sample status, and then by the index
        self._data = sort_all_rows(self._data)

        # Call function from dataframe_tools.py to standardize the names of the index and column axes
        self._data = standardize_axes_names(self._data)

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
