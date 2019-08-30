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
from .file_download import update_index
from .file_tools import validate_version, get_version_files_paths
from .dataframe_tools import *
from .exceptions import NoInternetError

class Brca(DataSet):

    def __init__(self, version="latest"):
        """Load all of the brca dataframes as values in the self._data dict variable, with names as keys, and format them properly."""

        # Call the parent DataSet __init__ function, which initializes self._data and other variables we need
        super().__init__("brca")

        # Update the index, if possible. If there's no internet, that's fine.
        try:
            update_index(self._cancer_type)
        except NoInternetError:
            pass

        # Validate the version
        self._version = validate_version(version, self._cancer_type, use_context="init")

        # Get the paths to all the data files
        data_files = [
            "prosp-brca-v3.1-acetylome-ratio-norm-NArm.gct.gz",
            "prosp-brca-v3.1-gene-level-cnv-gistic2-all_data_by_genes.gct.gz",
            "prosp-brca-v3.1-phosphoproteome-ratio-norm-NArm.gct.gz",
            "prosp-brca-v3.1-proteome-ratio-norm-NArm.gct.gz",
            "prosp-brca-v3.1-rnaseq-fpkm-log2-row-norm-2comp.gct.gz",
            "prosp-brca-v3.1-sample-annotation.csv.gz"] 
        data_files_paths = get_version_files_paths(self._cancer_type, self._version, data_files)

        # Load the data into dataframes in the self._data dict
        loading_msg = "Loading dataframes"
        for file_path in data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file

            if file_name == "prosp-brca-v3.1-acetylome-ratio-norm-NArm.gct.gz":
                df = pd.read_csv(file_path, sep='\t', skiprows=2, dtype=object) # First two rows of file aren't part of the dataframe. Also, due to extra metadata rows we're going to remove, all cols have mixed types, so we pass dtype=object for now.
                df = df[df["GeneSymbol"] != "na"] # There are several metadata rows at the beginning of the dataframe, which duplicate the clinical and derived_molecular dataframes. They all don't have a value for GeneSymbol, so we'll use that to filter them out.

                # Prepare some columns we'll need later for the multiindex
                df["variableSites"] = df["variableSites"].str.replace(r"[a-z\s]", "") # Get rid of all lowercase delimeters and whitespace in the sites
                df = df.rename(columns={
                    "GeneSymbol": "Name",
                    "variableSites": "Site",
                    "sequence": "Peptide", # We take this instead of sequenceVML, to match the other datasets' format
                    "accession_numbers": "Database_ID" # We take all accession numbers they have, instead of the singular accession_number column
                    })

                # Some rows have at least one localized acetylation site, but also have other acetylations that aren't localized. We'll drop those rows, if their localized sites are duplicated in another row, to avoid creating duplicates, because we only preserve information about the localized sites in a given row. However, if the localized sites aren't duplicated in another row, we'll keep the row.
                split_ids = df["id"].str.split('_', expand=True)
                unlocalized_to_drop = df.index[~split_ids[3].eq(split_ids[4]) & df.duplicated(["Name", "Site", "Peptide", "Database_ID"], keep=False)] # Column 3 of the split "id" column is number of phosphorylations detected, and column 4 is number of phosphorylations localized, so if the two values aren't equal, the row has at least one unlocalized site
                df = df.drop(index=unlocalized_to_drop)

                # Give it a multiindex
                df = df.set_index(["Name", "Site", "Peptide", "Database_ID"])                

                df = df.drop(columns=["id", "id.description", "geneSymbol", "numColumnsVMsiteObserved", "bestScore", "bestDeltaForwardReverseScore", 
                "Best_scoreVML", "sequenceVML", "accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA",
                "protein_mw", "species", "speciesMulti", "orfCategory", "accession_number", "protein_group_num", "entry_name"]) # We don't need these. The dropped columns include a "geneSymbol" column that is a duplicate of the original GeneSymbol.
                df = df.apply(pd.to_numeric) # Now that we've dropped all the extra metadata columns, convert everything to floats.
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                self._data["acetylproteomics"] = df

            elif file_name == "prosp-brca-v3.1-gene-level-cnv-gistic2-all_data_by_genes.gct.gz":
                df = pd.read_csv(file_path, sep='\t', skiprows=2, index_col=0, dtype=object) # First two rows of file aren't part of the dataframe. Also, due to extra metadata rows we're going to remove, all cols have mixed types, so we pass dtype=object for now.
                df = df[df["geneSymbol"] != "na"] # There are several metadata rows at the beginning of the dataframe, which duplicate the clinical and derived_molecular dataframes. They all don't have a value for geneSymbol, so we'll use that to filter them out.
                df = df.drop(columns="Cytoband")
                df["geneSymbol"] = df["geneSymbol"].str.rsplit('|', n=1, expand=True)[0] # Some of the geneSymbols have the gene IDs appended to them, to get rid of duplicates. We're going to create a multiindex with all the gene names and gene IDs, so we can drop the appended IDs.
                df = df.rename(columns={"geneSymbol": "Name", "Gene.ID": "Database_ID"})
                df = df.set_index(["Name", "Database_ID"])
                df = df.apply(pd.to_numeric) # Now that we've dropped all the extra metadata columns, convert everything to floats.
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                self._data["CNV"] = df

            elif file_name == "prosp-brca-v3.1-phosphoproteome-ratio-norm-NArm.gct.gz":
                df = pd.read_csv(file_path, sep='\t', skiprows=2, dtype=object) # First two rows of file aren't part of the dataframe. Also, due to extra metadata rows we're going to remove, all cols have mixed types, so we pass dtype=object for now.
                df = df[df["GeneSymbol"] != "na"] # There are several metadata rows at the beginning of the dataframe, which duplicate the clinical and derived_molecular dataframes. They all don't have a value for GeneSymbol, so we'll use that to filter them out.

                # Prepare some columns we'll need later for the multiindex
                df["variableSites"] = df["variableSites"].str.replace(r"[a-z\s]", "") # Get rid of all lowercase delimeters and whitespace in the sites
                df = df.rename(columns={
                    "GeneSymbol": "Name",
                    "variableSites": "Site",
                    "sequence": "Peptide", # We take this instead of sequenceVML, to match the other datasets' format
                    "accession_numbers": "Database_ID" # We take all accession numbers they have, instead of the singular accession_number column
                    })

                # Some rows have at least one localized phosphorylation site, but also have other phosphorylations that aren't localized. We'll drop those rows, if their localized sites are duplicated in another row, to avoid creating duplicates, because we only preserve information about the localized sites in a given row. However, if the localized sites aren't duplicated in another row, we'll keep the row.
                split_ids = df["id"].str.split('_', expand=True)
                unlocalized_to_drop = df.index[~split_ids[3].eq(split_ids[4]) & df.duplicated(["Name", "Site", "Peptide", "Database_ID"], keep=False)] # Column 3 of the split "id" column is number of phosphorylations detected, and column 4 is number of phosphorylations localized, so if the two values aren't equal, the row has at least one unlocalized site
                df = df.drop(index=unlocalized_to_drop)

                # Give it a multiindex
                df = df.set_index(["Name", "Site", "Peptide", "Database_ID"])                

                df = df.drop(columns=["id", "id.description", "geneSymbol", "numColumnsVMsiteObserved", "bestScore", "bestDeltaForwardReverseScore",
                "Best_scoreVML", "Best_numActualVMSites_sty", "Best_numLocalizedVMsites_sty", "sequenceVML",
                "accessionNumber_VMsites_numVMsitesPresent_numVMsitesLocalizedBest_earliestVMsiteAA_latestVMsiteAA", "protein_mw", "species",
                "speciesMulti", "orfCategory", "accession_number", "protein_group_num", "entry_name"]) # We don't need these. The dropped columns include a "geneSymbol" column that is a duplicate of the original GeneSymbol.
                df = df.apply(pd.to_numeric) # Now that we've dropped all the extra metadata columns, convert everything to floats.
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                self._data["phosphoproteomics"] = df

            elif file_name == "prosp-brca-v3.1-proteome-ratio-norm-NArm.gct.gz":
                df = pd.read_csv(file_path, sep='\t', skiprows=2, dtype=object) # First two rows of file aren't part of the dataframe. Also, due to extra metadata rows we're going to remove, all cols have mixed types, so we pass dtype=object for now.
                df = df[df["GeneSymbol"] != "na"] # There are several metadata rows at the beginning of the dataframe, which duplicate the clinical and derived_molecular dataframes. They all don't have a value for GeneSymbol, so we'll use that to filter them out.

                df = df.rename(columns={"GeneSymbol": "Name", "accession_numbers": "Database_ID"})
                df = df.set_index(["Name", "Database_ID"])
                df = df.drop(columns=["id", "id.description", "geneSymbol", "numColumnsProteinObserved", "numSpectraProteinObserved",
                "protein_mw", "percentCoverage", "numPepsUnique", "scoreUnique", "species", "orfCategory", "accession_number", 
                "subgroupNum", "entry_name"]) # We don't need these. The dropped columns include a "geneSymbol" column that is a duplicate of GeneSymbol.
                df = df.apply(pd.to_numeric) # Now that we've dropped all the extra metadata columns, convert everything to floats.
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                self._data["proteomics"] = df

            elif file_name == "prosp-brca-v3.1-rnaseq-fpkm-log2-row-norm-2comp.gct.gz":
                df = pd.read_csv(file_path, sep='\t', skiprows=2, index_col=0, dtype=object) # First two rows of file aren't part of the dataframe. Also, due to extra metadata rows we're going to remove, all cols have mixed types, so we pass dtype=object for now.
                df = df[df["geneSymbol"] != "na"] # There are several metadata rows at the beginning of the dataframe, which duplicate the clinical and derived_molecular dataframes. They all don't have a value for GeneSymbol, so we'll use that to filter them out.
                df = df.set_index("geneSymbol")
                df = df.drop(columns="description") # We don't need this.
                df = df.apply(pd.to_numeric) # Now that we've dropped all the extra metadata columns, convert everything to floats.
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.index.name = "Patient_ID"
                self._data["transcriptomics"] = df

            elif file_name == "prosp-brca-v3.1-sample-annotation.csv.gz":
                df = pd.read_csv(file_path, index_col=0)
                df = df.drop(columns="Participant") # This column is just a duplicate of the index
                df = df.rename(columns={"Sample.IDs": "Replicate_Measurement_IDs", "Type": "Sample_Tumor_Normal"})
                df.index.name = "Patient_ID"
                self._data["metadata"] = df

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')

        # Separate the clinical and derived_molecular dataframes
        metadata = self._data["metadata"]
        del self._data["metadata"] # We'll replace it, split into clinical and derived_molecular
        clinical = metadata[[
             "Replicate_Measurement_IDs", "Sample_Tumor_Normal", "Age.in.Month", "Gender", "Race", "Human.Readable.Label", "Experiment", "Channel", "Stage", 
             "PAM50", "NMF.v2.1", "ER", "PR", "ER.IHC.Score", "PR.IHC.Score", "Coring.or.Excision", "Ischemia.Time.in.Minutes", 
             "Ischemia.Decade", "Necrosis", "Tumor.Cellularity", "Total.Cellularity", "In.CR", "QC.status"]]
        self._data["clinical"] = clinical

        derived_molecular = metadata[[
            "HER2.IHC.Score", "HER2.FISH.Status", "HER2.original", "HER2.Amplified", "HER2.refined", "STARD3.ERBB2.GRB7.protein", 
            "HER2.class.Satpathy", "HER2.status.Satpathy", "PAM50.Her2.CNA", "PAM50.Her2.HER2.status", "CDH1.mutation", 
            "GATA3.mutation", "MAP3K1.mutation", "PIK3CA.mutation", "PTEN.mutation", "TP53.mutation", "CDH1.mutation.status", 
            "GATA3.mutation.status", "MAP3K1.mutation.status", "PIK3CA.mutation.status", "PTEN.mutation.status", "TP53.mutation.status", 
            "Number.of.Mutations", "Number.of.Mutated.Genes", "Chromosome.INstability.index.CIN.", "ESTIMATE.TumorPurity", 
            "ESTIMATE.ImmuneScore", "ESTIMATE.StromalScore", "xCell.ImmuneScore", "xCell.StromaScore", "Cibersort.Absolute.score", "Stemness.Score"]]
        self._data["derived_molecular"] = derived_molecular

        # Get a union of all dataframes' indices, with duplicates removed
        master_index = unionize_indices(self._data)

        # Use the master index to reindex the clinical dataframe, so the clinical dataframe has a record of every sample in the dataset. Rows that didn't exist before (such as the rows for normal samples) are filled with NaN.
        clinical = self._data["clinical"]
        clinical = clinical.reindex(master_index)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self._data['clinical'] = clinical

        # Generate a sample ID for each patient ID
        sample_id_dict = generate_sample_id_map(master_index)

        # Give all the dataframes Sample_ID indices
        dfs_to_delete = []
        for name in self._data.keys(): # Only loop over keys, to avoid changing the structure of the object we're looping over
            df = self._data[name]
            df.index.name = "Patient_ID"
            keep_old = (name == "clinical") # Keep the old Patient_ID index as a column in the clinical dataframe, so we have a record of it.
            try:
                df = reindex_dataframe(df, sample_id_dict, "Sample_ID", keep_old)
            except ReindexMapError:
                warnings.warn(f"Error mapping sample ids in {name} dataframe. At least one Patient_ID did not have corresponding Sample_ID mapped in clinical dataframe. {name} dataframe not loaded.", FailedReindexWarning, stacklevel=2) # stacklevel=2 ensures that the warning is registered as originating from the file that called this __init__ function, instead of from here directly, because the former is more useful information.
                dfs_to_delete.append(name)
                continue

            self._data[name] = df

        for name in dfs_to_delete: # Delete any dataframes that had issues reindexing
            del self._data[name]

        # Set name of column axis to "Name" for all dataframes
        for name in self._data.keys():
            df = self._data[name]
            df.columns.name = "Name"
            self._data[name] = df

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
