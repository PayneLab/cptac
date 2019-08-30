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
from .exceptions import NoInternetError, FailedReindexWarning


class Hnscc(DataSet):

    def __init__(self, version="latest"):
        """Load all of the hnscc dataframes as values in the self._data dict variable, with names as keys, and format them properly."""

        # Call the parent DataSet __init__ function, which initializes self._data and other variables we need
        super().__init__("hnscc")

        # Update the index, if possible. If there's no internet, that's fine.
        try:
            update_index(self._cancer_type)
        except NoInternetError:
            pass

        # Validate the version
        self._version = validate_version(version, self._cancer_type, use_context="init")

        # Get the paths to all the data files
        data_files = [
        "HNSCC.strelka.sorted.filtered.annovar.hg19_multianno_filtered.maf.txt.gz",
        "Proteomics_DIA_Gene_level_Normal.cct.gz",
        "Proteomics_DIA_Gene_level_Tumor.cct.gz",
        "RNAseq_RSEM_UQ_log2.cct.gz",
        "RNAseq_circ_RSEM_UQ_log2.cct.gz",
        "SCNA_gene_level.cct.gz",
        "clinic.tsi.gz"
        ]
        data_files_paths = get_version_files_paths(self._cancer_type, self._version, data_files)

        # Load the data into dataframes in the self._data dict
        loading_msg = "Loading dataframes"
        for file_path in data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            df_name = file_name.split(".")[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            if file_name == "SCNA_gene_level.cct.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.columns.name=None
                df.index.name = "Patient_ID"
                self._data["CNV"] = df

            elif file_name == "RNAseq_RSEM_UQ_log2.cct.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.columns.name=None
                df.index = df.index.str.replace(r'.', '-', 1)
                df.index.name = "Patient_ID"
                self._data["transcriptomics"] = df

            elif file_name == "RNAseq_circ_RSEM_UQ_log2.cct.gz":
                df = pd.read_csv(file_path, sep='\t')
                df = df.sort_index()
                df = df.transpose()
                df = df.sort_index()
                df.columns.name=None
                df.index = df.index.str.replace(r'.', '-', 1) #We want all the patientIDs to have the the format C3L-00977, and these have the form C3L.00977.N, so we need to replace the first "." with a "-"
                df.index.name = "Patient_ID"
                self._data["circular_RNA"] = df

            elif file_name == "HNSCC.strelka.sorted.filtered.annovar.hg19_multianno_filtered.maf.txt.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.rename(columns={"Tumor_Sample_Barcode":"Patient_ID","Hugo_Symbol_Annovar":"Gene","Variant_Classification_Annovar":"Mutation"}) #Rename the columns we want to keep to the appropriate names
                df['Location'] = df['Annovar_Info_protein'].str.extract(r'([^:]+$)') #The location that we care about is stored after the last colon
                keep = ['Gene', 'Mutation', 'Location', 'Patient_ID']
                df = df.drop(df.columns.difference(keep), 1)
                df = df.set_index("Patient_ID")
                df = df.sort_index()
                df = df.sort_index(axis='columns')
                df.columns.name=None
                self._data["somatic_mutation"] = df

            elif file_name == "clinic.tsi.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.set_index('CASE_ID')
                df.columns.name=None
                df.index.name="Patient_ID"
                #Split the clinicl data in to clincial data and derived molecular data
                derived_molecular_cols = ['P53GENE_ANALYSIS', 'EGFR_AMP_STATUS']
                derived_molecular_df = df[derived_molecular_cols]
                derived_molecular_df = derived_molecular_df.sort_index(axis='columns')
                derived_molecular_df = derived_molecular_df.sort_index()
                df = df.drop(columns=derived_molecular_cols)
                df = df.sort_index()
                df = df.sort_index(axis='columns')
                self._data["clinical"] = df
                self._data["derived_molecular"] = derived_molecular_df

            #the only files left are the proteomics files
            elif file_name == "Proteomics_DIA_Gene_level_Normal.cct.gz" or file_name == "Proteomics_DIA_Gene_level_Tumor.cct.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.transpose()
                df.columns.name=None
                df.index.name = "Patient_ID"

                #Once the files are formatted correctly load them into self._data
                if file_name == "Proteomics_DIA_Gene_level_Normal.cct.gz":
                    self._data["proteomics_normal"] = df
                elif file_name == "Proteomics_DIA_Gene_level_Tumor.cct.gz":
                    self._data["proteomics_tumor"] = df


        #get the proteomics data
        df_normal = self._data.get("proteomics_normal")
        df_tumor = self._data.get("proteomics_tumor")

        df_normal.index = "N" + df_normal.index #concatenate an ".N" onto the end of the normal data so we can identify it as normal after it's appended to tumor
        prot_combined = df_tumor.append(df_normal) #append the normal data onto the end of the tumor data
        prot_combined = prot_combined.sort_index(axis='columns') # Put all the columns in alphabetical order
        prot_combined = prot_combined.sort_index()
        self._data["proteomics"] = prot_combined
        del self._data["proteomics_normal"]
        del self._data["proteomics_tumor"]


        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = "Formatting dataframes..."
        print(formatting_msg, end='\r')



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
        dfs_to_delete = [] #
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

        # Drop name of column axis for all dataframes
        for name in self._data.keys():
            df = self._data[name]
            df.columns.name = None
            self._data[name] = df



        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
