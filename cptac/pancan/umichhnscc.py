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

from cptac.dataset import Dataset
from cptac.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError
from cptac.utils import get_boxnote_text


class UmichHnscc(Dataset):

    def __init__(self, no_internet, version):
        """Load all of the umichhnscc dataframes as values in the self._data dict variable, with names as keys, and format them properly.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest building. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        valid_versions = ["1.0"]

        data_files = {
            "1.0": ["Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv",
                    "Report_abundance_groupby=multi-site_protNorm=MD_gu=2.tsv",
                    "README_v3.boxnote", # proteomics 
                    "README.boxnote" # phosphoproteomics
                #"S039_BCprospective_observed_0920.tsv.gz",
                #"S039_BCprospective_imputed_0920.tsv.gz"
            ]
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="umichhnscc", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below
            
            
            if file_name == "Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv":
                df = pd.read_csv(file_path, sep = "\t") 
                df['Database_ID'] = df.Index.apply(lambda x: x.split('|')[0]) # get protein identifier 
                df['Name'] = df.Index.apply(lambda x: x.split('|')[6]) # get protein name 
                df = df.set_index(['Name', 'Database_ID']) # set multiindex
                df = df.drop(columns = ['Index', 'MaxPepProb', 'NumberPSM', 'Gene']) # drop unnecessary  columns
                df = df.transpose()
                ref_intensities = df.loc["ReferenceIntensity"] # get reference intensities to use to calculate ratios 
                df = df.subtract(ref_intensities, axis="columns") # subtract reference intensities from all the values 
                df = df.iloc[1:,:] # drop ReferenceIntensity row 
                df.index.name = 'Patient_ID'
                self._data["proteomics"] = df              
                
            elif file_name == "Report_abundance_groupby=multi-site_protNorm=MD_gu=2.tsv":
                df = pd.read_csv(file_path, sep = "\t") 
               # Parse a few columns out of the "Index" column that we'll need for our multiindex
                df[['Database_ID','Transcript_ID',"Gene_ID","Havana_gene","Havana_transcript","Transcript","Name", "Site"]] = df.Index.str.split("\\|",expand=True)
                df[['num1','start',"end","detected_phos","localized_phos","Site"]] = df.Site.str.split("_",expand=True) 

                # Some rows have at least one localized phosphorylation site, but also have other
                # phosphorylations that aren't localized. We'll drop those rows, if their localized 
                # sites are duplicated in another row, to avoid creating duplicates, because we only 
                # preserve information about the localized sites in a given row. However, if the localized 
                # sites aren't duplicated in another row, we'll keep the row.
                unlocalized_to_drop = df.index[~df["detected_phos"].eq(df["localized_phos"]) & \
                                               df.duplicated(["Name", "Site", "Peptide", "Database_ID"], keep=False)]
                # dectected_phos of the split "Index" column is number of phosphorylations detected, and 
                # localized_phos is number of phosphorylations localized, so if the two values aren't equal, the 
                # row has at least one unlocalized site
                df = df.drop(index=unlocalized_to_drop)
                df = df[df['Site'].notna()] # only keep columns with phospho site 
                df = df.set_index(['Name', 'Site', 'Peptide', 'Database_ID']) # create a multiindex, in this order.
                #drop columns not needed in df 
                df.drop(['Gene', "Index", "num1", "start", "end", "detected_phos", "localized_phos", "Havana_gene", 
                         "Havana_transcript", "MaxPepProb", "Gene_ID", "Transcript_ID", "Transcript"], axis=1, inplace=True)
                df = df.T #transpose df 
                ref_intensities = df.loc["ReferenceIntensity"]# Get reference intensities to use to calculate ratios 
                df = df.subtract(ref_intensities, axis="columns") # Subtract ref intensities from all the values, to get ratios
                df = df.iloc[1:,:] # drop ReferenceIntensity row 
                self._data["phosphoproteomics"] = df
                
            elif file_name == "README_v3.boxnote":
                self._readme_files["readme_proteomics"] = get_boxnote_text(file_path)
                
            elif file_name == "README.boxnote":
                self._readme_files["readme_phosphoproteomics"] = get_boxnote_text(file_path)
            
            '''
            if file_name == "S039_BCprospective_observed_0920.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.transpose()
                df.index.name = 'Patient_ID'
                df.columns.name = 'Name'
                df = average_replicates(df)
                df = df.sort_values(by=["Patient_ID"])
                self._data["proteomics"] = df  
                
            if file_name == "S039_BCprospective_imputed_0920.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.transpose()
                df.index.name = 'Patient_ID'
                df.columns.name = 'Name'
                df = average_replicates(df)
                df = df.sort_values(by=["Patient_ID"])
                self._data["proteomics_imputed"] = df'''
                
        
        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = f"Formatting {self.get_cancer_type()} dataframes..."
        print(formatting_msg, end='\r')               
        
        # There were 4 labels with "-duplicate" appended in proteomics and phosphoproteomics files.
        # I ran a pearson correlation to check how well the values from each duplicate correlated to 
        # the other duplicates for the same case ID. Three of the duplicates correlated well with their 
        # respective case IDs. C3L-02617-N-duplicate2 did not correlate well with the other C3L-02617 duplicates, 
        # so we dropped it and averaged the other two. I also created a scatterplot to compare each duplicate to 
        # the first occurence of its case ID. The linear scatterplots indicated similarity between the aliquots. 
        # We averaged the duplicates that correlated well together and were the same tissue type.        
        # A file containing the correlations can be downloaded at: 
        # https://byu.box.com/shared/static/jzsq69bd079oq0zbicw4w616hyicd5ev.xlsx
        
        drop_cols = ['128C', 'QC2', 'QC3', 'QC4', '129N', 'LungTumor1', 'Pooled-sample14',
                   'LungTumor2', 'QC6', 'LungTumor3', 'Pooled-sample17', 'QC7',
                   'Pooled-sample19', 'QC9', 'RefInt_pool01', 'RefInt_pool02',
                   'RefInt_pool03', 'RefInt_pool04', 'RefInt_pool05', 'RefInt_pool06',
                   'RefInt_pool07', 'RefInt_pool08', 'RefInt_pool09', 'RefInt_pool10',
                   'RefInt_pool11', 'RefInt_pool12', 'RefInt_pool13', 'RefInt_pool14',
                   'RefInt_pool15', 'RefInt_pool16', 'RefInt_pool17', 'RefInt_pool18',
                   'RefInt_pool19', 'RefInt_pool20']    
        
        # Proteomics
        prot = self._data["proteomics"]        
        prot = prot.drop(drop_cols, axis = 'index') # drop quality control and ref intensity cols        
        prot = prot.drop(['C3L-02617-N-duplicate2'], axis = 'index') # drop duplicate that did not correlate well  
        # These IDs had a high correlation with their respective duplicates, so we average them
        # duplicates: 'C3L-02617-T-duplicate', 'C3L-00994-N-duplicate', 'C3L-02617-N-duplicate'
        prot = average_replicates(prot, ['C3L-02617-T','C3L-02617-N','C3L-00994-N'], normal_identifier = '-N') 
        prot.index = prot.index.str.replace('-T$','', regex = True)
        prot.index = prot.index.str.replace('-N$','.N', regex = True)
        prot.index = prot.index.str.replace('-C$','.C', regex = True) # 6 cored normal samples in Hnscc
        self._data["proteomics"] = prot
        
        # Phosphoproteomics
        phos = self._data["phosphoproteomics"]
        phos = phos.drop(drop_cols, axis = 'index') # drop quality control and ref intensity cols        
        phos = phos.drop(['C3L-02617-N-duplicate2'], axis = 'index') # drop duplicate that did not correlate well
        # average IDs that correlated well to their respective duplicates
        phos = average_replicates(phos, ['C3L-02617-T','C3L-02617-N','C3L-00994-N'], normal_identifier = '-N') 
        phos.index = phos.index.str.replace('-T$','', regex = True)
        phos.index = phos.index.str.replace('-N$','.N', regex = True)
        phos.index = phos.index.str.replace('-C$','.C', regex = True) # 6 cored normal samples in Hnscc
        self._data["phosphoproteomics"] = phos

        
        # Sort rows (tumor first then normal) and columns by first level (protein/gene name)
        self._data = sort_all_rows_pancan(self._data) 
        

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
