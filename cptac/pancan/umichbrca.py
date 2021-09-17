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


class UmichBrca(Dataset):

    def __init__(self, no_internet, version):
        """Load all of the umichbrca dataframes as values in the self._data dict variable, with names as keys, and format them properly.

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
                    "prosp-brca-all-samples.txt",
                    "README_v3.boxnote", # proteomics 
                    "README.boxnote" # phosphoproteomics 
            ]
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="umichbrca", version=version, valid_versions=valid_versions, data_files=data_files, no_internet=no_internet)

        # Load the data into dataframes in the self._data dict
        loading_msg = f"Loading {self.get_cancer_type()} v{self.version()}"
        for file_path in self._data_files_paths: # Loops through files variable

            # Print a loading message. We add a dot every time, so the user knows it's not frozen.
            loading_msg = loading_msg + "."
            print(loading_msg, end='\r')

            path_elements = file_path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file. We'll use this to identify files for parsing in the if/elif statements below
                
            #Proteomics
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
                
                # drop ending of CPT retrospective samples to match cptac 
                df = df.rename(index={'CPT0008140004':'CPT000814', 'CPT0018460005': 'CPT001846', 
                                      '604':'CPT000814'}) # 604 mapped to CPT000814 in the proteomics file from the PDC pipeline 
    
                drop_cols = ['RetroIR', 'RetroIR.1',
                   'RefInt_Pool01', 'RefInt_Pool02', 'RefInt_Pool03', 'RefInt_Pool04',
                   'RefInt_Pool05', 'RefInt_Pool06', 'RefInt_Pool07', 'RefInt_Pool08',
                   'RefInt_Pool09', 'RefInt_Pool10', 'RefInt_Pool11', 'RefInt_Pool12',
                   'RefInt_Pool13', 'RefInt_Pool14', 'RefInt_Pool15', 'RefInt_Pool16',
                   'RefInt_Pool17']                
                df = df.drop(drop_cols, axis = 'index') # drop quality control and ref intensity cols
                self._data["proteomics"] = df
                
            #Phosphoproteomics    
            elif file_name == "Report_abundance_groupby=multi-site_protNorm=MD_gu=2.tsv":
                df = pd.read_csv(file_path, sep = "\t") 
                 # Parse a few columns out of the "Index" column that we'll need for our multiindex
                df[['Database_ID','Transcript_ID',"Gene_ID","Havana_gene","Havana_transcript","Transcript","Name","Site"]] = df.Index.str.split("\\|",expand=True)
                df[['num1','start',"end","detected_phos","localized_phos","Site"]] = df.Site.str.split("_",expand=True) 

                # Some rows have at least one localized phosphorylation site, but also have other 
                # phosphorylations that aren't localized. We'll drop those rows, if their localized 
                # sites are duplicated in another row, to avoid creating duplicates, because we only 
                # preserve information about the localized sites in a given row. However, if the localized 
                # sites aren't duplicated in another row, we'll keep the row.
                unlocalized_to_drop = df.index[~df["detected_phos"].eq(df["localized_phos"]) & df.duplicated(["Name", "Site", "Peptide", "Database_ID"], keep=False)]# dectected_phos of the split "Index" column is number of phosphorylations detected, and localized_phos is number of phosphorylations localized, so if the two values aren't equal, the row has at least one unlocalized site
                df = df.drop(index=unlocalized_to_drop)

                df = df[df['Site'].notna()] # only keep columns with phospho site 
                df = df.set_index(['Name', 'Site', 'Peptide', 'Database_ID']) # Create a multiindex in this order.
                #drop columns not needed in df 
                df.drop(['Gene', "Index", "num1", "start", "end", "detected_phos", "localized_phos", "Havana_gene", 
                         "Havana_transcript", "MaxPepProb", "Gene_ID", "Transcript_ID", "Transcript"], axis=1, inplace=True)
                df = df.T #transpose df 
                ref_intensities = df.loc["ReferenceIntensity"]# Get reference intensities to use to calculate ratios 
                df = df.subtract(ref_intensities, axis="columns") # Subtract ref intensities from all the values, to get ratios
                df = df.iloc[1:,:] # drop ReferenceIntensity row 
                # drop ending of CPT retrospective samples to match cptac
                df = df.rename(index={'CPT0008140004':'CPT000814', 'CPT0018460005': 'CPT001846', 
                                      '604':'CPT000814'}) # 604 mapped to CPT000814 in pdc index
                
                drop_cols = ['RetroIR','RetroIR.1','RefInt_Pool01-1','RefInt_Pool02-1',
                             'RefInt_Pool03-1','RefInt_Pool04-1','RefInt_Pool05-1','RefInt_Pool06-1',
                             'RefInt_Pool07-1','RefInt_Pool08-1','RefInt_Pool09-1','RefInt_Pool10-1',
                             'RefInt_Pool11-1','RefInt_Pool12-1','RefInt_Pool13-1','RefInt_Pool14-1',
                             'RefInt_Pool15-1','RefInt_Pool16-1','RefInt_Pool17-1']
                # Drop quality control and ref intensity cols
                df = df.drop(drop_cols, axis = 'index')
                self._data["phosphoproteomics"] = df
            
            
            # prosp-brca-all-samples.txt shows which patient IDs have normal samples and which have replicates 
            # This file can be found on Box under CPTAC/cptac/pancan/helper_files
            if file_name == "prosp-brca-all-samples.txt":
                df = pd.read_csv(file_path, sep = "\t")
                df = df[['Participant', 'id', 'Type']]
                self._helper_tables["map_ids"] = df
                
            elif file_name == "README_v3.boxnote":
                self._readme_files["readme_proteomics"] = get_boxnote_text(file_path)
                
            elif file_name == "README.boxnote":
                self._readme_files["readme_phosphoproteomics"] = get_boxnote_text(file_path)

        print(' ' * len(loading_msg), end='\r') # Erase the loading message
        formatting_msg = f"Formatting {self.get_cancer_type()} dataframes..."
        print(formatting_msg, end='\r')
        
        # Get patient IDs with normal samples or replicates (from mapping file)         
        # 7 IDs with replicates: '11BR031', '11BR053', '11BR036', '11BR060', '14BR005', '11BR011', '21BR010'
        # 18 IDs with normals: '11BR074', '11BR073', '20BR007', '21BR010', '11BR017', '05BR029', '18BR003', '11BR030',
        #   '01BR027','11BR025', '11BR047', '11BR028', '11BR020', '20BR008', '11BR024', '11BR023', '11BR015', '11BR006'
        map_df = self._helper_tables["map_ids"]
        # Get IDs with replicates 
        replicate_list = list(set(map_df.loc[map_df.id.str.contains('REP')].Participant))
        replicate_list.remove('RetroIR') 
        replicate_list = [x[1:] for x in replicate_list]
        # Get IDs with normals 
        norm_df = map_df.loc[map_df.Type == 'Adjacent_Normal'] # get all patient_IDs with normal samples
        norm_df.index = norm_df.Participant.apply(lambda x: x[1:]+'.1') #remove initial 'X' and add '.1' (did not correlate well)
        not_tumor = norm_df.index.to_list()
        
        # Drop samples that don't correlate well with the original cptac tumor values
        # The proteomics and phosphoproteomics file contained duplicate patient IDs with unique values. 18 of these were IDs 
        # with both tumor and normal samples (according to the mapping file). To find which ID contained values nearest to the 
        # original cptac tumor measurements, we calculated pearson correlations and created scatterplots to compare 
        # flagship values to each of the duplicate patient ID values. The first occurrance of every patient ID in the 
        # file with a normal sample showed a high correlation with the flagship values (average around 0.9). 
        # The second occurrence showed a low correlation with the flagship tumor values (average around  0.2). 
        # In Brca, normal samples were dropped in downstream analysis because of quality control issues. Therefore, we dropped 
        # the second occurrence of the 18 patient IDs with a normal sample because they did not correlate well with the 
        # flagship tumor values and are likely normal samples. There were 7 IDs with replicates shown in the mapping file 
        # (prosp-brca-all-samples.txt). The same method was used to check that these correlated well with their respective 
        # flagship cptac tumor values. Replicates were averaged (consistent with the handling of other replicates in the pancan 
        # module).
        # note:  21BR010.1 had a correlation of 0.275 (so dropped), and 21BR010.2 had correlation of 0.848 (so averaged)        
        # A file containing the correlations can be downloaded at: 
        # https://byu.box.com/shared/static/jzsq69bd079oq0zbicw4w616hyicd5ev.xlsx
        
        if self._version == "1.0":
            # Proteomics
            prot = self._data["proteomics"]
            prot = prot.loc[ ~ prot.index.isin(not_tumor)] # drop rows that don't correlate well with respective cptac tumor 
            prot = average_replicates(prot, replicate_list) # average 7 IDs with replicates  
            self._data["proteomics"] = prot

            # Phosphoproteomics
            phos = self._data["phosphoproteomics"]
            phos = phos.loc[ ~ phos.index.isin(not_tumor)] # drop rows that don't correlate well with respective cptac tumor 
            phos = average_replicates(phos, replicate_list) # average 7 IDs with replicates
            self._data["phosphoproteomics"] = phos
        
        
        # Sort rows (tumor first then normal) and columns by first level (protein/gene name)
        self._data = sort_all_rows_pancan(self._data)  
        

        print(" " * len(formatting_msg), end='\r') # Erase the formatting message
