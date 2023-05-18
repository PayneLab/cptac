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
from cptac.cancers.source import Source

class UmichPdac(Source):
    def __init__(self, no_internet=False):
        """Define which dataframes as are available in the self.load_functions dictionary variable, with names as keys.

        Parameters:
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        self.data_files = {
            "proteomics" : "Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv",
            "phosphoproteomics" : "Report_abundance_groupby=multi-site_protNorm=MD_gu=2.tsv",
            "mapping" : "aliquot_to_patient_ID.tsv",
            # "README_v3.boxnote" is proteomics
            # "README.boxnote" is phosphoproteomics 
            "readme" : ["README_v3.boxnote", "README.boxnote"],
        }
        
        self.load_functions = {
            'phosphoproteomics' : self.load_phosphoproteomics,
            'proteomics' : self.load_proteomics,
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="pdac", source="umich", data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

    def load_mapping(self):
        df_type = 'mapping'

        if not self._helper_tables:
            file_path = self.locate_files(df_type)
            
            # aliquot_to_patient_ID.tsv contains only unique aliquots (no duplicates), 
            # so there is no need to slice out cancer specific aliquots
            # This file can be found on Box under CPTAC/cptac/pancan/helper_files
            df = pd.read_csv(file_path, sep = "\t", index_col = 'aliquot_ID', usecols = ['aliquot_ID', 'patient_ID'])
            map_dict = df.to_dict()['patient_ID'] # create dictionary with aliquots as keys and patient IDs as values
            self._helper_tables["map_ids"] = map_dict

    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep = "\t") 
            # Parse a few columns out of the "Index" column that we'll need for our multiindex
            df[['Database_ID','Transcript_ID',"Gene_ID","Havana_gene","Havana_transcript","Transcript","Name","Site"]] = df.Index.str.split("\\|",expand=True)
            df[['num1','start',"end","detected_phos","localized_phos","Site"]] = df.Site.str.split("_",expand=True) 

             # Some rows have at least one localized phosphorylation site, but also have other phosphorylations 
            # that aren't localized. We'll drop those rows, if their localized sites are duplicated in another row, 
            # to avoid creating duplicates, because we only preserve information about the localized sites in a 
            # given row. However, if the localized sites aren't duplicated in another row, we'll keep the row.
            unlocalized_to_drop = df.index[~df["detected_phos"].eq(df["localized_phos"]) & \
                                           df.duplicated(["Name", "Site", "Peptide", "Database_ID"], keep=False)]
            # dectected_phos of the split "Index" column is number of phosphorylations detected, and 
            # localized_phos is number of phosphorylations localized, so if the two values aren't equal, 
            # the row has at least one unlocalized site
            df = df.drop(index=unlocalized_to_drop)

            df = df[df['Site'].notna()] # only keep columns with phospho site 
            df = df.set_index(['Name', 'Site', 'Peptide', 'Database_ID']) # create a multiindex in this order
            #drop columns not needed in df 
            df.drop(['Gene',  "Index", "num1", "start", "end", "detected_phos", "localized_phos", "Havana_gene", "Havana_transcript", "MaxPepProb", "Gene_ID", "Transcript_ID", "Transcript"], axis=1, inplace=True)

            df = df.T #transpose df 
            ref_intensities = df.loc["ReferenceIntensity"]# Get reference intensities to use to calculate ratios 
            df = df.subtract(ref_intensities, axis="columns") # Subtract ref intensities from all the values, to get ratios
            # Drop qauality control and ref intensity 
            drop_cols = ['ReferenceIntensity', 'QC1', 'QC2', 'QC3', 'QC4', 'QC5', 'QC6', 'KoreanReference1',
                       'KoreanReference2', 'KoreanReference3', 'Pool-24-2', 'WU-PDA1', 'WU-Pool-25','RefInt_pool-01',
                     'RefInt_pool-02','RefInt_pool-03','RefInt_pool-04','RefInt_pool-05','RefInt_pool-06','RefInt_pool-07',
                     'RefInt_pool-08','RefInt_pool-09', 'RefInt_pool-10','RefInt_pool-11','RefInt_pool-12','RefInt_pool-13',
                     'RefInt_pool-14','RefInt_pool-15','RefInt_pool-16','RefInt_pool-17','RefInt_pool-18','RefInt_pool-19',
                     'RefInt_pool-20','RefInt_pool-21','RefInt_pool-22','RefInt_pool-23','RefInt_pool-24','RefInt_pool-25']
            df = df.drop(drop_cols, axis = 'index')

            # These 8 aliquots were not in the mapping file. Yize said they are all normal samples.
            manually_mapped = {'CPT0347760002': 'C3L-07032.N', 'CPT0347790002': 'C3L-07033.N',
                'CPT0347820002': 'C3L-07034.N', 'CPT0347850002': 'C3L-07035.N', 'CPT0347880002': 'C3L-07036.N',
                'CPT0355180003': 'C3L-03513.N', 'CPT0355190003': 'C3L-03515.N', 'CPT0355200003': 'C3L-03514.N'}

            # Get dictionary to map aliquots to patient IDs
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]

            df = df.rename(index = mapping_dict) # replace aliquots with patient IDs (normals have .N) 
            df = df.rename(index = manually_mapped) # map 8 aliquots that were not in the mapping file

            # save df in self._data
            self.save_df(df_type, df)

    def load_proteomics(self):
        df_type = 'proteomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            
            df = pd.read_csv(file_path, sep = "\t") 
            df['Database_ID'] = df.Index.apply(lambda x: x.split('|')[0]) # get protein identifier 
            df['Name'] = df.Index.apply(lambda x: x.split('|')[6]) # get protein name 
            df = df.set_index(['Name', 'Database_ID']) # set multiindex
            df = df.drop(columns = ['Index', 'MaxPepProb', 'NumberPSM', 'Gene']) # drop unnecessary  columns
            df = df.transpose()
            ref_intensities = df.loc["ReferenceIntensity"] # get reference intensities to use to calculate ratios 
            df = df.subtract(ref_intensities, axis="columns") # subtract reference intensities from all the values
            df.index.name = 'Patient_ID'
            # Drop quality control and ref intensity 
            drop_cols = ['ReferenceIntensity', 'QC1', 'QC2', 'QC3', 'QC4', 'QC5', 'QC6', 'KoreanReference1',
               'KoreanReference2', 'KoreanReference3', 'Pool-24-2', 'WU-PDA1', 'WU-Pool-25']
            df = df.drop(drop_cols, axis = 'index')
            
            # These 8 aliquots were not in the mapping file. Yize said they are all normal samples.
            manually_mapped = {'CPT0347760002': 'C3L-07032.N', 'CPT0347790002': 'C3L-07033.N',
                'CPT0347820002': 'C3L-07034.N', 'CPT0347850002': 'C3L-07035.N', 'CPT0347880002': 'C3L-07036.N',
                'CPT0355180003': 'C3L-03513.N', 'CPT0355190003': 'C3L-03515.N', 'CPT0355200003': 'C3L-03514.N'}

            # Get dictionary to map aliquots to patient IDs
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]

            df = df.rename(index = mapping_dict) # replace aliquots with patient IDs (normals have .N)
            df = df.rename(index = manually_mapped) # map 8 aliquots that were not in the mapping file

            # save df in self._data
            self.save_df(df_type, df)
        
        
#############################################

#             elif file_name == "README_v3.boxnote":
#                 self._readme_files["readme_proteomics"] = get_boxnote_text(file_path)
                
#             elif file_name == "README.boxnote":
#                 self._readme_files["readme_phosphoproteomics"] = get_boxnote_text(file_path)
