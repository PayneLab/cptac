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
from cptac import CPTAC_BASE_DIR

class UmichCcrcc(Source):
    """
    This class handles the loading of the University of Michigan's
    ccRCC dataset into the CPTAC data structure.
    """

    def __init__(self, no_internet=False):
        """
        Initializes the class with a set of predefined dataframes
        that are available in the self.load_functions dictionary
        variable, with names as keys

        Parameters:
        no_internet (bool, optional): Whether to skip the index update step 
        because it requires an internet connection. This will be skipped 
        automatically if there is no internet at all, but you may want to 
        manually skip it if you have a spotty internet connection. Default is False.
        """

        # Define the data files associated with this dataset
        self.data_files = {
            "proteomics" : "Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz",
            "phosphoproteomics" : "Report_abundance_groupby=multi-site_protNorm=MD_gu=2.tsv.gz",
            "mapping" : "aliquot_to_patient_ID.tsv.gz",
            # "README_v3.boxnote" is proteomics
            # "README.boxnote" is phosphoproteomics 
            # "readme" : ["README_v3.boxnote", "README.boxnote"],
            #"not_used": "S039_BCprospective_observed_0920.tsv.gz",
            #"not_used": "S039_BCprospective_imputed_0920.tsv.gz"
        }
        
        # Define the functions used to load each type of data
        self.load_functions = {
            'phosphoproteomics' : self.load_phosphoproteomics,
            'proteomics' : self.load_proteomics,
        }
        
        # Call the parent class __init__ function
        super().__init__(cancer_type="ccrcc", source="umich", 
                         data_files=self.data_files, 
                         load_functions=self.load_functions, 
                         no_internet=no_internet)

    def load_mapping(self):
        """
        Load the mapping file to map aliquot IDs to patient IDs.
        """
        df_type = 'mapping'
        if not self._helper_tables:
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep = "\t", index_col = 'aliquot_ID', 
                             usecols = ['aliquot_ID', 'patient_ID'])
            map_dict = df.to_dict()['patient_ID'] # create dictionary with aliquot_ID as keys and patient_ID as values
            self._helper_tables["map_ids"] = map_dict
            self._helper_tables["drop_cols"] = ['NCI7-1','NCI7-2','NCI7-3','NCI7-4','NCI7-5', 'QC1', 'QC2', 'QC3', 'QC4', 'QC5', 'QC6', 'QC7', 
                    'QC8', 'RefInt_pool01', 'RefInt_pool02', 'RefInt_pool03', 'RefInt_pool04', 'RefInt_pool05', 
                    'RefInt_pool06', 'RefInt_pool07', 'RefInt_pool08', 'RefInt_pool09', 'RefInt_pool10', 'RefInt_pool11', 
                    'RefInt_pool12', 'RefInt_pool13', 'RefInt_pool14', 'RefInt_pool15', 'RefInt_pool16', 'RefInt_pool17', 
                    'RefInt_pool18', 'RefInt_pool19', 'RefInt_pool20', 'RefInt_pool21', 'RefInt_pool22', 'RefInt_pool23']

    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'
        
        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep = "\t")
            
            # Parse a few columns out of the "Index" column that we'll need for our multiindex
            df[['Database_ID','Transcript_ID',"Gene_ID","Havana_gene","Havana_transcript","Transcript","Name","Site"]] = df.Index.str.split("\\|",expand=True)
            df[['num1','start',"end","detected_phos","localized_phos","Site"]] = df.Site.str.split("_",expand=True) 

            # Some rows have at least one localized phosphorylation site, but also have other
            # phosphorylations that aren't localized. We'll drop those rows, if their localized sites
            # are duplicated in another row, to avoid creating duplicates, because we only preserve 
            # information about the localized sites in a given row. However, if the localized sites aren't 
            # duplicated in another row, we'll keep the row.
            unlocalized_to_drop = df.index[~df["detected_phos"].eq(df["localized_phos"]) & \
                                           df.duplicated(["Name", "Site", "Peptide", "Database_ID"], keep=False)]
            # dectected_phos of the split "Index" column is number of phosphorylations detected, and 
            # localized_phos is number of phosphorylations localized, so if the two values aren't equal, 
            # the row has at least one unlocalized site
            df = df.drop(index=unlocalized_to_drop)

            df = df[df['Site'].notna()] # only keep columns with phospho site 
            df = df.set_index(['Name', 'Site', 'Peptide', 'Database_ID']) # create a multiindex in this order
            #drop columns not needed in df 
            df.drop(['Gene', "Index", "num1", "start", "end", "detected_phos", "localized_phos", "Havana_gene", 
                     "Havana_transcript", "MaxPepProb", "Gene_ID","Transcript_ID", "Transcript"], axis=1, inplace=True)
            df = df.transpose()
            ref_intensities = df.loc["ReferenceIntensity"]# Get reference intensities (prep to calculate ratios) 
            df = df.subtract(ref_intensities, axis="columns") # Subtract refintensities from all the values, to get ratios
            df = df.iloc[1:,:] # drop ReferenceIntensity row 
            
            self.load_mapping()
            # see mapping for what exactly this is
            drop_cols = self._helper_tables["drop_cols"]
            map_ids = self._helper_tables["map_ids"]
            df = df.drop(drop_cols, axis = 'index') # drop quality control and ref intensity cols
            df = df.rename(index = map_ids) # replace aliquot_IDs with Patient_IDs (normal samples have .N appended)
            
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
            df = df.iloc[1:,:] # drop ReferenceIntensity row 
            df.index.name = 'Patient_ID'

            self.load_mapping()
            # see mapping for what exactly this is
            drop_cols = self._helper_tables["drop_cols"]
            map_ids = self._helper_tables["map_ids"]
            df = df.drop(drop_cols, axis = 'index') # drop quality control and ref intensity cols
            df = df.rename(index = map_ids) # replace aliquot_IDs with Patient_IDs (normal samples have .N appended)

            # save df in self._data
            self.save_df(df_type, df)


#############################################

# TODO: Readmes
#             elif file_name == "README_v3.boxnote":
#                 self._readme_files["readme_proteomics"] = get_boxnote_text(file_path)
                
#             elif file_name == "README.boxnote":
#                 self._readme_files["readme_phosphoproteomics"] = get_boxnote_text(file_path)
      
            '''
            elif file_name == "S039_BCprospective_observed_0920.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.transpose()
                df.index.name = 'Patient_ID'
                df.columns.name = 'Name'
                df = average_replicates(df)
                df = df.sort_values(by=["Patient_ID"])
                self._data["proteomics"] = df  
                
            elif file_name == "S039_BCprospective_imputed_0920.tsv.gz":
                df = pd.read_csv(file_path, sep="\t")
                df = df.transpose()
                df.index.name = 'Patient_ID'
                df.columns.name = 'Name'
                df = average_replicates(df)
                df = df.sort_values(by=["Patient_ID"])
                self._data["proteomics_imputed"] = df
            '''