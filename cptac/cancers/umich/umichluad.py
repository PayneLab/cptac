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

from asyncio import CancelledError
import pandas as pd
import numpy as np
import os
import warnings

from cptac.cancers.source import Source
from cptac.tools.dataframe_tools import *
from cptac.exceptions import FailedReindexWarning, PublicationEmbargoWarning, ReindexMapError
from cptac.utils import get_boxnote_text


class UmichLuad(Source):

    def __init__(self, version="latest", no_internet=False):
        """Define which dataframes as are available in the self.load_functions dictionary variable, with names as keys.

        Parameters:
        version (str, optional): The version number to load, or the string "latest" to just load the latest datafreeze. Default is "latest".
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        # This keeps a record of all versions that the code is equipped to handle. That way, if there's a new data release but they didn't update their package, it won't try to parse the new data version it isn't equipped to handle.
        self.valid_versions = ["1.0"]

        self.data_files = {
            "1.0": {
                "proteomics" : "Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv",                    
                "phosphoproteomics" : "Report_abundance_groupby=multi-site_protNorm=MD_gu=2.tsv",
                "mapping" : "aliquot_to_patient_ID.tsv",
                # "README_v3.boxnote" is proteomics
                # "README.boxnote" is phosphoproteomics 
                "readme" : ["README_v3.boxnote", "README.boxnote"],            
            }
        }

        self.load_functions = {
            'phosphoproteomics' : self.load_phosphoproteomics,
            'proteomics' : self.load_proteomics,
        }

        if version == "latest":
            version = sorted(self.valid_versions)[-1]

        # Call the parent class __init__ function
        super().__init__(cancer_type="luad", source="umich", version=version, valid_versions=self.valid_versions, data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

        
    def load_mapping(self):
        df_type = 'mapping'

        if not self._helper_tables:
            file_path = self.locate_files(df_type)
            
            # aliquot_to_patient_ID.tsv contains only unique aliquots (no duplicates), 
            # so no need to slice out cancer specific aliquots
            # This file can be found on Box under CPTAC/cptac/pancan/helper_files
            df = pd.read_csv(file_path, sep = "\t", index_col = 'aliquot_ID', usecols = ['aliquot_ID', 'patient_ID'])
            map_dict = df.to_dict()['patient_ID'] # create dictionary with aliquot_ID as keys and patient_ID as values
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
            # that aren't localized. We'll drop those rows, if their localized sites are duplicated in another 
            # row, to avoid creating duplicates, because we only preserve information about the localized sites 
            # in a given row. However, if the localized sites aren't duplicated in another row, we'll keep the row.
            unlocalized_to_drop = df.index[~df["detected_phos"].eq(df["localized_phos"]) \
                                           & df.duplicated(["Name", "Site", "Peptide", "Database_ID"], keep=False)]
            # dectected_phos of the split "Index" column is number of phosphorylations detected, and localized_phos
            # is number of phosphorylations localized, so if the two values aren't equal, the row has at least 
            # one unlocalized site
            df = df.drop(index=unlocalized_to_drop)

            df = df[df['Site'].notna()] # only keep columns with phospho site 
            df = df.set_index(['Name', 'Site', 'Peptide', 'Database_ID']) # create a multiindex in this order.
            #drop columns not needed in df 
            df.drop(['Gene', "Index", "num1", "start", "end", "detected_phos", "localized_phos", "Havana_gene", 
                     "Havana_transcript", "MaxPepProb", "Gene_ID", "Transcript_ID", "Transcript"], axis=1, inplace=True)
            df = df.transpose()
            ref_intensities = df.loc["ReferenceIntensity"]# Get reference intensities to use to calculate ratios 
            df = df.subtract(ref_intensities, axis="columns") # Subtract ref intensities from all the values, to get ratios
            df = df.iloc[1:,:] # drop ReferenceIntensity row 
            
            # There were 2 duplicate IDs (aliquot mapped to the same tissue type and case ID) in the proteomic 
            # and phosphoproteomic data. I used the Payne lab mapping file "aliquot_to_patient_ID.tsv" to determine 
            # the tissue type for these duplicates. They were all tumor samples. Next, I ran a pearson correlation 
            # to check how well the values from each aliquot correlated to its respective tumor flagship sample. 
            # Each aliquot had a high correlation with the flagship values which indicates that they are replicates. 
            # I also created a scatterplot for each aliquot and flagship pair. The linear scatterplots indicated
            # similarity between the aliquot and flagship values. As the duplicate IDs were both tumor samples and 
            # correlated well with the flagship values, we averaged them.
            # A file containing the correlations can be downloaded at: 
            # https://byu.box.com/shared/static/jzsq69bd079oq0zbicw4w616hyicd5ev.xlsx

            # Drop quality control and ref intensity cols
            drop_cols = ['TumorOnlyIR01', 'NormalOnlyIR02', 'TumorOnlyIR03', 
                       'NormalOnlyIR04','NormalOnlyIR', 'TumorOnlyIR14',
                       'TaiwaneseIR19', 'TumorOnlyIR21', 'TaiwaneseIR22',
                       'NormalOnlyIR25', 'RefInt_pool01', 'RefInt_pool02', 'RefInt_pool03',
                       'RefInt_pool04', 'RefInt_pool05', 'RefInt_pool06', 'RefInt_pool07',
                       'RefInt_pool08', 'RefInt_pool09', 'RefInt_pool10', 'RefInt_pool11',
                       'RefInt_pool12', 'RefInt_pool13', 'RefInt_pool14', 'RefInt_pool15',
                       'RefInt_pool16', 'RefInt_pool17', 'RefInt_pool18', 'RefInt_pool19',
                       'RefInt_pool20', 'RefInt_pool21', 'RefInt_pool22', 'RefInt_pool23',
                       'RefInt_pool24', 'RefInt_pool25']
            df = df.drop(drop_cols, axis = 'index')
            
            # Get dictionary with aliquots as keys and patient IDs as values
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]
            df = df.rename(index = mapping_dict) # replace aliquots with patient IDs (normals have .N appended) 
            # manually map duplicates - these aliquots are in the mapping file, but they didn't map because of the appended ".1" 
            df = df.rename(index = {'CPT0146580004.1':'C3N-02379.1', 'CPT0148080004.1':'C3N-02587.1'}) 
            # these duplicates correlated well with their tumor flagship samples, so we average them
            df = average_replicates(df, ['C3N-02379', 'C3N-02587'])
            
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

            # There were 2 duplicate IDs (aliquot mapped to the same tissue type and case ID) in the proteomic
            # and phosphoproteomic data. I used the Payne lab mapping file "aliquot_to_patient_ID.tsv" to determine
            # the tissue type for these duplicates. They were all tumor samples. Next, I ran a pearson correlation
            # to check how well the values from each aliquot correlated to its respective tumor flagship sample.
            # Each aliquot had a high correlation with the flagship values which indicates that they are replicates.
            # I also created a scatterplot for each aliquot and flagship pair. The linear scatterplots indicated
            # similarity between the aliquot and flagship values. As the duplicate IDs were both tumor samples and
            # correlated well with the flagship values, we averaged them.
            # A file containing the correlations can be downloaded at:
            # https://byu.box.com/shared/static/jzsq69bd079oq0zbicw4w616hyicd5ev.xlsx

            # Drop quality control and ref intensity cols
            drop_cols = ['TumorOnlyIR01', 'NormalOnlyIR02', 'TumorOnlyIR03', 
                       'NormalOnlyIR04','NormalOnlyIR', 'TumorOnlyIR14',
                       'TaiwaneseIR19', 'TumorOnlyIR21', 'TaiwaneseIR22',
                       'NormalOnlyIR25', 'RefInt_pool01', 'RefInt_pool02', 'RefInt_pool03',
                       'RefInt_pool04', 'RefInt_pool05', 'RefInt_pool06', 'RefInt_pool07',
                       'RefInt_pool08', 'RefInt_pool09', 'RefInt_pool10', 'RefInt_pool11',
                       'RefInt_pool12', 'RefInt_pool13', 'RefInt_pool14', 'RefInt_pool15',
                       'RefInt_pool16', 'RefInt_pool17', 'RefInt_pool18', 'RefInt_pool19',
                       'RefInt_pool20', 'RefInt_pool21', 'RefInt_pool22', 'RefInt_pool23',
                       'RefInt_pool24', 'RefInt_pool25']
            df = df.drop(drop_cols, axis='index')
            
            # Get dictionary with aliquots as keys and patient IDs as values
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]
            df = df.rename(index = mapping_dict) # replace aliquots with patient IDs (normals have .N appended)
            # manually map duplicates - these aliquots are in the mapping file but they didn't map because of the appended ".1"
            df = df.rename(index = {'CPT0146580004.1':'C3N-02379.1', 'CPT0148080004.1':'C3N-02587.1'})
            # these duplicates correlated well with their tumor flagship samples, so we average them
            df = average_replicates(df, ['C3N-02379', 'C3N-02587'])

            # save df in self._data
            self.save_df(df_type, df)
        
#############################################



            # TODO take care of readmes
#             elif file_name == "README_v3.boxnote":
#                 text = get_boxnote_text(file_path)
#                 self._readme_files["readme_proteomics"] = text
                
#             elif file_name == "README.boxnote":
#                 text = get_boxnote_text(file_path)
#                 self._readme_files["readme_phosphoproteomics"] = text
