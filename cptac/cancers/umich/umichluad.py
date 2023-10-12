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
import cptac.tools.dataframe_tools as df_tools
from cptac import CPTAC_BASE_DIR

class UmichLuad(Source):
    """
    This class handles data loading for the Lung Adenocarcinoma (LUAD) dataset from University of Michigan.
    """
    def __init__(self, no_internet=False):
        """
        Initialize the UmichLuad object. It initializes data_files and load_functions which contains the name of the
        dataframes and their corresponding load functions. This information will be passed to the parent Source class.

        Parameters:
        no_internet (bool, optional): Determines whether to skips the index update step, useful when there's no internet connection.
        """

        # Define the available data files
        self.data_files = {
            "proteomics" : "Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz",                    
            "phosphoproteomics" : "Report_abundance_groupby=multi-site_protNorm=MD_gu=2.tsv.gz",
            "acetylproteomics" : "abundance_multi-site_MD.tsv.gz",
            "mapping" : "aliquot_to_patient_ID.tsv.gz",          
        }

        # Define the load functions for the different data types
        self.load_functions = {
            'phosphoproteomics' : self.load_phosphoproteomics,
            'proteomics' : self.load_proteomics,
            'acetylproteomics' : self.load_acetylproteomics,
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="luad", source="umich", data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

    def load_mapping(self):
        """
        Loads the mapping file which maps aliquot_ID to patient_ID. This file is needed for correct
        identification and correlation of data between different types of dataframes (e.g., proteomics,
        phosphoproteomics, etc.)
        """
        df_type = 'mapping'

        if not self._helper_tables:
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep = "\t", index_col = 'aliquot_ID', usecols = ['aliquot_ID', 'patient_ID'])
            map_dict = df.to_dict()['patient_ID'] # create dictionary with aliquot_ID as keys and patient_ID as values
            self._helper_tables["map_ids"] = map_dict

    def load_phosphoproteomics(self):
        """
        Loads the phosphoproteomics data. It processes the data to a clean, usable format and stores it in 
        the _data attribute
        """
        df_type = 'phosphoproteomics'

        if df_type not in self._data:
            # perform initial checks and get file path 
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
            df = df_tools.average_replicates(df, ['C3N-02379', 'C3N-02587'])
            
            # save df in self._data
            self.save_df(df_type, df)

    def load_proteomics(self):
        """
        Loads the proteomics data. It processes the data to a clean, usable format and stores it in
        the _data attribute.
        """
        df_type = 'proteomics'

        if df_type not in self._data:
            # perform initial checks and get file path 
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
            df = df_tools.average_replicates(df, ['C3N-02379', 'C3N-02587'])

            # save df in self._data
            self.save_df(df_type, df)

    def load_acetylproteomics(self):
        """
        Loads the acetylproteomics data. It processes the data to a clean, usable format and stores it in
        the _data attribute.
        """
        df_type = 'acetylproteomics'

        if df_type not in self._data:
            # perform initial checks and get file path
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep = "\t")
            # Parse a few columns out of the "Index" column that we'll need for our multiindex
            df[['Database_ID','Site1',"Site2","Int1","Int2", "Site"]] = df.Index.str.split("_",expand=True)
            df = df[df['Site'].notna()] # only keep columns with phospho site

            # Load the gene names and merge them with the current dataframe based on 'Database_ID'
            df_gene_names = pd.read_csv(f"{CPTAC_BASE_DIR}/data/cptac_genes.csv")
            df_gene_names = df_gene_names.rename(columns={'Gene_Name': 'Name'}) # Renaming 'Gene_Name' to 'Name'
            df = pd.merge(df, df_gene_names, on='Database_ID', how='left')

            # Move 'Name' into the multiindex
            df = df.set_index(['Name', 'Site', 'Peptide', 'Database_ID']) # This will create a multiindex from these columns
            df.drop(["Gene", "Int1", "Int2", "MaxPepProb", "ProteinID", "ReferenceIntensity", "Site1", "Site2", "Index"], axis=1, inplace=True)
            df = df.T # transpose
            #ref_intensities = df.loc["ReferenceIntensity"]# Get reference intensities to use to calculate ratios
            #df = df.iloc[1:,:] # drop ReferenceIntensity row

            # There was 1 duplicate ID (C3N-01825) in the proteomic and phosphoproteomic data.
            # I used the Payne lab mapping file "aliquot_to_patient_ID.tsv" to determine the tissue type
            # for these duplicates, and they were both tumor samples. Next, I ran a pearson correlation
            # to check how well the values from each duplicate correlated to its tumor flagship sample.
            # The first occurrence in the file had a higher correlation with the flagship sample
            # than the second occurrence. I also created scatterplots comparing each duplicate to its flagship sample.
            # We dropped the second occurrence of the duplicate because it didn't correlate very well to its flagship sample.

            # Get dictionary with aliquots as keys and patient IDs as values
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]
            df = df.rename(index = mapping_dict) # replace aliquots with patient IDs (normal samples have .N appended)
            # Add '.N' to enriched normal samples ('NX')
            df.index.name = 'Patient_ID'
            df = df.reset_index()
            df['Patient_ID'] = df['Patient_ID'].apply(lambda x: x+'.N' if 'NX' in x else x) # 'NX' are enriched normals
            df = df.set_index('Patient_ID')
            df = df_tools.rename_duplicate_labels(df, 'index') # add ".1" to the second ocurrence of the ID with a duplicate
            if 'C3N-01825.1' in df.index:
                df = df.drop('C3N-01825.1', axis = 'index') # drop the duplicate that didn't correlate well with flagship

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
