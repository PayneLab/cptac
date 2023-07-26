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

class UmichHnscc(Source):
    def __init__(self, no_internet=False):
        """
        This class loads the dataset for the University of Michigan's Head and Neck Squamous Cell Carcinoma (HNSCC) study.

        Parameters:
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        self.data_files = {
            "proteomics" : "Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz",
            "phosphoproteomics" : "Report_abundance_groupby=multi-site_protNorm=MD_gu=2.tsv.gz",
            # "README_v3.boxnote" is proteomics
            # "README.boxnote" is phosphoproteomics 
            # "readme" : ["README_v3.boxnote", "README.boxnote"],
            #"not_used": "S039_BCprospective_observed_0920.tsv.gz",
            #"not_used": "S039_BCprospective_imputed_0920.tsv.gz"
        }
        
        self.load_functions = {
            'phosphoproteomics' : self.load_phosphoproteomics,
            'proteomics' : self.load_proteomics,
        }
        
        # Call the parent class __init__ function
        super().__init__(cancer_type="hnscc", source="umich", data_files=self.data_files, load_functions=self.load_functions, no_internet=no_internet)

    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

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
            
            # There were 4 labels with "-duplicate" appended in proteomics and phosphoproteomics files.
            # I ran a pearson correlation to check how well the values from each duplicate correlated to 
            # the other duplicates for the same case ID. Three of the duplicates correlated well with their 
            # respective case IDs. C3L-02617-N-duplicate2 did not correlate well with the other C3L-02617 duplicates, 
            # so we dropped it and averaged the other two. I also created a scatterplot to compare each duplicate to 
            # the first occurence of its case ID. The linear scatterplots indicated similarity between the aliquots. 
            # We averaged the duplicates that correlated well together and were the same tissue type.
            drop_cols = ['128C', 'QC2', 'QC3', 'QC4', '129N', 'LungTumor1', 'Pooled-sample14',
                       'LungTumor2', 'QC6', 'LungTumor3', 'Pooled-sample17', 'QC7',
                       'Pooled-sample19', 'QC9', 'RefInt_pool01', 'RefInt_pool02',
                       'RefInt_pool03', 'RefInt_pool04', 'RefInt_pool05', 'RefInt_pool06',
                       'RefInt_pool07', 'RefInt_pool08', 'RefInt_pool09', 'RefInt_pool10',
                       'RefInt_pool11', 'RefInt_pool12', 'RefInt_pool13', 'RefInt_pool14',
                       'RefInt_pool15', 'RefInt_pool16', 'RefInt_pool17', 'RefInt_pool18',
                       'RefInt_pool19', 'RefInt_pool20']
            phos = df
            phos = phos.drop(drop_cols, axis = 'index') # drop quality control and ref intensity cols        
            phos = phos.drop(['C3L-02617-N-duplicate2'], axis = 'index') # drop duplicate that did not correlate well
            # average IDs that correlated well to their respective duplicates
            phos = df_tools.average_replicates(phos, ['C3L-02617-T','C3L-02617-N','C3L-00994-N'], normal_identifier = '-N') 
            phos.index = phos.index.str.replace('-T$','', regex = True)
            phos.index = phos.index.str.replace('-N$','.N', regex = True)
            phos.index = phos.index.str.replace('-C$','.C', regex = True) # 6 cored normal samples in Hnscc
            df = phos
            
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
            
            # There were 4 labels with "-duplicate" appended in proteomics and phosphoproteomics files.
            # I ran a pearson correlation to check how well the values from each duplicate correlated to 
            # the other duplicates for the same case ID. Three of the duplicates correlated well with their 
            # respective case IDs. C3L-02617-N-duplicate2 did not correlate well with the other C3L-02617 duplicates, 
            # so we dropped it and averaged the other two. I also created a scatterplot to compare each duplicate to 
            # the first occurence of its case ID. The linear scatterplots indicated similarity between the aliquots. 
            # We averaged the duplicates that correlated well together and were the same tissue type.
            drop_cols = ['128C', 'QC2', 'QC3', 'QC4', '129N', 'LungTumor1', 'Pooled-sample14',
                       'LungTumor2', 'QC6', 'LungTumor3', 'Pooled-sample17', 'QC7',
                       'Pooled-sample19', 'QC9', 'RefInt_pool01', 'RefInt_pool02',
                       'RefInt_pool03', 'RefInt_pool04', 'RefInt_pool05', 'RefInt_pool06',
                       'RefInt_pool07', 'RefInt_pool08', 'RefInt_pool09', 'RefInt_pool10',
                       'RefInt_pool11', 'RefInt_pool12', 'RefInt_pool13', 'RefInt_pool14',
                       'RefInt_pool15', 'RefInt_pool16', 'RefInt_pool17', 'RefInt_pool18',
                       'RefInt_pool19', 'RefInt_pool20']
            prot = df
            prot = prot.drop(drop_cols, axis = 'index') # drop quality control and ref intensity cols        
            prot = prot.drop(['C3L-02617-N-duplicate2'], axis = 'index') # drop duplicate that did not correlate well  
            # These IDs had a high correlation with their respective duplicates, so we average them
            # duplicates: 'C3L-02617-T-duplicate', 'C3L-00994-N-duplicate', 'C3L-02617-N-duplicate'
            prot = df_tools.average_replicates(prot, ['C3L-02617-T','C3L-02617-N','C3L-00994-N'], normal_identifier = '-N') 
            prot.index = prot.index.str.replace('-T$','', regex = True)
            prot.index = prot.index.str.replace('-N$','.N', regex = True)
            prot.index = prot.index.str.replace('-C$','.C', regex = True) # 6 cored normal samples in Hnscc
            df = prot
            
            # save df in self._data
            self.save_df(df_type, df)

    def load_acetylproteomics(self):
        df_type = 'acetylproteomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep = "\t")
            # Parse a few columns out of the "Index" column that we'll need for our multiindex
            df[['Database_ID', "Site"]] = df.Index.str.split("_",expand=True)
            df = df[df['Site'].notna()] # only keep columns with phospho site

            # Load the gene names and merge them with the current dataframe based on 'Database_ID'
            df_gene_names = pd.read_csv(f"{CPTAC_BASE_DIR}/data/cptac_genes.csv")
            df_gene_names = df_gene_names.rename(columns={'Gene_Name': 'Name'}) # Renaming 'Gene_Name' to 'Name'
            df = pd.merge(df, df_gene_names, on='Database_ID', how='left')

            # Move 'Name' into the multiindex
            df = df.set_index(['Name', 'Site', 'Peptide', 'Database_ID']) # This will create a multiindex from these columns
            df = df.T # transpose
            ref_intensities = df.loc["ReferenceIntensity"]# Get reference intensities to use to calculate ratios
            df = df.iloc[1:,:] # drop ReferenceIntensity row

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
            df = df.drop('C3N-01825.1', axis = 'index') # drop the duplicate that didn't correlate well with flagship

            # save df in self._data
            self.save_df(df_type, df)

        
#############################################
# TODO: Readmes
#             elif file_name == "README_v3.boxnote":
#                 self._readme_files["readme_proteomics"] = get_boxnote_text(file_path)
                
#             elif file_name == "README.boxnote":
#                 self._readme_files["readme_phosphoproteomics"] = get_boxnote_text(file_path)
            
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
                self._data["proteomics_imputed"] = df
            '''