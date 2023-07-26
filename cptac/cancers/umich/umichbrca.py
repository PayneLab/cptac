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


class UmichBrca(Source):
    def __init__(self, no_internet=False):
        """Define which dataframes as are available in the self.load_functions dictionary variable, with names as keys.

        Parameters:
        no_internet (bool, optional): Whether to skip the index update step because it requires an internet connection. This will be skipped automatically if there is no internet at all, but you may want to manually skip it if you have a spotty internet connection. Default is False.
        """

        # Set some needed variables, and pass them to the parent Dataset class __init__ function

        self.data_files = {
            "proteomics": "Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz",
            "phosphoproteomics": "Report_abundance_groupby=multi-site_protNorm=MD_gu=2.tsv.gz",
            # prosp-brca-all-samples.txt shows which patient IDs have normal samples and which have replicates
            "mapping": "prosp-brca-all-samples.txt.gz",
            "acetylproteomics": "abundance_multi-site_MD.tsv.gz",
            # "README_v3.boxnote" is proteomics
            # "README.boxnote" is phosphoproteomics 
            # "readme" : ["README_v3.boxnote", "README.boxnote"],
        }

        self.load_functions = {
            'phosphoproteomics': self.load_phosphoproteomics,
            'proteomics': self.load_proteomics,
            'acetylproteomics': self.load_acetylproteomics,
        }

        # Call the parent class __init__ function
        super().__init__(cancer_type="brca", source="umich", data_files=self.data_files,
                         load_functions=self.load_functions, no_internet=no_internet)

    def load_mapping(self):
        df_type = 'mapping'

        if not self._helper_tables:
            file_path = self.locate_files(df_type)

            map_df = pd.read_csv(file_path, sep="\t")
            map_df = map_df[['Participant', 'id', 'Type']]
            self._helper_tables["map_ids"] = map_df

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
            # flagship cptac tumor values. Replicates were averaged (consistent with the handling of other replicates in the pancan module).
            # note:  21BR010.1 had a correlation of 0.275 (so dropped), and 21BR010.2 had correlation of 0.848 (so averaged)
            # Get patient IDs with normal samples or replicates (from mapping file)
            # 7 IDs with replicates: '11BR031', '11BR053', '11BR036', '11BR060', '14BR005', '11BR011', '21BR010'
            # 18 IDs with normals: '11BR074', '11BR073', '20BR007', '21BR010', '11BR017', '05BR029', '18BR003', '11BR030',
            #   '01BR027','11BR025', '11BR047', '11BR028', '11BR020', '20BR008', '11BR024', '11BR023', '11BR015', '11BR006'

            # Get IDs with replicates
            replicate_list = list(set(map_df.loc[map_df.id.str.contains('REP')].Participant))
            replicate_list.remove('RetroIR')
            replicate_list = [x[1:] for x in replicate_list]
            self._helper_tables["replicate_list"] = replicate_list

            # Get IDs with normals
            norm_df = map_df.loc[map_df.Type == 'Adjacent_Normal']  # get all patient_IDs with normal samples
            norm_df.index = norm_df.Participant.apply(
                lambda x: x[1:] + '.1')  # remove initial 'X' and add '.1' (did not correlate well)
            not_tumor = norm_df.index.to_list()
            self._helper_tables["not_tumor"] = not_tumor

    def load_phosphoproteomics(self):
        df_type = 'phosphoproteomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep="\t")

            df_mapping = pd.read_csv(f"{CPTAC_BASE_DIR}/data/brca_mapping.csv")
            patient_dict = dict(zip(df_mapping['Hash'], df_mapping['Patient_ID']))

            # Parse a few columns out of the "Index" column that we'll need for our multiindex
            df[['Database_ID', 'Transcript_ID', "Gene_ID", "Havana_gene", "Havana_transcript", "Transcript", "Name",
                "Site"]] = df.Index.str.split("\\|", expand=True)
            df[['num1', 'start', "end", "detected_phos", "localized_phos", "Site"]] = df.Site.str.split("_",
                                                                                                        expand=True)

            # Some rows have at least one localized phosphorylation site, but also have other
            # phosphorylations that aren't localized. We'll drop those rows, if their localized
            # sites are duplicated in another row, to avoid creating duplicates, because we only
            # preserve information about the localized sites in a given row. However, if the localized
            # sites aren't duplicated in another row, we'll keep the row.

            # dectected_phos of the split "Index" column is number of phosphorylations detected,
            # and localized_phos is number of phosphorylations localized, so if the two values aren't equal, the row has at least one unlocalized site
            unlocalized_to_drop = df.index[~df["detected_phos"].eq(df["localized_phos"]) & df.duplicated(
                ["Name", "Site", "Peptide", "Database_ID"], keep=False)]
            df = df.drop(index=unlocalized_to_drop)

            df = df[df['Site'].notna()]  # only keep columns with phospho site
            df = df.set_index(['Name', 'Site', 'Peptide', 'Database_ID'])  # Create a multiindex in this order.
            # drop columns not needed in df
            df.drop(['Gene', "Index", "num1", "start", "end", "detected_phos", "localized_phos", "Havana_gene",
                     "Havana_transcript", "MaxPepProb", "Gene_ID", "Transcript_ID", "Transcript"], axis=1, inplace=True)
            df = df.transpose()
            ref_intensities = df.loc["ReferenceIntensity"]  # Get reference intensities to use to calculate ratios
            df = df.subtract(ref_intensities,
                             axis="columns")  # Subtract ref intensities from all the values, to get ratios
            df = df.iloc[1:, :]  # drop ReferenceIntensity row
            # drop ending of CPT retrospective samples to match cptac
            df = df.rename(index={'CPT000814-0004': 'CPT000814', 'CPT001846-0005': 'CPT001846'})

            drop_cols = ['RetroIR-07', 'RetroIR-13', 'RefInt_Pool01', 'RefInt_Pool02',
                         'RefInt_Pool03', 'RefInt_Pool04', 'RefInt_Pool05', 'RefInt_Pool06',
                         'RefInt_Pool07', 'RefInt_Pool08', 'RefInt_Pool09', 'RefInt_Pool10',
                         'RefInt_Pool11', 'RefInt_Pool12', 'RefInt_Pool13', 'RefInt_Pool14',
                         'RefInt_Pool15', 'RefInt_Pool16', 'RefInt_Pool17']
            # Drop quality control and ref intensity cols
            df = df.drop(drop_cols, axis='index')

            self.load_mapping()
            not_tumor = self._helper_tables["not_tumor"]
            replicate_list = self._helper_tables["replicate_list"]
            df = df.loc[~ df.index.isin(not_tumor)]  # drop rows that don't correlate well with respective cptac tumor
            df = df_tools.average_replicates(df, replicate_list)  # average 7 IDs with replicates

            df.index.name = 'Patient_ID'
            df = df.reset_index()

            df['Patient_ID'] = df['Patient_ID'].apply(lambda x: patient_dict.get(x, x))

            df['Patient_ID'] = df['Patient_ID'].apply(
                lambda x: x + '.N' if 'NX' in x else x)  # 'NX' are enriched normals
            df = df.set_index('Patient_ID')

            # save df in self._data
            self.save_df(df_type, df)

    def load_proteomics(self):
        df_type = 'proteomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)

            df = pd.read_csv(file_path, sep="\t")
            df['Database_ID'] = df["Index"].str.split('|').str[0]  # get protein identifier
            df['Name'] = df["Index"].str.split('|').str[6]  # get protein name
            df = df.set_index(['Name', 'Database_ID'])  # set multiindex
            df = df.drop(columns=['Index', 'MaxPepProb', 'NumberPSM', 'Gene'])  # drop unnecessary  columns
            df = df.transpose()
            ref_intensities = df.loc["ReferenceIntensity"]  # get reference intensities to use to calculate ratios
            df = df.subtract(ref_intensities, axis="columns")  # subtract reference intensities from all the values
            df = df.iloc[1:, :]  # drop ReferenceIntensity row
            df.index.name = 'Patient_ID'
            # drop ending of CPT retrospective samples to match cptac
            df = df.rename(index={'CPT0008140004': 'CPT000814', 'CPT0018460005': 'CPT001846',
                                  '604': 'CPT000814'})  # 604 mapped to CPT000814 in the proteomics file from the PDC pipeline

            drop_cols = ['RetroIR', 'RetroIR.1',
                         'RefInt_Pool01', 'RefInt_Pool02', 'RefInt_Pool03', 'RefInt_Pool04',
                         'RefInt_Pool05', 'RefInt_Pool06', 'RefInt_Pool07', 'RefInt_Pool08',
                         'RefInt_Pool09', 'RefInt_Pool10', 'RefInt_Pool11', 'RefInt_Pool12',
                         'RefInt_Pool13', 'RefInt_Pool14', 'RefInt_Pool15', 'RefInt_Pool16',
                         'RefInt_Pool17']
            df = df.drop(drop_cols, axis='index')  # drop quality control and ref intensity cols

            self.load_mapping()
            not_tumor = self._helper_tables["not_tumor"]
            replicate_list = self._helper_tables["replicate_list"]
            df = df.loc[~ df.index.isin(not_tumor)]  # drop rows that don't correlate well with respective cptac tumor
            df = df_tools.average_replicates(df, replicate_list)  # average 7 IDs with replicates

            # save df in self._data
            self.save_df(df_type, df)

    def load_acetylproteomics(self):
        df_type = 'acetylproteomics'

        if df_type not in self._data:
            # perform initial checks and get file path (defined in source.py, the parent class)
            file_path = self.locate_files(df_type)
            df = pd.read_csv(file_path, sep="\t")

            df_mapping = pd.read_csv(f"{CPTAC_BASE_DIR}/data/brca_mapping.csv")
            patient_dict = dict(zip(df_mapping['Hash'], df_mapping['Patient_ID']))

            # Parse a few columns out of the "Index" column that we'll need for our multiindex
            df[['Database_ID', 'Site1', "Site2", "Int1", "Int2", "Site"]] = df.Index.str.split("_", expand=True)
            df = df[df['Site'].notna()]  # only keep columns with phospho site

            # Load the gene names and merge them with the current dataframe based on 'Database_ID'
            df_gene_names = pd.read_csv(f"{CPTAC_BASE_DIR}/data/cptac_genes.csv")
            df_gene_names = df_gene_names.rename(columns={'Gene_Name': 'Name'})  # Renaming 'Gene_Name' to 'Name'
            df = pd.merge(df, df_gene_names, on='Database_ID', how='left')

            # Move 'Name' into the multiindex
            df = df.set_index(
                ['Name', 'Site', 'Peptide', 'Database_ID'])  # This will create a multiindex from these columns
            df = df.T  # transpose
            ref_intensities = df.loc["ReferenceIntensity"]  # Get reference intensities to use to calculate ratios
            df = df.iloc[1:, :]  # drop ReferenceIntensity row

            # Get dictionary with aliquots as keys and patient IDs as values
            self.load_mapping()
            mapping_dict = self._helper_tables["map_ids"]
            # df = df.rename(index = mapping_dict) # replace aliquots with patient IDs (normal samples have .N appended)
            # Add '.N' to enriched normal samples ('NX')
            df.index.name = 'Patient_ID'
            df = df.reset_index()

            df['Patient_ID'] = df['Patient_ID'].apply(lambda x: patient_dict.get(x, x))

            df['Patient_ID'] = df['Patient_ID'].apply(
                lambda x: x + '.N' if 'NX' in x else x)  # 'NX' are enriched normals
            df = df.set_index('Patient_ID')
            df = df_tools.rename_duplicate_labels(df,
                                                  'index')  # add ".1" to the second ocurrence of the ID with a duplicate
            # save df in self._data
            self.save_df(df_type, df)
#############################################

# TODO: Readmes
#             elif file_name == "README_v3.boxnote":
#                 self._readme_files["readme_proteomics"] = get_boxnote_text(file_path)

#             elif file_name == "README.boxnote":
#                 self._readme_files["readme_phosphoproteomics"] = get_boxnote_text(file_path)