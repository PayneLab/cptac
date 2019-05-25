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

import numpy as np
import pandas as pd
import os
import re
import glob
from .dataset import DataSet

class Colon(DataSet):

    def __init__(self):
        """Load all of the endometrial dataframes as values in the self.data dict variable, with names as keys, and format them properly."""

        # Call the parent DataSet __init__ function, which initializes self.data and other variables we need
        super().__init__()

        # Print welcome message
        message = "You have loaded the cptac colon dataset. To view available dataframes, call the dataset's list_data() method. To view available functions for accessing and manipulating the dataframes, call its list_api() method."
        wrapped_list = textwrap.wrap(message)
        for line in wrapped_list:
            print(line)

        # Print data version
        data_version = "Most recent release"
        print("colon data version: {}\n".format(data_version))

        # Get the path to the data files
        path_here = os.path.dirname(os.path.realpath(__file__))
        data_path = os.path.join(path_here, "data_colon", "*.*")
        files = glob.glob(data_path) # Put all files into a list

        # Load the data into dataframes in the self.data dict
        print("Loading cptac colon data:")
        for file in files: # Loops through files variable
            path_elements = path.split(os.sep) # Get a list of the levels of the path
            file_name = path_elements[-1] # The last element will be the name of the file
            file_name_split = file_name.split(".")
            df_name = file_name_split[0] # Our dataframe name will be the first section of file name (i.e. proteomics.txt.gz becomes proteomics)

            # Load the file, based on what it is
            print("Loading {} data...".format(df_name))
            if df_name == "mutation":
                df = pd.read_csv(file, sep="\t")
                df = df.sort_values(by="SampleID")
                df = df[["SampleID","Gene","Variant_Type","Protein_Change"]]
                df = df.rename({"Variant_Type":"Mutation","Protein_Change":"Location"},axis="columns")
                df.name = "somatic_" + df_name
            else:
                df = pd.read_csv(file, sep="\t",index_col=0)
                df = df.transpose()
                df.name = df_name
            self.data[df.name] = df # Maps dataframe name to dataframe. self.data was initialized when we called the parent class __init__()

        # Rename mutation_binary dataframe to somatic_mutation_binary
        df = self.data["mutation_binary"]
        df.name = "somatic_mutation_binary"
        self.data["somatic_mutation_binary"] = df
        del self.data["mutation_binary"]

        # Separate clinical and derived molecular dataframes
        all_clinical_data = self.data.get("clinical")
        clinical_df = all_clinical_data.drop(columns=['StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'TumorPurity','immuneSubtype', 'CIN', 'Integrated.Phenotype'])
        clinical_df.name = "clinical"
        derived_molecular_df = all_clinical_data[['StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'TumorPurity', 'immuneSubtype', 'CIN', 'Integrated.Phenotype']]
        derived_molecular_df.name = "derived_molecular"

        # Put them in our data dictionary
        self.data["clinical"] = clinical_df # Replaces original clinical dataframe
        self.data["derived_molecular"] = derived_molecular_df

        # Combine the two proteomics dataframes
        prot_tumor = self.data.get("proteomics_tumor")
        prot_normal = self.data.get("proteomics_normal") # Normal entries are marked with 'N' on the end of the ID
        prot_combined = prot_tumor.append(prot_normal)
        prot_combined.name = "proteomics"
        self.data[prot_combined.name] = prot_combined
        del self.data["proteomics_tumor"]
        del self.data["proteomics_normal"]

        # Get phosphoproteomics dataframes, so we can process and combine them
        phos_tumor = self.data.get("phosphoproteomics_tumor")
        phos_normal = self.data.get("phosphoproteomics_normal") # Normal entries are not marked

        # Mark entries in phosphoproteomics_normal dataframe with an N at the end of the ID, to match proteomics_normal
        phos_normal_indicies = phos_normal.index.values.tolist()
        for i in range(len(phos_normal_indicies)):
            index = phos_normal_indicies[i]
            index_marked = index + 'N'
            phos_normal_indicies[i] = index_marked
        new_phos_index = pd.Index(phos_normal_indicies)
        phos_normal = phos_normal.set_index(new_phos_index)

        # Combine the two phosphoproteomics dataframes into one dataframe.
        phos_combined = phos_tumor.append(phos_normal)
        phos_combined.name = 'phosphoproteomics'
        self.data[phos_combined.name] = phos_combined
        del self.data["phosphoproteomics_tumor"]
        del self.data["phosphoproteomics_normal"]

        # Rename the somamtic_mutation dataframe's "SampleID" column to "PatientID", then set that as the index, to match the other dataframes
        new_somatic = self.data["somatic_mutation"]
        new_somatic = new_somatic.rename(columns={"SampleID":"Patient_ID"})
        new_somatic = new_somatic.set_index("Patient_ID")
        new_somatic.name = "somatic_mutation"
        self.data["somatic_mutation"] = new_somatic

        # Get a union of all dataframes' indicies, with duplicates removed
        indicies = [df.index for df in self.data.values()]
        master_index = pd.Index([])
        for index in indicies:
            master_index = master_index.union(index)
            master_index = master_index.drop_duplicates()

        # Sort this master_index so all the samples with an N suffix are last. Because the N is a suffix, not a prefix, this is kind of messy.
        status_df = pd.DataFrame(master_index, columns=['Patient_ID']) # Create a new dataframe with the master_index as a column called "Patient_ID"
        status_col = []
        for index in master_index:
            if index[-1] == 'N':
                status_col.append("Normal")
            else:
                status_col.append("Tumor")
        status_df = status_df.assign(Status=status_col)
        status_df = status_df.sort_values(by=['Status', 'Patient_ID'], ascending=[False, True]) # Sorts first by status, and in descending order, so "Tumor" samples are first
        master_index = status_df["Patient_ID"].tolist()

        # Generate a sample ID for each patient ID
        sample_id_dict = {}
        for i in range(len(master_index)):
            patient_id = master_index[i]
            sample_id_dict[patient_id] = "S{:0>3}".format(i + 1) # Use string formatter to give each sample id the format S*** filled with zeroes, e.g. S001, S023, or S112

        # Put a mapping in the clinical dataframe of all patient ids to their sample ids, including patient ids for samples not originally in the clinical dataframe. 
        master_df = pd.DataFrame(index=master_index)
        master_clinical = self.data['clinical'].join(master_df, how='outer') # Do an outer join with the clinical dataframe, so that clinical has a row for every sample in the dataset
        master_clinical.name = self.data["clinical"].name

        # Add a column to clinical, Sample_Tumor_Normal, indicating whether each sample is a tumor or normal sample. Samples with a Patient_ID ending in N are normal.
        clinical_status_col = []
        for sample in master_clinical.index:
            if sample[-1] == 'N':
                clinical_status_col.append("Normal")
            else:
                clinical_status_col.append("Tumor")
        master_clinical.insert(1, "Sample_Tumor_Normal", clinical_status_col)

        # Replace the clinical dataframe in the data dictionary with our new and improved version!
        self.data['clinical'] = master_clinical 

        # Give the other dataframes Sample_ID indicies
        for name in self.data.keys(): # Only loop over keys, to avoid changing the structure of the object we're looping over
            df = self.data[name]
            sample_id_column = []
            for row in df.index:
                if row in sample_id_dict.keys():
                    sample_id_column.append(sample_id_dict[row])
                else:
                    print("Error mapping sample ids in {0} dataframe. Patient_ID {1} did not have corresponding Sample_ID mapped in clinical dataframe. {0} dataframe not loaded.".format(df.name, row))
                    continue
            df = df.assign(Sample_ID=sample_id_column)
            old_index_name = df.index.name
            if old_index_name is None:
                old_index_name = 'index' # If the current index doesn't have a name, the column it's put in when we do reset_index will have the default name of 'index'
            df = df.reset_index() # This gives the dataframe a default numerical index and makes the old index a column, which prevents it from being dropped when we set Sample_ID as the index.
            df = df.rename(columns={old_index_name:'Patient_ID'}) # Rename the old index as Patient_ID
            df = df.set_index('Sample_ID') # Make the Sample_ID column the index, which also drops the default numerical index from reset_index
            df = df.sort_index()
            df.name = df.name
            self.data[name] = df

        # Drop Patient_ID column from dataframes other than clinical. Keep it in clinical, so we have a mapping of Sample_ID to Patient_ID for every sample
        for name in self.data.keys():
            if name != "clinical":
                df = self.data[name]
                df = df.drop(columns="Patient_ID")
                self.data[name] = df
