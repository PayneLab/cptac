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

import os
import webbrowser
import textwrap
import glob
import pandas as pd
from .dataframe import DataFrameLoader
from .utilities import Utilities

dir_path = os.path.dirname(os.path.realpath(__file__))
data_directory = dir_path + os.sep + "Data" + os.sep
path = data_directory + "*.*"
files = glob.glob(path) #puts all files into iterable variable
data = {}
print("Loading Colon CPTAC data:")
for file in files: #loops through files variable
    try:
        df = DataFrameLoader(file).createDataFrame()
        data[df.name] = df #maps dataframe name to dataframe
    except IOError:
        print("Error reading", file)
        print("Check that all file names coincide with DataFrameLoader specs")

def list():
	"""
	Parameters: 
	None

	Prints a list of available dataframes and dimensions

	Returns:
	None
	"""
	print("Below are the available colon data frames contained in this package:")
	#dataframes = [clinical, miRNA, mutation, proteomics, transcriptomics, phosphoproteomics]
	for dataframe in data:
		print("\t", data[dataframe].name)
		print("\t", "\t", "Dimensions:", data[dataframe].shape)
	print("To access the data, use a get function with the data frame name, i.e. colon.get_proteomics()")

def get_clinical():
	"""
	Parameters: 
	None

	Returns:
	Clinical dataframe
	"""
	return data.get("clinical")
def get_miRNA():
	"""
	Parameters: 
	None

	Returns:
	miRNA dataframe
	"""
	return data.get("miRNA")
def get_mutation():
	"""
	Parameters: 
	None

	Returns:
	Mutation dataframe
	"""
	return data.get("mutation")
def get_phosphoproteomics():
	"""
	Parameters: 
	None

	Returns:
	Phosphoproteomics dataframe (both normal and tumor entries combined in one dataframe)
	"""
	tumor = data.get("phosphoproteomics_tumor")
	normal = data.get("phosphoproteomics_normal") #normal entries are not marked
	combined = tumor.append(normal)
	return combined
def get_proteomics():
	"""
	Parameters: 
	None

	Returns:
	Proteomics dataframe (both normal and tumor entries combined in one dataframe)
	"""
	tumor = data.get("proteomics_tumor")
	normal = data.get("proteomics_normal") #normal entries are marked with 'N' on the end of the ID
	combined = tumor.append(normal)
	return combined
def get_transcriptomics():
	"""
	Parameters: 
	None

	Returns:
	Transcriptomics dataframe
	"""
	return data.get("transcriptomics")
#TODO: add wrapper functions
"""
def compare_gene():
def compare_genes():
def compare_clinical():
def get_phosphosites():
def compare_phosphosites():
def add_mutation_hierarchy():
def merge_somatic():
def merge_mutations():
def merge_mutations_trans():
"""