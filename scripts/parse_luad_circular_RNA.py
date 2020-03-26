#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# Purpose of this script: The LUAD dataset's circular RNA file takes a long time to parse into the format
# we want, and it uses a prohibitive amount of RAM to do so. In order to avoid these problems, we have
# written this script that will parse the LUAD circular RNA file separately and save the table to a tsv.
# This pre-parsed tsv can then be included with the data release instead of the unparsed file from the
# data freeze, so people don't have to spend all that time and RAM re-parsing the file every time they
# load the dataset. Normally we don't make any changes to the files from data freezes, and just
# distribute the original files through the package and handle parsing completely when we load the
# tables, but the extreme amount of time and RAM needed to parse this file made an exception to the rule
# seem like the best option.
#
# To use this script, download the unparsed circular RNA file from the data freeze and pass this script the
# path to the file. This script will read in the raw file, parse it, and save the parsed table to a tsv in
# the same directory as script. The output file will have the same name as the input file, but with
# "_parsed" appended to the file name. The output file will be a tsv, even if the input file was a csv.
# You can then distribute the parsed file with the package. Here is an example:
#
# python parse_luad_circular_RNA.py luad-v3.0-rnaseq-circ-rna.csv.gz
#
# Note: This script will automatically gzip the file it creates. Do not attempt to manually compress it again.

import pandas as pd
import numpy as np
import sys
import os

if len(sys.argv) < 2:
    raise ValueError("Must pass path to file to parse.")
elif len(sys.argv) > 2:
    raise ValueError("Too many arguments passed.")
else:
    file_path = sys.argv[1]

print(f"Parsing file '{file_path}'...")

split_path = file_path.split(os.sep)
file_name = split_path[-1]
split_name = file_name.split(".")

if split_name[-1] in ["gz", "zip"]:
    just_name = ".".join(split_name[0:-2])
    if split_name[-2] == "csv":
        cell_delim = ","
        print("Parsing as a csv...")
        
    elif split_name[-2] == "tsv":
        cell_delim = "\t"
        print("Parsing as a tsv...")
        
    else:
        cell_delim = "," # If all else fails, we'll assume it's a csv
        print("Parsing as a csv...")
else:
    just_name = ".".join(split_name[0:-1])
    if split_name[-1] == "csv":
        cell_delim = ","
        print("Parsing as a csv...")
        
    elif split_name[-1] == "tsv":
        cell_delim = "\t"
        print("Parsing as a tsv...")
        
    else:
        cell_delim = "," # If all else fails, we'll assume it's a csv
        print("Parsing as a csv...")

df = pd.read_csv(file_path, sep=cell_delim, dtype={"spanning.reads": np.int16}, engine="c")

junct_3_split = df['junction.3'].str.split(':', n=2, expand=True)
chrm = junct_3_split[0] # Get the chromosome
three_prime = junct_3_split[1] # Get the nucleotide coordinate of the last base of the acceptor

junct_5_split = df['junction.5'].str.split(':', n=2, expand=True)
five_prime = junct_5_split[1] # Get the nucleotide coordinates of the first base of the donor

# Now we need the gene name
diff = df['gene.5'] != df['gene.3'] # Create a boolean filter where genes are different
temp = df['gene.5'].where(diff, other="") # Replace the ones that are the same with an empty string
gene_name = temp + '_' + df["gene.3"] # Concatentate the temp column(which only has the genes from gene.5 that are different) to gene.3

# Put all those pieces of information together
df = df.assign(geneID=chrm + '_' + five_prime + '_' + three_prime + '_' + gene_name)

# Slice out the columns we want
df = df[['Sample.ID', 'geneID', 'spanning.reads']]

# There are about 3,000 duplicates in the file. Duplicate meaning that they have identical Sample IDs and identical geneID, but different spanning reads.
# Marcin Cieslik said to drop the one with the lowest spanning read.
df = df.sort_values(by='spanning.reads', ascending=False).drop_duplicates(['Sample.ID','geneID']).sort_index()
df = df.set_index("Sample.ID") # It's more efficient to save the table with the sample id index instead of a range index

new_file_name = just_name + "_parsed.tsv.gz"
df.to_csv(new_file_name, sep='\t', compression="gzip")

print(f"Parse complete. Parsed file saved to {new_file_name} in current directory.")
