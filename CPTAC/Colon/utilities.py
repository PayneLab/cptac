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
class Utilities:
    def __init__(self):
        pass

    def compare_gene(self, df1, df2, gene):
        pass

    def compare_genes(self, df1, df2, genes):
        pass

    def compare_clinical(self, clinical, data, clinical_col):
        pass

    def compare_phosphosites(self, proteomics, phosphoproteomics, gene):
        pass

    def add_mutation_hierarchy(self, somatic): # private
        pass

    def merge_somatic(self, somatic, gene, df_gene, multiple_mutations = False): # private
        pass

    def merge_mutations(self, omics, somatic, gene, duplicates = False):
        pass

    def merge_mutations_trans(self, omics, omics_gene, somatic, somatic_gene, duplicates = False):
        pass
