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
import cptac.ovarian as ov

class Basic:
    def __init__(self):
        pass
    def evaluate_getters(self):
        print("Evaluating getters...")
        dataframes = []
        file_names = {}
        dataframes.append(ov.get_clinical()); file_names[len(dataframes)] = "clinical"
        dataframes.append(ov.get_cnv()); file_names[len(dataframes)] = "cnv"
        dataframes.append(ov.get_phosphoproteomics()); file_names[len(dataframes)] = "phosphoproteomics"
        dataframes.append(ov.get_proteomics()); file_names[len(dataframes)] = "proteomics"
        dataframes.append(ov.get_somatic_mutations()); file_names[len(dataframes)] = "somatic_19"
        dataframes.append(ov.get_transcriptomics()); file_names[len(dataframes)] = "transcripmics"
        PASS = True
        for x in range(0,len(dataframes)):
            if dataframes[x] is None:
                print("Error reading", file_names[x+1], "data")
                PASS = False
        if PASS:
            print("PASS")
        else:
            print("FAIL")

    def evaluate_comparers(self):

        # Load our data
        prot = ov.get_proteomics()
        phos = ov.get_phosphoproteomics()

        # Set our variables
        gene_1 = 'TP53'
        gene_2 = 'AURKA'

        # Test compare_mutations with one gene
        comp_1 = ov.compare_mutations(prot, gene_1)
        comp_2 = ov.compare_mutations(phos, gene_1)

        # Test compare_mutations with two genes. Should return omics data for first gene, with somatic mutation data for second gene
        comp_3 = ov.compare_mutations(prot, gene_1, gene_2)
        comp_4 = ov.compare_mutations(phos, gene_1, gene_2)

        # Test compare_mutations_full with one gene
        comp_5 = ov.compare_mutations_full(prot, gene_1)
        comp_6 = ov.compare_mutations_full(phos, gene_1)

        # Test compare_mutations_full with two genes. Should return omics data for first gene, with somatic mutation data for second gene
        comp_7 = ov.compare_mutations_full(prot, gene_1, gene_2)
        comp_8 = ov.compare_mutations_full(phos, gene_1, gene_2)

        print(comp_1)
        print(comp_2)
        print(comp_3)
        print(comp_4)
        print(comp_5)
        print(comp_6)
        print(comp_7)
        print(comp_8)

        print('PASS')

print("\nRunning tests:\n")

Basic().evaluate_getters()
Basic().evaluate_comparers()
