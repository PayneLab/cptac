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
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import CPTAC

class Basic:
    def __init__(self):
        pass
    def evaluate_getters(self):
        print("Evaluating getters...")
        dataframes = []
        file_names = {}
        dataframes.append(CPTAC.get_clinical()); file_names[len(dataframes)] = "clinical"
        dataframes.append(CPTAC.get_proteomics()); file_names[len(dataframes)] = "proteomics"
        dataframes.append(CPTAC.get_phosphoproteomics()); file_names[len(dataframes)] = "site phosphoproteomics"
        dataframes.append(CPTAC.get_phosphoproteomics(gene_level = True)); file_names[len(dataframes)] = "gene phosphoproteomics"
        dataframes.append(CPTAC.get_transcriptomics()); file_names[len(dataframes)] = "linear transcriptomics"
        dataframes.append(CPTAC.get_transcriptomics(circular = True)); file_names[len(dataframes)] = "circular transcriptomics"
        dataframes.append(CPTAC.get_transcriptomics(miRNA = True)); file_names[len(dataframes)] = "miRNA"
        dataframes.append(CPTAC.get_CNA()); file_names[len(dataframes)] = "CNA"
        dataframes.append(CPTAC.get_somatic()); file_names[len(dataframes)] = "parsed somatic maf"
        dataframes.append(CPTAC.get_somatic(binary=True)); file_names[len(dataframes)] = "binary somatic"
        dataframes.append(CPTAC.get_somatic(unparsed=True)); file_names[len(dataframes)] = "unparsed somatic maf"
        PASS = True
        for x in range(0,len(dataframes)):
            if dataframes[x] is None:
                print("Error reading", file_names[x+1], "data")
                PASS = False
        if PASS:
            print("PASS")
        else:
            print("FAIL")
    def evaluate_special_getters(self):
        print("Evaluating special getters...")
        results = []
        functions = {}
        results.append(CPTAC.get_clinical_cols()); functions[len(results)] = "clinical_cols"
        results.append(CPTAC.get_cohort_clinical(["Diabetes","BMI"])); functions[len(results)] = "cohort_meta"
        results.append(CPTAC.get_proteomics_quant(["S018","S100"])); functions[len(results)] = "proteomics_quant"
        results.append(CPTAC.get_proteomics_cols()); functions[len(results)] = "proteomics_cols"
        results.append(CPTAC.get_transcriptomics_cols()); functions[len(results)] = "transcriptomics_cols"
        results.append(CPTAC.get_cohort_proteomics(["A1BG","TP53"])); functions[len(results)] = "cohort_proteomics"
        results.append(CPTAC.get_cohort_transcriptomics(["A1BG","TP53"])); functions[len(results)] = "cohort_transcriptomics"
        results.append(CPTAC.get_cohort_cna(["SASS6","TTTY22"])); functions[len(results)] = "cohort_cna"
        results.append(CPTAC.get_cohort_phosphoproteomics(["TP53-S315","AAAS-S541"])); functions[len(results)] = "cohort_phosphoproteomics"
        results.append(CPTAC.get_patient_mutations("C3L-00157")); functions[len(results)] = "patient_mutations(Patient_Id)"
        results.append(CPTAC.get_patient_mutations("S013")); functions[len(results)] = "patient_mutations(Clinical_Patient_Key)"
        results.append(CPTAC.get_phosphosites("TP53")); functions[len(results)] = "phosphosites"
        PASS = True
        for x in range(0,len(results)):
            if results[x] is None:
                print("Error with get",functions[x+1], "function")
                PASS = False
        if PASS:
            print("PASS")
        else:
            print("FAIL")
    def evaluate_utilities(self): #compare_**** functions
        print("Evaluating utilities...")
        results = []
        functions = {}
        results.append(CPTAC.compare_gene(CPTAC.get_proteomics(), CPTAC.get_transcriptomics(), "A1BG")); functions[len(results)] = "compare_gene"
        results.append(CPTAC.compare_gene(CPTAC.get_proteomics(), CPTAC.get_transcriptomics(), ["A1BG","RPL11"])); functions[len(results)] = "compare_genes"
        results.append(CPTAC.compare_clinical(CPTAC.get_proteomics(), "BMI")); functions[len(results)] = "compare_clinical"
        results.append(CPTAC.compare_mutations(CPTAC.get_proteomics(),"TP53")); functions[len(results)] = "compare_mutations(Proteomics)"
        results.append(CPTAC.compare_mutations(CPTAC.get_proteomics(),"TP53","AURKA")); functions[len(results)] = "compare_mutations(Proteomics with Somatic)"
        results.append(CPTAC.compare_mutations(CPTAC.get_phosphoproteomics(), "IRS2")); functions[len(results)] = "compare_mutations(Phosphoproteomics)"
        results.append(CPTAC.compare_mutations(CPTAC.get_phosphoproteomics(), "IRS2","PIK3CA")); functions[len(results)] = "compare_mutations(Phosphoproteomics with Somatic)"
        results.append(CPTAC.compare_mutations_full(CPTAC.get_proteomics(),"TP53")); functions[len(results)] = "compare_mutations_full(Proteomics)"
        results.append(CPTAC.compare_mutations_full(CPTAC.get_proteomics(),"TP53","AURKA")); functions[len(results)] = "compare_mutations_full(Proteomics with Somatic)"
        results.append(CPTAC.compare_mutations_full(CPTAC.get_phosphoproteomics(), "IRS2")); functions[len(results)] = "compare_mutations_full(Phosphoproteomics)"
        results.append(CPTAC.compare_mutations_full(CPTAC.get_phosphoproteomics(), "IRS2","PIK3CA")); functions[len(results)] = "compare_mutations_full(Phosphoproteomics with Somatic)"
        results.append(CPTAC.compare_phosphosites("TP53")); functions[len(results)] = "compare_phosphosites"
        PASS = True
        for x in range(0,len(results)):
            if results[x] is None:
                print("Error with",functions[x+1],"function")
                PASS = False
        if PASS:
            print("PASS")
        else:
            print("FAIL")


class Stats:
    def __init__(self):
        pass
    def evaluate(self, data, trait):
        data_trait = CPTAC.compare_clinical(data, trait)
        threshold = .05 / len(data.columns)
        tscutoff = .5
        significantTests = []
        significantGenes = []
        for num in range(1,len(data_trait.columns)):
            gene = data_trait.columns[num]
            oneGene = data_trait[[trait, gene]]
            oneGene = oneGene.dropna(axis=0)
            spearmanrTest = stats.spearmanr(oneGene[trait], oneGene[gene])
            if (abs(spearmanrTest[0]) >= tscutoff) and (spearmanrTest[1] <= threshold):
                significantTests.append(spearmanrTest)
                significantGenes.append(gene)
        if len(significantGenes) > 0:
            return "PASS"
        else:
            return "FAIL"
class Plotter:
    def __init__(self):
        pass
    def plot(self, data, column1, column2, method):
        if method == "scatterplot":
            plot = sns.relplot(x = column1, y = column2, data = data)
        elif method == "barplot":
            plot = sns.barplot(x = column1, y = column2, data = data)
        elif method == "boxplot":
            plot = sns.boxplot(x = column1, y = column2, data = data)
        else:
            message = method + " not a recognized method"
            print(message)
            return ""
        plt.show()

print("\nRunning tests:\n")

Basic().evaluate_getters()
Basic().evaluate_special_getters()
Basic().evaluate_utilities()

print("Plotting...")
Plotter().plot(CPTAC.get_proteomics(), "A1BG","PTEN","scatterplot")
Plotter().plot(CPTAC.get_clinical(), "Diabetes","BMI","barplot")
Plotter().plot(CPTAC.get_clinical(), "Diabetes","BMI","boxplot")
print("PASS")

print("Running statistics...")
message = Stats().evaluate(CPTAC.get_proteomics(), "Pathway_activity_p53")
print(message)

print("Version:",CPTAC.version())
