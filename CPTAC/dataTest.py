import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import CPTAC

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

print("Running tests:")

print("Plotting...")
Plotter().plot(CPTAC.get_proteomics(), "A1BG","PTEN","scatterplot")
Plotter().plot(CPTAC.get_clinical(), "Diabetes","BMI","barplot")
Plotter().plot(CPTAC.get_clinical(), "Diabetes","BMI","boxplot")
print("PASS")

print("Running statistics...")
message = Stats().evaluate(CPTAC.get_proteomics(), "Pathway_activity_p53")
print(message)
