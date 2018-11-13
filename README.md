# CPTAC
This project is intended to facilitate accessing and interacting with cancer data from the NIH. Currently, the dataset available is from a series of endometrial cancer studies. These cancer studies are downloadable via our Python package as native dataframe objects from the pandas package and can therefore be integrated very quickly and easily with other Python-based data analysis tools. Follow our walkthrough tutorials for a basic cookbook of ways to use our system.

Setup instructions can be found in <code>doc/setup.md</code>

## Tutorials
Tuturials for this package are located under the <code>doc</code> folder. They are all in jupyter notebooks; if you are unfamiliar with jupyter, follow the instructions given at jupyter.org/install to begin running jupyter on your machine. You will then be able to run our documentation as interactive, jupyter-based tutorials.
Use Case 1: Comparing transcriptomics and proteomics for a given gene
Use Case 2: Comparing various clinical factors, such as BMI, diabetes, and cancer stage, to look for correlation
Use Case 3: Using Spearman correlation to find genes significantly correlated with a given clinical attribute
Use Case 4: Investigating protein abundance for different types of genetic mutations
Use Case 5: Running gene set enrichment analysis - which gene pathways are significantly expressed based on microsatellite instability

## Requirements
This package is intended to run on Python 3.6 with pandas 0.23.4 and the seaborn package for data visualization. 
