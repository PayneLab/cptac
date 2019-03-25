# CPTAC
This project is intended to facilitate accessing and interacting with cancer data from the National Cancer Institute CPTAC consortium, which characterizes and studies the proteogenomic landscape of tumors. Currently, the datasets available are: endometrial cancer, ovarian cancer and colon cancer. These cancer studies are downloadable via our Python package as native dataframe objects and can therefore be integrated very quickly and easily with other Python-based data analysis tools. Follow our walkthrough tutorials for a basic cookbook of ways to use our system.

Setup instructions can be found in <code>doc/setup.md</code>

## Tutorials
Tutorials for this package are located under the <code>doc</code> folder and describe how to use the package for research. All the tutorials are written in python using the interactive jupyter notebooks. If you are unfamiliar with jupyter, follow the instructions given at jupyter.org/install to begin running jupyter on your machine. You will then be able to run our documentation as interactive, jupyter-based tutorials.
<ul>
  <li>Use Case 1: Comparing transcriptomics and proteomics for a single gene</li>
<li>Use Case 2: Looking for correlation between various clinical factors, such as BMI, diabetes, and cancer stage</li>
<li>Use Case 3: Using statistical test to find genes significantly correlated with a clinical attribute</li>
<li>Use Case 4: Investigating how genetic mutation affects protein abundance</li>
<li>Use Case 5: Running gene set enrichment analysis - which pathways are differentially expressed based on clinical attribute</li>
</ul>

## Requirements
This package is intended to run on Python 3.6 with pandas 0.23.4. In the tutorials, we use seaborn 0.9.0 for data visualization. 

## License
This package contains LICENSE.md document which describes the license for use. Please note the difference between the license as it applies to code versus data.
