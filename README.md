# cptac
This project is intended to facilitate accessing and interacting with cancer data from the National Cancer Institute CPTAC consortium, which characterizes and studies the proteogenomic landscape of tumors. Currently, the datasets available are: endometrial cancer, ovarian cancer and colon cancer. These cancer studies are downloadable via our Python package as native dataframe objects and can therefore be integrated very quickly and easily with other Python-based data analysis tools. Follow our walkthrough tutorials for a demonstration of ways to use our system.

Setup instructions can be found in <code>doc/setup.md</code>

## Tutorials
Tutorials for this package describe how to use the package functions for research with the provided data. All the tutorials are written in Python using the interactive Jupyter notebooks. If you are unfamiliar with Jupyter, follow the instructions given at jupyter.org/install. You will then be able to run our tutorials as interactive, exploratory data analyses.
<ul>
  <li>Use Case 0: Exploring the data</li>
  <li>Use Case 1: Comparing transcriptomics and proteomics for a single gene</li>
  <li>Use Case 2: Looking for correlation between clinical factors</li>
  <li>Use Case 3: Find genes significantly correlated with a clinical attribute</li>
  <li>Use Case 4: Investigating how genetic mutation affects protein abundance</li>
  <li>Use Case 5: Running gene set enrichment analysis</li>
  <li>Use Case 6: Comparing derived molecular features with protein abundance</li>
</ul>

## Requirements
This package is intended to run on Python 3.6 with pandas 0.23.4. It also requires numpy and wget for Python. In the tutorials, we use seaborn 0.9.0 for data visualization. 

## License
This package contains LICENSE.md document which describes the license for use. Please note the difference between the license as it applies to code versus data.
