# CPTAC v.1.5 Release Candidate is now available
Please report any bugs or possible improvements. Enjoy!

# Easy access to CPTAC data
This software provides easy access to cancer data from the National Cancer Institute's CPTAC program, which characterizes and studies the proteogenomic landscape of tumors. We implement the software as a Python package called `cptac`, but you can seamlessly use it in an R environment with the help of the `reticulate` package (demonstrated in [Tutorial 6](https://paynelab.github.io/cptac/tutorial06_cptac_in_R.html)). Our package is installed in one step with `pip`:
```
pip install cptac
```
See the [Installation section below](https://paynelab.github.io/cptac/#installation) if you have further questions.

The package gives you the data as `pandas` `DataFrame` objects in Python. If you are using R, `reticulate` converts the tables to `data.frame` objects. By providing the tables natively in your programming environment, we eliminate the need for parsing and formatting, allowing you to quickly feed the data into whatever analysis code you have written. Follow our walkthrough tutorials and use cases for examples of how to use the software.

Additionally, the software automatically handles data downloading, storage, and updates. You need only to tell it which datasets you want downloaded, and it will automatically get the data without requiring you to write any HTTP requests or database queries.

## Installation
This package is intended to run on Python 3.6 or greater. If you plan on interfacing with it from R via `reticulate`, you must still have Python installed on your computer, and download the package into that Python environment.

### Installing Python
If you do not already have Python installed on your computer, we suggest using either the [standard Python distribution](https://www.python.org/downloads/) or the [Anaconda distribution](https://www.anaconda.com/distribution/). Follow the installation instructions at the respective links. The Anaconda distribution allows you to set up multiple distinct Python environments and comes with many useful Python packages pre-installed. For more information, see the Ananconda documentation.

### Installing the cptac package
We distribute the package [through the Python Package Index (PyPI)](https://pypi.org/project/cptac/), so regardless of which Python distribution you are using, you install the package using the `pip` program:
```
pip install cptac
```
If you are using the Anaconda distribution of Python, this will install `cptac` to the currently active environment as long as `pip` is available in that environment, which it would be by default. If `pip` is not installed in your environment, you can install it with `conda install -c anaconda pip`. Then, you can use `pip` to install the `cptac` package. We plan on making `cptac` directly available through `conda` in the near future.

The package depends on several other Python libraries including `numpy`, `pandas`, `requests`, and others. Normally, `pip` will automatically handle these dependencies when it installs `cptac` and you don't have to worry about any of it. However, if you have a special use case or are interested in exactly which versions of which packages are needed, you can consult the `install_requires` list in the [setup.py file](https://github.com/PayneLab/cptac/blob/master/setup.py).

## Documentation
Our goal is that our documentation will make this software and data accessible both to people without a computer science background, and people without a biology background. We provide two types of documentation to accomplish this: tutorials and use cases. The tutorials give a basic introduction to the software as well as conventions for storing and accessing the data. The use cases are short examples focused on a biological question and show practical uses of the software and data for biological discovery. Each use case works with a different combination of data types and explores meaningful cancer research hypotheses. 

You can access the tutorials and use cases as static webpages using the links below. They were originally written in Python as interactive Jupyter notebooks, so if you want to run them interactively with Jupyter you can download the notebooks from the [notebooks folder](https://github.com/PayneLab/cptac/tree/master/notebooks) on the GitHub repository. If you are unfamiliar with Jupyter, follow the installation and usage instructions given [here](https://jupyter.org/install) on the Jupyter website. You will then be able to run our tutorials as interactive, exploratory data analyses. If you want to run them interactively without installing anything, please visit our Binder site which hosts the notebooks [here](https://mybinder.org/v2/gh/PayneLab/cptac/master?filepath=%2Fnotebooks).

### Tutorials
- [Tutorial 1: CPTAC data introduction](https://paynelab.github.io/cptac/tutorial01_data_intro.html)
- [Tutorial 2: Using pandas to work with cptac dataframes](https://paynelab.github.io/cptac/tutorial02_pandas.html)
- [Tutorial 3: Joining dataframes with cptac](https://paynelab.github.io/cptac/tutorial03_joining_dataframes.html)
- [Tutorial 4: Understanding multi-indexes](https://paynelab.github.io/cptac/tutorial04_multiindex.html)
- [Tutorial 5: How to keep up to date with new package and data releases](https://paynelab.github.io/cptac/tutorial05_updates.html)
- [Tutorial 6: Easy integration with R](https://paynelab.github.io/cptac/tutorial06_cptac_in_R.html)

### Use cases
- [Use Case 1: Comparing transcriptomics and proteomics](https://paynelab.github.io/cptac/usecase01_omics.html)
- [Use Case 2: Correlation between clinical attributes](https://paynelab.github.io/cptac/usecase02_clinical_attributes.html)
- [Use Case 3: Associating clinical variables with omics data](https://paynelab.github.io/cptac/usecase03_clinical_and_acetylation.html)
- [Use Case 4: How Do Mutations Affect Protein Abundance?](https://paynelab.github.io/cptac/usecase04_mutations_and_omics.html)
- [Use Case 5: Gene Set Enrichment Analysis](https://paynelab.github.io/cptac/usecase05_enrichment_analysis.html)
- [Use Case 6: Comparing Derived Molecular Data with Proteomics](https://paynelab.github.io/cptac/usecase06_derived_molecular.html)
- [Use Case 7: Trans Genetics Effects](https://paynelab.github.io/cptac/usecase07_trans_genetic_effect.html)
- [Use Case 8: Outliers](https://paynelab.github.io/cptac/usecase08_outliers.html)
- [Use Case 9: Clinical Outcomes](https://paynelab.github.io/cptac/usecase09_clinical_outcomes.html)
- [Use Case 10: Pathway diagram overlay](https://paynelab.github.io/cptac/usecase10_pathway_overlay.html)

## Developer documentation
Documentation for anyone wanting to understand the internal workings of the package is available on the GitHub repository in the [devdocs folder](https://github.com/PayneLab/cptac/tree/master/devdocs).

## License
See the [LICENSE.md](https://github.com/PayneLab/cptac/blob/master/LICENSE.md) document on the GitHub repository. Please note the difference between the license as it applies to code versus data.

## Contact
This package is maintained by [the Payne lab](https://payne.byu.edu) at Brigham Young University.
