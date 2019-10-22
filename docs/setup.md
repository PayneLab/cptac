
# Setup instructions

## Python
This data API is a Python package and is therefore only accessible via Python software. If you do not already have Python installed on your computer, we suggest either a download from Python (https://www.python.org/downloads/) or Anaconda (https://www.anaconda.com/distribution/). 

## Getting the API and data

The package is accessible <a href="https://pypi.org/project/cptac/">on the Python Package Index (PyPI)</a> and therefore installs with pip. You can download it from there with

```
pip install cptac
```

## Setting up for Jupyter notebooks
If you want the package available within Jupyter, you can install it in your base conda environment. To do this, open a bash shell and run the following:
```
conda activate
pip install cptac
```

## Opening the data in Python
After installation, open Python anywhere on your system and import our package by entering:

```
import cptac
```

The package should load automatically and print progress to the screen. See the tutorial examples in the doc folder.
