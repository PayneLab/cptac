# cptac
This project is intended to facilitate accessing and interacting with cancer data from the National Cancer Institute CPTAC consortium, which characterizes and studies the proteogenomic landscape of tumors. These cancer studies are downloadable via our Python package as native dataframe objects and can therefore be integrated very quickly and easily with other Python-based data analysis tools. Follow our walkthrough tutorials for a demonstration of ways to use our system.

Our package is publicly available <a href="https://pypi.org/project/cptac/">on the Python Package Index (PyPI)</a>. Setup instructions can be found in <a href="https://github.com/PayneLab/cptac/blob/master/docs/setup.md">docs/setup.md</a>. The quick answer is that this Python package is installed via pip and then is available to all your python software. At a command line enter the following:

```
pip install cptac
```

## Documentation and Tutorials
All documentation for the project is contained in the <a href="https://github.com/PayneLab/cptac/tree/master/docs">docs/</a> folder. Documentation comes in two flavors: 1 - tutorials which are written to explain the data formats and conventions, 2 - use cases which are written to demonstrate exploring the data for biological discovery. All documentation is written in Python using the interactive Jupyter notebooks. If you are unfamiliar with Jupyter, follow the instructions given at <a href = "https://jupyter.org/install">jupyter.org/install</a>. You will then be able to run our tutorials as interactive, exploratory data analyses.

Separate documentation for software developers is available in the <a href="https://github.com/PayneLab/cptac/tree/master/devdocs">devdocs</a> folder.

## Requirements
This package is intended to run on Python 3.6 or greater with pandas version 0.25.0 or greater. In the tutorials, we use seaborn 0.9.0 for data visualization. 

## License
This package contains LICENSE.md document which describes the license for use. Please note the difference between the license as it applies to code versus data.

## Contact
This package is maintained by the Payne lab at Brigham Young University, https://payne.byu.edu
