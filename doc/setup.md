
# Setup instructions

## Downloading the data
### From GitHub
To download the data, open bash, navigate to the directory you want to store CPTAC under, and enter:

```
git clone 'http://github.com/PayneLab/CPTAC'
```

To install our package, we use pip. If you want the package available to command line interpreted python, open bash in the CPTAC directory and enter:

```
pip install .
```

### Directly from PyPI


The package will soon be accessible on pip itelf. For now, it's on the pip test site for debugging. You can download it from there with

```
pip install -i https://test.pypi.org/simple/ CPTAC
```

CPTAC is available in several versions with minor differences at the moment for testing purposes. Note that this is likely not the most recent distribution of the CPTAC data; for the latest release, please download the data from the GitHub page.
To download a different version of the data, add the version number after the CPTAC argument to pip install, like so:

```
pip install -i https://test.pypi.org/simple/ CPTAC==0.1.1
``` 

## Setting up for Jupyter notebooks
If you want the package available within Jupyter, you can install it in your base conda environment. To do this, open a bash shell and do the following:
```
conda activate
cd /path/to/CPTAC
pip install .
```

## Opening the data in Python
After installation, open Python anywhere on your system and use our method by entering:

```
import CPTAC
```

The data should load automatically and print progress to the screen. See the tutorial examples in the doc folder.
