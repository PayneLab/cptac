# CPTAC
Python packaging for CPTAC data

To download the data, open bash, navigate to the directory you want to store CPTAC under, and enter:

```
git clone 'http://github.com/PayneLab/CPTAC'
```

To install our package, we use pip. IF you want the package available to command line interpreted python, open bash in the CPTAC directory and enter:

```
pip install .
```

If you want the package available within Jupyter, you should open a bash shell and do the following:
```
conda activate
cd /path/to/CPTAC
pip install .
```

After installation, open Python anywhere on your system and use our method by entering:

```
import CPTAC
```

The data should load automatically and print progress to the screen. See the tutorial examples in the doc folder.
