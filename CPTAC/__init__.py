import webbrowser
import textwrap
import os
def list():
    """Parameters
    None

    Lists all available datasets

    Returns
    None
    """
    print("Available datasets:")
    print("Endometrial")
    print("Ovarian")
    print("Colon")
def help():
    """
    Parameters
    None

    Opens github help page

    Returns
    None
    """
    print("Opening help.txt in web browser...")
    webbrowser.open("https://github.com/PayneLab/CPTAC/blob/master/doc/help.txt")
def embargo():
    """
    Parameters
    None

    Opens CPTAC embargo details in web browser

    Returns
    None
    """
    print("Opening embargo details in web browser...")
    webbrowser.open("https://proteomics.cancer.gov/data-portal/about/data-use-agreement")
def version():
    """
    Parameters
    None

    Prints version number of CPTAC package

    Returns
    Version number
    """
    version = {}
    with open(dir_path + os.sep + "version.py") as fp:
    	exec(fp.read(), version)
    return(version['__version__'])

dir_path = os.path.dirname(os.path.realpath(__file__))
message = "Welcome to the CPTAC data service package. Available datasets may be viewed using CPTAC.list(). In order to access a specific data set, import a CPTAC subfolder using either \'import CPTAC.Dataset\' or \'from CPTAC import Dataset\'.\n"
wrapped_list = textwrap.wrap(message)
for line in wrapped_list:
    print(line)
print("******")
print("Version:",version())
print("******")
