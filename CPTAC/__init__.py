import webbrowser
import textwrap
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
message = "Welcome to the CPTAC data service package. This import contains information about the package. In order to access a specific data set, import a CPTAC subfolder by either \'import CPTAC.DataName\' or \'from CPTAC import DataName\'."
wrapped_list = textwrap.wrap(message)
for line in wrapped_list:
    print(line)

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
    with open(dir_path + os.sep + "version.py") as fp: #.. required to navigate up to CPTAC folder from Endometrial folder, TODO: how to navigate from dataTest.py?
    	exec(fp.read(), version)
    return(version['__version__'])
