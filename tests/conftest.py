import pytest
import cptac
import sys
# TODO: figure out what cwd of pytest is
sys.path.insert(0, "cptac/tests/")
from tests.cancer import Cancer
import curses


# Setting autouse=True here makes it so that this method always runs before any tests
@pytest.fixture(scope="session", autouse=True)
def get_datasets_lists():
    '''
    Returns: a dict of dataset lists
        keys = ["public", "private"]
    '''
    #curses.wrapper
    print(f"Getting dataset lists (public and private)...", end='\r')
    data = cptac.list_datasets()["Data reuse status"]

    public_datasets = []
    restricted_datasets = []
    for i in data.index:
        if data[i] == "no restrictions":
            public_datasets.append(i)
        else:
            restricted_datasets.append(i)
    
    dataset_lists = {}
    dataset_lists["public"] = public_datasets
    dataset_lists["restricted"] = restricted_datasets
    
    return dataset_lists

### Download all datasets
# Must have autouse=True or else this never gets called
@pytest.fixture(scope="session", autouse=True)
def download_datasets(get_datasets_lists):
    # Download public datasets
    for cancer in get_datasets_lists["public"]:
        try:
            print(f"Downloading {cancer}...", end='\r')
            cptac.download(cancer, redownload=True)
        except:
            pytest.fail(f"Unable to download data for {cancer} dataset.")

    # TODO: Download restricted datasets
            
    return True


@pytest.fixture(scope="session", autouse=True)
def get_cancer_test_units(get_datasets_lists):
    '''
    Returns: a dict of this format:
        { "cancer name" : <Cancer Object>, ... }
    '''
    cancer_wrappers = list()
    for cancer_name in get_datasets_lists["public"]:
        c = getattr(cptac, cancer_name)
        try:
            print(f"Creating {c} object...", end='\r')
            cancer_wrappers.append(Cancer(cancer_name, c()))
        except:
            pytest.fail(f"unable to create {c} object")
        
    return cancer_wrappers
