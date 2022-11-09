import pytest
import cptac
import sys
from itertools import product

#sys.path.insert(0, "cptac/tests/") #pytest cwd is ?
from tests.test_cancer import TestCancer

# @pytest.fixture(scope="session")
# def cancer_source_dtype_combos():
#     # get options df from cptac
#     options = cptac.get_options()
#     options.drop(["Loadable datatypes"],  axis=1, inplace=True)
#     options.set_index("Cancers")

#     # convert Datatypes from type string to type list
#     options.Datatypes = options.Datatypes.apply(lambda x: x.split(", "))

#     # make a list of [(cancer, {source:datatype})] items
#     combos = options.apply(lambda x: [(x[0], {source: datatype}) for (source, datatype) in list(product([x[1]], x[2]))], axis=1)
    
#     # return a flattened list of (cancer, {source:datatype}) items
#     return [item for sublist in combos.to_list() for item in sublist]

# Setting autouse=True here makes it so that this method always runs before any tests
@pytest.fixture(scope="session", autouse=True)
def get_datasets_lists():
    '''
    Returns: a dict of dataset lists
        keys = ["public", "private"]
    '''
    #curses.wrapper
    print(f"Getting dataset lists (public and private)...", end='\r')
    # TODO list_datasets() no longer works as it used to.
    # Need to understand what it used to get, and what it should
    # return now.
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




# @pytest.fixture
# def get_all_datasets(scope="class"):
#     return cptac.list_datasets()

# @pytest.fixture
# def get_public_datasets(get_all_datasets):
#     public_datasets = []
#     for dataset in get_all_datasets:
#         if dataset.loc[dataset, "Data reuse stats"] == "no restricions":
#             public_datasets.append(dataset.lower())
#     return public_datasets

# @pytest.fixture
# def get_restricted_datasets(get_all_datasets):
#     restricted_datasets = []
#     for dataset in get_all_datasets:
#         if dataset.loc[dataset, "Data reuse stats"] == "password access only":
#             restricted_datasets.append(dataset.lower())
#     return restricted_datasets

# @pytest.fixture
# def test_public_datasets(get_public_datasets):
#     for dataset in get_public_datasets:
#         assert cptac.download(dataset)

# @pytest.fixture
# def test_protected_datasets(get_private_datasets):
#     for dataset in get_private_datasets:
#         assert cptac.download(dataset)