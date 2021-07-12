import pytest
import cptac

# returns a pandas dataframe of datasets plus extra info
@pytest.fixture(autouse=True, scope="session")
def get_datasets_list():
    return cptac.list_datasets()

# returns a list of dataset strings
@pytest.fixture(scope="session")
def get_public_datasets(get_datasets_list):
    data = get_datasets_list["Data reuse status"]
    public_datasets = []
    for i in data.index:
        if data[i] == "no restrictions":
            public_datasets.append(i)
    return public_datasets

# returns list of dataset strings
@pytest.fixture(scope="session")
def get_restricted_datasets(get_datasets_list):
    data = get_datasets_list["Data reuse status"]
    restricted_datasets = []
    for i in data.index:
        if data[i] == "password access only":
            restricted_datasets.append(i)
    return restricted_datasets

# return a list of dataset strings
@pytest.fixture(scope="session")
def get_all_datasets(get_public_datasets, get_restricted_datasets):
    return get_public_datasets + get_restricted_datasets