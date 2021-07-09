import pytest
import cptac

@pytest.fixture(autouse=True, scope="session")
def get_datasets_list():
    return cptac.list_datasets()

@pytest.fixture(scope="session")
def get_public_datasets(get_datasets_list):
    data = get_datasets_list["Data reuse status"]
    public_datasets = []
    for i in data.index:
        if data[i] == "no restrictions":
            public_datasets.append(i.lower())
    return public_datasets

@pytest.fixture(scope="session")
def get_restricted_datasets(get_datasets_list):
    data = get_datasets_list["Data reuse status"]
    restricted_datasets = []
    for i in data.index:
        if data[i] == "password access only":
            restricted_datasets.append(i.lower())
    return restricted_datasets

