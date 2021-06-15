import pytest
import cptac
from requests.api import get

@pytest.fixture(scope="session")
def get_all_datasets():
    return cptac.list_datasets()

@pytest.fixture(scope="session")
def get_public_datasets(get_all_datasets):
    data = get_all_datasets
    public_datasets = []
    for dataset in data.index:
        if data.loc[dataset, "Data reuse stats"] == "no restrictions":
            public_datasets.append(dataset.lower())
    return public_datasets

@pytest.fixture(scope="session")
def get_restricted_datasets(get_all_datasets):
    data = get_all_datasets
    restricted_datasets = []
    for dataset in data.index:
        if data.loc[dataset, "Data reuse stats"] == "password access only":
            restricted_datasets.append(dataset.lower())
    return restricted_datasets

