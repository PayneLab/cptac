import pytest
import cptac

@pytest.fixture
def get_all_datasets():
    return cptac.list_datasets()

@pytest.fixture
def get_public_datasets(get_all_datasets):
    public_datasets = []
    for dataset in get_all_datasets:
        if dataset.loc[dataset, "Data reuse stats"] == "no restrictions":
            public_datasets.append(dataset.lower())
    return public_datasets

@pytest.fixture
def get_restricted_datasets(get_all_datasets):
    restricted_datasets = []
    for dataset in get_all_datasets:
        if dataset.loc[dataset, "Data reuse stats"] == "password access only":
            restricted_datasets.append(dataset.lower())
    return restricted_datasets

@pytest.fixture
def test_public_datasets(get_public_datasets):
    for dataset in get_public_datasets:
        assert cptac.download(dataset)

@pytest.fixture
def test_protected_datasets(get_private_datasets):
    for dataset in get_private_datasets:
        assert cptac.download(dataset)