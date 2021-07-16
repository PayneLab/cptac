import pytest
import cptac

# returns a dict of dataset lists
# key is accessibility: public or private
@pytest.fixture(scope="session")
def get_datasets_list():
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

    