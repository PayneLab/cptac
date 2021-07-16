import pytest
import cptac
'''
Setting autouse=True here makes it so that this method always runs before any tests
# returns a dict of dataset lists
# key is accessibility: "public" or "private"
'''
@pytest.fixture(scope="session", autouse=True)
def get_datasets_lists():
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

'''Download all datasets'''
@pytest.fixture(scope="session")
def download_datasets(get_datasets_lists):
    # Download public datasets
    for cancer in get_datasets_lists["public"]:
        try:
            cptac.download(cancer, redownload=True)
        except:
            pytest.fail(f"Unable to download data for {cancer} dataset.")

    # TODO: Download restricted datasets
            
    return True

@pytest.fixture(scope="session", autouse=True)
def get_public_dataset_objects(get_datasets_lists):
    cancer_dict = {}
    for cancer_data_set in get_datasets_lists["public"]:
        cancer = getattr(cptac, cancer_data_set)
        try:
            cancer_dict[cancer] = cancer()
        except:
            pytest.fail(f"unable to create {cancer} object")
        
    return cancer_dict, True
