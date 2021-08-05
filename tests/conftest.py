import pytest
import cptac
from cancer import Cancer
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

'''
Return a dict of this format:
{
    "cancer name" : <Cancer Object>,
    ...
}
'''
@pytest.fixture(scope="session", autouse=True)
def get_cancer_test_units(get_datasets_lists):
    cancer_wrappers = list()
    for cancer_name in get_datasets_lists["public"]:
        c = getattr(cptac, cancer_name)
        try:
            cancer_wrappers.append(Cancer(cancer_name, c()))
        except:
            pytest.fail(f"unable to create {c} object")
        
    return cancer_wrappers, True
