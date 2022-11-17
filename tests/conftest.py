import cptac

def get_cancer_inputs(include_restricted=False):
    """
    Lists all available combinations of cancer-datatype-source.
    @param include_restricted (bool): Include restricted datasets as well.

    @return: List of tuples in the form [(cancer, dtype, source)].
    """
    datasets = cptac.list_datasets()
    datasets['Datatypes'] = datasets['Datatypes'].str.split(', ')
    datasets = datasets.explode('Datatypes')
    if not include_restricted:
        datasets = datasets[~datasets['Sources'].isin(['awgconf'])]
    return list(datasets.itertuples(index=False, name=None))


#
## Setting autouse=True here makes it so that this method always runs before any tests
#@pytest.fixture(scope="session", autouse=True)
#def get_datasets_lists():
#    '''
#    Returns: a dict of dataset lists
#        keys = ["public", "private"]
#    '''
#    #curses.wrapper
#    print(f"Getting dataset lists (public and private)...", end='\r')
#    # TODO list_datasets() no longer works as it used to.
#    # Need to understand what it used to get, and what it should
#    # return now.
#    data = cptac.list_datasets()["Data reuse status"]
#    public_datasets = []
#    restricted_datasets = []
#    for i in data.index:
#        if data[i] == "no restrictions":
#            public_datasets.append(i)
#        else:
#            restricted_datasets.append(i)
#            dataset_lists = {}
#    dataset_lists["public"] = public_datasets
#    dataset_lists["restricted"] = restricted_datasets
#
#    return dataset_lists
#
#### Download all datasets
## Must have autouse=True or else this never gets called
#@pytest.fixture(scope="session", autouse=True)
#def download_datasets(get_datasets_lists):
#    # Download public datasets
#    for cancer in get_datasets_lists["public"]:
#        try:
#            print(f"Downloading {cancer}...", end='\r')
#            cptac.download(cancer, redownload=True)
#        except:
#            pytest.fail(f"Unable to download data for {cancer} dataset.")
#
#    # TODO: Download restricted datasets
#
#    return True
#
