# import cptac
# import logging
# import pytest
# import pandas as pd
# import os.path as path

# logging.basicConfig(level=logging.INFO)

# def get_cancer_inputs(CPTAC_BASE_DIR):
#     # Get the index file
#     INDEX = pd.read_csv(path.join(CPTAC_BASE_DIR, 'data', 'index.tsv'), sep='\t')
    
#     # Split description column
#     INDEX['description_split'] = INDEX['description'].str.split('-')
    
#     # Get cancer, source and datatype from description
#     INDEX[['Source', 'Cancer', 'Datatype']] = pd.DataFrame(INDEX['description_split'].tolist(), index= INDEX.index)
    
#     # Create tuples (cancer, source, datatype, filename) for each row in the dataframe
#     cancer_inputs = list(INDEX[['Cancer', 'Source', 'Datatype', 'filename']].to_records(index=False))
    
#     return cancer_inputs

# #@pytest.fixture(scope="session", autouse=True)
# def get_datasets_lists():
#     '''
#     Gets the public and private dataset lists.
    
#     Returns: a dict of dataset lists with keys = ["public", "restricted"]
#     '''
#     logging.info("Getting datset lists (public and restricted)...")
#     data = cptac.list_datasets()["Data reuse status"]
#     public_datasets = [i for i in data.index if data[i] == "no restrictions"]
#     restricted_datasets = [i for i in data.index if data[i] != "no restrictions"]

#     return {"public": public_datasets, "restricted": restricted_datasets}

# #@pytest.fixture(scope="session", autouse=True)
# def download_datasets(get_datasets_lists):
#     """
#     Downloads all public and restricted datsets.

#     Returns: True upon successful completion.
#     """
#     # Download public datsets
#     for cancer in get_datasets_lists["public"]:
#         try:
#             logging.info(f"Downloading public {cancer} dataset...")
#             cptac.download(cancer, redownload=False)
#         except Exception as e:
#             logging.error(f"Unable to download data for {cancer} dataset. Error: {e}")
#             continue
    
#     # Download restricted datasets
#     for cancer in get_datasets_lists["restricted"]:
#         try:
#             logging.info(f"Downloading restricted {cancer} dataset...")
#             # Do we need additional authentification step or credentials to download restricted datasets?
#             cptac.download(cancer, redownload=False)
#         except Exception as e:
#             logging.error(f"Unable to download data for restricted {cancer} dataset. Error: {e}")
#     return True

