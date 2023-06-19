# import pytest
# import cptac
# from .utils import get_cancer_class

# datasets = cptac.list_datasets()

# @pytest.mark.parametrize("dataset", datasets.iterrows())
# def test_all_datasets(dataset):
#     """Test all datasets available in cptac"""
#     index, row = dataset
#     cancer_class = get_cancer_class(row['Cancer'])
#     cancer_instance = cancer_class()
#     try:
#         data = cancer_instance.get_dataframe(row['Datatype'], row['Source'])
#         assert not data.empty, f"Data for {row['Cancer']} - {row['Source']} - {row['Datatype']} is empty."

#         # Add more assertions here if desired

#     except cptac.exceptions.DataFrameNotIncludedError:
#         pass  # Expected for some combinations
#     except cptac.exceptions.DataSourceNotFoundError:
#         pass  # Expected if the data source doesn't exist
#     except AttributeError:
#         pass  # Expected if the get_dataframe function doesn't exist
#     except Exception as e:
#         pytest.fail(f"Unexpected error occurred: {str(e)}")


