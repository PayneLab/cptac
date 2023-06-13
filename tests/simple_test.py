import pytest
import cptac
import pandas as pd

def test_import():
    """Test that cptac is importable"""
    try:
        import cptac
    except ImportError:
        pytest.fail("Could not import cptac")

def test_get_dataset():
    """Test that the dataset instantiation works as expected"""
    try:
        lung = cptac.Luad()
    except Exception as e:
        pytest.fail(f"Could not instantiate Luad dataset: {e}")

@pytest.mark.parametrize("dataset", ['Luad', 'Brca', 'Ccrcc'])
def test_multiple_datasets(dataset):
    """Test that multiple datasets can be instantiated"""
    try:
        data = getattr(cptac, dataset)()
    except Exception as e:
        pytest.fail(f"Could not instantiate {dataset} dataset: {e}")
def test_dataset_format():
    """Test the dataset format"""
    brca = cptac.Brca()
    data = brca.get_proteomics(source='umich')
    assert isinstance(data, pd.DataFrame)

def test_expected_output():
    """Test the dataset output"""
    brca = cptac.Brca()
    data = brca.get_proteomics(source='umich')
    expected_row_num, expected_column_num = 125, 12922
    assert data.shape == (expected_row_num, expected_column_num)

def test_error_handling():
    """Test the error handling for erroneous function calls"""
    with pytest.raises(AttributeError):
        brca = cptac.Brca()
        brca.get_nonexistent_data()

@pytest.mark.parametrize("dataset", ['Luad', 'Brca', 'Ccrcc'])
def test_consistency_across_calls(dataset):
    """Test for consistency across function calls"""
    data = getattr(cptac, dataset)()
    first_call = data.get_proteomics(source='umich')
    second_call = data.get_proteomics(source='umich')
    pd.testing.assert_frame_equal(first_call, second_call)

# @pytest.mark.parametrize("cancer", cptac.list_cancers())
# def test_all_cancers(cancer):
#     """Test that all cancers can be instantiated"""
#     try:
#         data = getattr(cptac, cancer)()
#     except Exception as e:
#         pytest.fail(f"Could not instantiate {cancer} dataset: {e}")

# @pytest.mark.parametrize("cancer", cptac.list_cancers())
# @pytest.mark.parametrize("source", cptac.list_sources())
# @pytest.mark.parametrize("datatype", cptac.list_datatypes())
# def test_all_sources_and_datatypes(cancer, source, datatype):
#     """Test that all sources and datatypes can be accessed for all cancers"""
#     try:
#         data = getattr(cptac, cancer)()
#         result = data.get_dataframe(datatype, source)
#         assert isinstance(result, pd.DataFrame)
#     except Exception as e:
#         pytest.fail(f"Could not get {datatype} data from source {source} for {cancer} dataset: {e}")