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

def test_all_datasets():
    datasets = cptac.list_datasets()
    for index, row in datasets.iterrows():
        cancer_class = get_cancer_class(row['Cancer'])
        cancer_instance = cancer_class()
        try:
            data = cancer_instance.get_dataframe(row['Datatype'], row['Source'])
            # Perform your checks here. For example, check that data is not empty:
            assert not data.empty, f"Data for {row['Cancer']} - {row['Source']} - {row['Datatype']} is empty."
        except cptac.exceptions.DataFrameNotIncludedError:
            # The dataset does not include this data. This is expected for some combinations, so we just pass.
            pass

def get_cancer_class(cancer_str):
    """Converts a string to a corresponding cancer class."""
    # This dictionary should be updated as necessary
    mapping = {
        "brca": cptac.Brca,
        "ccrcc": cptac.Ccrcc,
        "coad": cptac.Coad,
        "gbm": cptac.Gbm,
        "hnscc": cptac.Hnscc,
        "lscc": cptac.Lscc,
        "luad": cptac.Luad,
        "ov": cptac.Ov,
        "pdac": cptac.Pdac,
        "ucec": cptac.Ucec,
        "all_cancers" : cptac.Ucec
    }
    return mapping[cancer_str]