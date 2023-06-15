import pytest
import cptac

def test_get_dataset():
    """Test that the dataset instantiation works as expected"""
    try:
        lung = cptac.Luad()
    except Exception as e:
        pytest.fail(f"Could not instantiate Luad dataset: {str(e)}")

@pytest.mark.parametrize("dataset", ['Luad', 'Brca', 'Ccrcc'])
def test_multiple_datasets(dataset):
    """Test that multiple datasets can be instantiated"""
    try:
        data = getattr(cptac, dataset)()
    except Exception as e:
        pytest.fail(f"Could not instantiate {dataset} dataset: {str(e)}")
