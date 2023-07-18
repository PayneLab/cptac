import pytest
import cptac
import pandas as pd

@pytest.mark.parametrize("dataset", ['Luad', 'Brca', 'Ccrcc'])
def test_consistency_across_calls(dataset):
    """Test for consistency across function calls"""
    data = getattr(cptac, dataset)()
    first_call = data.get_proteomics(source='umich')
    second_call = data.get_proteomics(source='umich')
    pd.testing.assert_frame_equal(first_call, second_call, check_dtype=True)
