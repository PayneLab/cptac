import pytest
import cptac
import pandas as pd

def test_dataset_format():
    """Test the dataset format"""
    brca = cptac.Brca()
    data = brca.get_proteomics(source='umich')
    assert isinstance(data, pd.DataFrame), "Data is not in DataFrame format"

def test_expected_output():
    """Test the dataset output"""
    brca = cptac.Brca()
    data = brca.get_proteomics(source='umich')
    expected_shape = (125, 12922)
    assert data.shape == expected_shape, f"Data shape mismatch. Expected {expected_shape}, got {data.shape}"
