import pytest
import cptac
from cptac.exceptions import DataSourceNotFoundError

def test_error_handling():
    """Test the error handling for erroneous function calls"""
    brca = cptac.Brca()

    # Test a non-existent method
    with pytest.raises(AttributeError):
        brca.get_nonexistent_data()

    # Test with incorrect datatype
    with pytest.raises(DataSourceNotFoundError):
        brca.get_proteomics('invalid_datatype')

    # Test with incorrect source
    with pytest.raises(DataSourceNotFoundError):
        brca.get_dataframe("proteomics", source='nonexistent_source')

    # Test non-existent attribute
    with pytest.raises(AttributeError):
        brca.nonexistent_attribute