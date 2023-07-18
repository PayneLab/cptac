import pytest
import sys

def test_import():
    """Test that cptac is importable"""
    try:
        import cptac
    except ImportError:
        pytest.fail("Could not import cptac")