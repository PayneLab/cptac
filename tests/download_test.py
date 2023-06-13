# #   Copyright 2018 Samuel Payne sam_payne@byu.edu
# #   Licensed under the Apache License, Version 2.0 (the "License");
# #   you may not use this file except in compliance with the License.
# #   You may obtain a copy of the License at
# #       http://www.apache.org/licenses/LICENSE-2.0
# #   Unless required by applicable law or agreed to in writing, software
# #   distributed under the License is distributed on an "AS IS" BASIS,
# #   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# #   See the License for the specific language governing permissions and
# #   limitations under the License.

# import os
# import time
# import cptac
# from cptac import CPTAC_BASE_DIR
# import logging
# import pytest
# from cptac.exceptions import InvalidParameterError

# from .conftest import get_cancer_inputs

# logging.basicConfig(level=logging.INFO)

# def validate_download(cancer, source):
#     """Helper function to check if the directory exists."""
#     download_path = f"{cptac.CPTAC_BASE_DIR}/data/data_{source}_{cancer}"
#     return os.path.isdir(download_path)

# @pytest.mark.parametrize("cancer, source, datatype, data_file", get_cancer_inputs(CPTAC_BASE_DIR))
# def test_redownload(cancer, source, datatype, data_file):
#     logging.info(f"Testing download of {cancer}, {source}, {datatype}, {data_file}")
#     assert cptac.download(cancer=cancer, source=source, dtype=datatype, data_file=data_file)

# @pytest.mark.parametrize("cancer, source, datatype, data_file", get_cancer_inputs(CPTAC_BASE_DIR))
# def test_invalid_parameter(cancer, source, datatype, data_file):
#     """Test the function with invalid parameters. This should raise an exception."""
#     with pytest.raises(InvalidParameterError):
#         cptac.download(sources={source: datatype}, cancers="invalid_cancer", redownload=True)
    
