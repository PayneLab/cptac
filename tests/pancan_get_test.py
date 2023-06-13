#   Copyright 2018 Samuel Payne sam_payne@byu.edu
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# import pytest
# import pandas as pd
# import cptac
# from cptac.exceptions import DataFrameNotIncludedError, InvalidParameterError

# class TestPancanGet:

#     @pytest.fixture(autouse=True)
#     def setup(self):
#         self.cancer = cptac.Cancer() # Cancer class
#         pass

#     @pytest.mark.parametrize("valid_input,expected_result", [("Input1", "Expected1"), ("Input2", "Expected2")])
#     def test_valid_getters(self, valid_input, expected_result):
#         """
#         Test case for valid getters.
#         """
#         result = self.cancer.get(valid_input)
#         assert result == expected_result, f"For input: {valid_input}, expected{expected_result}, but got: {result}"

#     @pytest.mark.parametrize("invalid_input", ["Invalid1", "Invalid2"])
#     def test_invalid_getters(self, invalid_input):
#         """
#         Test for invalid getters.
#         """
#         with pytest.raises((DataFrameNotIncludedError, InvalidParameterError)):
#             self.cancer.get(invalid_input)