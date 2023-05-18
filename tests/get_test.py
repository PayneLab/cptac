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

import sys
import pytest
import pandas as pd
import cptac
from cptac.exceptions import DataFrameNotIncludedError, InvalidParameterError

#from conftest import get_cancer_inputs
#
class TestGet:

    # from cancer test units, grab the dict of function names and function pointers for valid getters and test them to make sure they return a dataframe object
    @pytest.mark.parametrize("cancer", get_cancer_inputs())
    def test_valid_getters(self, cancer):
        for getter_name, getter in cancer.valid_getters.items():
            try:
                if cancer.cancer_type == "UcecConf" and getter_name == "get_CNV":
                    dataframe = getter("log2ratio")
                    assert isinstance(dataframe, pd.DataFrame)
                    
                    dataframe = getter("gistic")
                    assert isinstance(dataframe, pd.DataFrame)
                else:
                    dataframe = getter()
                    assert isinstance(dataframe, pd.DataFrame)
                    
            except (DataFrameNotIncludedError, InvalidParameterError) as error:
                pytest.fail(f"Calling {getter} resulted raised {error}")
                
            except Exception as e:
                pytest.fail(f"Calling {getter} caused error:\n\t{type(e).__name__}")

    @pytest.mark.parametrize("cancer", get_cancer_inputs())
    def test_invalid_getters(self, cancer):
        for getter_name, getter in cancer.invalid_getters.items():
            # verify the correct error is thrown
            with pytest.raises(DataFrameNotIncludedError):
                getter()