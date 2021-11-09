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
import pandas
import cptac
from cptac.exceptions import DataFrameNotIncludedError, InvalidParameterError

class TestGet:

    def test_valid_getters(self, get_cancer_test_units):
        test_units = get_cancer_test_units
        for cancer in test_units:
            for (getter_name, getter) in cancer.valid_getters.items():
                try:
                    if cancer.cancer_type == "UcecConf" and getter_name == "get_CNV":
                        dataframe = getter("log2ratio")
                        assert type(dataframe) == pandas.DataFrame
                        
                        dataframe = getter("gistic")
                        type(dataframe) == pandas.DataFrame
                    else:
                        dataframe = getter()
                        assert type(dataframe) == pandas.DataFrame
                
                except (DataFrameNotIncludedError, InvalidParameterError) as error:
                    pytest.fail(f"Calling {getter} resulted raised {error}")
                
                except:
                    pytest.fail(f"Calling {getter} caused error:\n\t{sys.exc_info()[0]}")

    def test_invalid_getters(self, get_cancer_test_units):
        test_units = get_cancer_test_units
        for cancer in test_units:
            for (getter_name, getter) in cancer.invalid_getters.items():
            # verify the correct error is thrown
                with pytest.raises(DataFrameNotIncludedError) as exception:
                    getter()