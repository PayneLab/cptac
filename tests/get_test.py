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

    @pytest.fixture(scope='class')
    def all_getters(self):
        getters = set()
        for attribute in dir(cptac.dataset.Dataset):
            if attribute.startswith("get_"):
                getters.add(attribute)
        return getters
    
    @pytest.fixture(scope="class")
    def valid_getters(self, get_public_dataset_objects):
        '''@return a dict of str(cancers) : set(valid getter strings)'''
        cancer_sets = get_public_dataset_objects[0]
        valid_cancer_getters = {}
        for (cancer_type, cancer_object) in cancer_sets.items():
            valid_getters = set()
            for attribute in dir(cancer_type):
                if attribute.startswith("get_"):
                    for dataset in cancer_object.get_data_list():
                        function_name = "get_" + dataset
                        if attribute == function_name:
                            valid_getters.add(attribute)
            valid_cancer_getters[cancer_type] = valid_getters

        return valid_cancer_getters


    ''' Test Valid Getters '''

    def test_valid_getters(self, valid_getters, get_public_dataset_objects):
        # use cancer_sets dict {cptac.Cancer : cptac.Cancer instance}
        cancer_sets = get_public_dataset_objects[0]
        for (cancer_type, valid_getter_set) in valid_getters.items():
            for getter in valid_getter_set:
                g = getattr(cancer_sets[cancer_type], getter)
                try:
                    if cancer_type == cptac.UcecConf and getter == "get_CNV":
                        dataframe = g("log2ratio")
                        assert type(dataframe) == pandas.DataFrame
                        
                        dataframe = g("gistic")
                        type(dataframe) == pandas.DataFrame
                    else:
                        dataframe = g()
                        assert type(dataframe) == pandas.DataFrame
                
                except (DataFrameNotIncludedError, InvalidParameterError) as error:
                    pytest.fail(f"Calling {g} resulted raised {error}")
                
                except:
                    pytest.fail(f"Calling {g} caused error:\n\t{sys.exc_info()[0]}")


    ''' Test Invalid Getters '''
    
    @pytest.fixture(scope="class")
    def invalid_getters(self, all_getters, valid_getters):
        '''
        @return dict of str(cancers) : list(invalid getter strings)
        '''
        invalid_cancer_getters = {}
        for (cancer_type, valid_getter_list) in valid_getters.items():
            invalid_cancer_getters[cancer_type] = all_getters.difference(valid_getter_list)

        return invalid_cancer_getters

    def test_invalid_getters(self, invalid_getters, get_public_dataset_objects):
        # use cancer_sets dict {cptac.Cancer : cptac.Cancer instance}
        cancer_sets = get_public_dataset_objects[0]
        for (cancer_type, invalid_getter_set) in invalid_getters.items():
            for getter in invalid_getter_set:
                g = getattr(cancer_sets[cancer_type], getter)
            # verify the correct error is thrown
            with pytest.raises(DataFrameNotIncludedError) as exception:
                g()