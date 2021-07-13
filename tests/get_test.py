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

import pytest
import cptac
from cptac.exceptions import DataFrameNotIncludedError

class TestGet:

    '''
    return an iterable dict of "cancer name" : <cptac.[Cancer Type]> entries
        for the cancer types that have publicly available data
    '''
    @pytest.fixture(scope="class")
    def make_public_cancer_type_objects(self, get_public_datasets):
        cancer_type_objects = {}

        for cancer in get_public_datasets:
            try: 
                method = getattr(cptac, cancer)
                cancer_type_objects[cancer] = method()
            except:
                pytest.fail(f"cptac.{cancer} object was unable to be created. Method: {method}")
        
        return cancer_type_objects

    '''
    returns dict of cancer types and their getters
    format:
    {
        "cancer type name" : {
            cancer type object : [list of getters]
        },
    }
    '''
    @pytest.fixture(scope="class")
    def get_cancers_and_dataset_getters(self, make_public_cancer_type_objects):
        cancer_types_and_getters = {}

        # start with dict of "name" : obj pairs
        cancer_types = make_public_cancer_type_objects.items()
        # loop through cancer objects
        for (name, cancer_obj) in cancer_types:
            single_cancer_and_getters = {} # will have <cptac.Cancer> : [cptac.Cancer.get_dataset] pairs
            getters = self.get_getters(cancer_obj)
            single_cancer_and_getters[cancer_obj] = getters
            cancer_types_and_getters[name] = single_cancer_and_getters

        return cancer_types_and_getters

    def get_getters(self, cancer_obj, data_list):
        # make list for dataset_obj's getters
        getters = []
        # loop through datasets in the current cancer type
        for dataset in cancer_obj.get_data_list():
            function_name = "get_" + dataset
            if function_name.__contains__("CNV"):
                function_name = "get_CNV"
            try: 
                g = getattr(cancer_obj, function_name)
                getters.append(g)
            except:
                pytest.fail(f"Unable to produce {dataset} getter for the {cancer_obj} cancer type.")

        return getters

    def test_all_getters(self, get_cancers_and_dataset_getters):
        cancers_and_getters = get_cancers_and_dataset_getters # rename for readability

        for cancer in cancers_and_getters.values():
            for (cancer_obj, getters) in cancer.items():
                for g in getters:
                    try:
                        if cancer_obj == cptac.UcecConf and g == cancer_obj.get_CNV:
                            g("log2ratio")
                            g("gistic")
                        else:
                            g()
                    except DataFrameNotIncludedError:
                        pytest.fail(DataFrameNotIncludedError)

    '''
    verify that no getters return a dataframe for datasets not in the cancer class
    '''
    def test_DataFrameNotIncludedError(self, make_public_cancer_type_objects):
        # use public cancer type objects
        cancer_types = make_public_cancer_type_objects.items()
        for (name, cancer_obj) in cancer_types:


        # find all getters not applicable to each cancer type
        # verify the correct error is thrown
        with DataFrameNotIncludedError:
            pass # call getter here
        pass