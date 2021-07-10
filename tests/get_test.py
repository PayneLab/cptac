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

class TestGet:

    # return an iterable list of <cptac.[Dataset]> objects
    @pytest.fixture(scope="class")
    def make_dataset_objects(self, get_all_datasets):
        """figure out how to pull specific attributes out of cptac to call them"""
        datasets_set = set(dir(cptac)).intersection(get_all_datasets)
        for dataset in datasets_set:
            try: 
                cptac.dataset()

            except exception:
                pytest.fail(f"cptac.{dataset} object was unable to be created")
        pass

    def test_all_getters(self):
        pass