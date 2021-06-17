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

'''class for testing the loading of datasets'''
class TestLoad:

    def test_public_datasets(self, get_public_datasets):
        for dataset in get_public_datasets:
            # TODO: add way to see dataset-specific failures
            assert cptac.download(dataset)

    
    def test_protected_datasets(self, get_restricted_datasets):
        for dataset in get_restricted_datasets:
            # TODO: figure out how to handle passwords
            # assert cptac.download(dataset)
            assert False
    
