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

import cptac

'''class for testing the loading of datasets'''
class TestLoad:
    # example test
    def test_brca(self):
        success = cptac.download("brca")
        assert success == True

    # trying to get automate the testing of all available datasets without having to update this test as more are added
    def test_all_cancer_types(self):
        datasets = cptac.list_datasets()
        # TODO get list of datasets out of datasets pandas dataframe
        # TODO convert all dataset names to lowercase
        # TODO figure out how to handle password protected datasets
        
        # How about something like this?
        for dataset in datasets.index:
            if datasets.loc[dataset, "Data reuse status"] == "password access only":
                # Still not sure how to test password protected files, but we can handle that here
                pass
            else:
                # Assert all non-password protected datasets can download
                # TODO: how do we identify individual failures dynamically?
                assert cptac.download(dataset.lower()) == True