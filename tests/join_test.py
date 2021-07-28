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
import itertools
import cptac

"""
TODO: Things that should happen in a join test:
    - check for correct error throwing with non-existant data (columns, tables, etc)
    - check for correct indexing and table sizing
    - check slicing for correct table sizing
    - check non-overlapping rows
    - 
"""
# TODO: Look through Caleb's stuff on tests
# TODO: Check use cases for standard usage and then try to mess that up
class TestJoin:

    # @pytest.fixture(scope="class")
    # def all_joiners(self):
    #     joiners = set()
    #     for attribute in dir(cptac.dataset.Dataset):
    #         if attribute.startswith("join"):
    #             joiners.add(attribute)
    #     return joiners


    '''Change this method to a non fixture that just combines two types of datasets'''
    @pytest.fixture(scope="class")
    def all_omics_combos(self, get_public_dataset_objects):
        omics_lists_by_cancer = dict()
        
        # iterate through cancer objects
        cancer_sets = get_public_dataset_objects[0]
        for (cancer_type, cancer_object) in cancer_sets.items():
            omics_list = list() # list of omics datasets for current cancer
            # iterate through datasets of each cancer
            for (dataset, dimensions) in cancer_object.get_data_list().items():
                # check for 'omics' datasets and add to omics list for current cancer
                if dataset.__contains__("omics"):
                    omics_list.append((dataset, dimensions)) #including dimensions
            # make list of omics pairings to be used in join tests
            omics_combos = [(a, b) for a, b in itertools.combinations(omics_list, 2)]
            omics_lists_by_cancer[cancer_object] = omics_combos

        return omics_lists_by_cancer

    def test_join_omics_to_omics(self, all_omics_combos):
        for (cancer_object, omics_combos) in all_omics_combos.items():
            for (omics_1, omics_2) in omics_combos:
                dimensions = 1
                dataset = 0
                expected_columns = omics_1[dimensions]["columns"] + omics_2[dimensions]["columns"]
                joiner = getattr(cancer_object, "join_omics_to_omics")
                df = joiner(omics_1[dataset], omics_2[dataset])
                assert df.shape[1] == expected_columns

    def test_join_omics_to_mutations():
        pass

    def test_join_metadata_to_metadata():
        pass

    def test_join_metadata_to_omics():
        pass

    def test_join_metadata_to_mutations():
        pass

    def test_multi_join():
        pass
