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

    def _combinations(self, get_cancer_test_units, dict1, dict2):
       #combo_list = [ (a, b) for a, b in itertools.product(list1, list2) ]
       pass

    def test_join_omics_to_omics(self, get_cancer_test_units):
        for (cancer_object, omics_combos) in self._combinations(self._list_omics()).items():
            for (omics_1, omics_2) in omics_combos:
                dimensions = 1
                dataset = 0
                expected_columns = omics_1[dimensions]["columns"] + omics_2[dimensions]["columns"]
                joiner = getattr(cancer_object, "join_omics_to_omics")
                df = joiner(omics_1[dataset], omics_2[dataset])
                assert df.shape[1] == expected_columns

    def test_join_omics_to_mutations(self, get_cancer_test_units):
        pass

    def test_join_metadata_to_metadata(self, get_cancer_test_units):
        pass

    def test_join_metadata_to_omics(self, get_cancer_test_units):
        pass

    def test_join_metadata_to_mutations(self, get_cancer_test_units):
        pass

    def test_multi_join(self, get_cancer_test_units):
        pass
