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

    def _combinations(self, list1):
       combo_list = [ (a, b) for a, b in itertools.combinations(list1, 2) ]
       return combo_list

    # verify that each join produces a dataframe with the expected number of rows and columns
    def test_join_omics_to_omics(self, get_cancer_test_units):
        # loop through cancers
        for cancer in get_cancer_test_units:
            # generate omics combos per cancer and iterate
            omics_combos = self._combinations(cancer.omics)
            for o in omics_combos:
                # test the join
                omics1 = o[0]
                omics2 = o[1]
                expected_size = (omics1.shape[0] + omics2.shape[0], omics1.shape[1] + omics2.shape[1])
                df = cancer.join_omics(o[0], o[2])
                assert df.shape == expected_size

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
