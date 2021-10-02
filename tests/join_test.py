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

# TODO: Check use cases for standard usage and then try to mess that up
class TestJoin:
    '''
    Tests to verify that each join produces a working dataframe
    '''
   
    def test_join_omics_to_omics(self, get_cancer_test_units):
        # loop through cancers
        for cancer in get_cancer_test_units:
            # generate omics combos per cancer
            combos = self._combinations(cancer.omics)
            # test each combo
            self._run_combos(cancer, combos, cancer.cancer_object.join_omics_to_omics)

    def test_join_omics_to_mutations(self, get_cancer_test_units):
        # loop through cancers
        for cancer in get_cancer_test_units:
            # generate omics-mutation genes combos
            combos = self._combinations(cancer.omics, [cancer.mutation_genes])
            # test each combo
            self._run_combos(cancer, combos, cancer.cancer_object.join_omics_to_mutations)

    def test_join_metadata_to_metadata(self, get_cancer_test_units):
        # loop through cancers
        for cancer in get_cancer_test_units:
            # generate metadata combos per cancer
            combos = self._combinations(cancer.metadata)
            # test each combo
            self._run_combos(cancer, combos, cancer.cancer_object.join_metadata_to_metadata)

    def test_join_metadata_to_omics(self, get_cancer_test_units):
        # loop through cancers
        for cancer in get_cancer_test_units:
            # generate metadata combos per cancer
            combos = self._combinations(cancer.metadata, cancer.omics)
            # test each combo
            self._run_combos(cancer, combos, cancer.cancer_object.join_metadata_to_omics)

    def test_join_metadata_to_mutations(self, get_cancer_test_units):
        # loop through cancers
        for cancer in get_cancer_test_units:
            # generate metadata-mutation genes combos
            combos = self._combinations(cancer.metadata, [cancer.mutation_genes])
            # test each combo
            self._run_combos(cancer, combos, cancer.cancer_object.join_metadata_to_mutations)

    def test_multi_join(self, get_cancer_test_units):
        for cancer in get_cancer_test_units:
            try:
                df = cancer.cancer_object.multi_join(cancer.multi_joinables)
            except:
                pytest.fail(f"Unable to perform multijoin on {cancer.multi_joinables} for {cancer.cancer_type}")


    def _combinations(self, list1, list2=None):
        '''
        Returns:
            unique list of combo tuples
        '''
        if not list2:
            combo_list = [ (a, b) for a, b in itertools.combinations(list1, 2) ] 
        else:
            combo_list = [ (a, b) for a, b in itertools.product(list1, list2) ]
        
        return combo_list

    def _run_combos(self, cancer, combos, func):
        '''
        Test each combo

        Args:
            cancer: the cancer whose data will be tested
            combos: the combos to be run
            func: the function the combos will be tested on
        '''
        for c in combos:
                # ds is for dataset
                ds1 = c[0]
                # ds1_df = cancer.get_dataset(ds1)
                ds2 = c[1]
                # ds2_df = cancer.get_dataset(ds2)
                # expected_columns = ds1_df.shape[1] + ds2_df.shape[1]
                try:
                    df = func(ds1, ds2)
                except (Exception) as e:
                    pytest.fail(f"Unable to perform join on {ds1} and {ds2} for {cancer.cancer_type}.\n{e}")
                # verify the join worked based on column counts
                #assert df.shape[1] == expected_columns
                
