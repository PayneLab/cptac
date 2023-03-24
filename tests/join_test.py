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
import itertools as it

from .conftest import get_cancer_inputs

metadata_types = { 'clinical', 'derived_molecular', 'experimental_design'}
    # See dataset.py for why these aren't included: 'medical_history','treatment','followup'

omics_types = { 'acetylproteomics', 'circular_RNA', 'CNV', 'lincRNA',
        'lipidomics', 'metabolomics', 'miRNA', 'phosphoproteomics',
        'phosphoproteomics_gene', 'proteomics', 'somatic_mutation_binary',
        'transcriptomics', 'CNV_log2ratio', 'CNV_gistic'
}

multi_join_types = { "acetylproteomics", "CNV", "CNV_gistic", "CNV_log2ratio",
    "phosphoproteomics", "phosphoproteomics_gene", "proteomics", 
    "somatic_mutation_binary", "somatic_mutation", "transcriptomics", 
    "clinical", "derived_molecular", "experimental_design"
}

def get_join_inputs():
    """ 
    Get nested dictionary of each cancer's valid datatypes, grouped by category."
    @return: Nested dict in form {cancer_name:{'omics':[...], 'metadata':[...], 'multi':[...]}
    """
    join_inputs = dict()
    for cancer, source, dtype in get_cancer_inputs():
        if cancer not in join_inputs:
            join_inputs[cancer] = {'omics':[], 'metadata':[], 'multi':[]}
        if dtype in omics_types:
            join_inputs[cancer]['omics'].append((source, dtype))
        if dtype in metadata_types:
            join_inputs[cancer]['metadata'].append((source, dtype))
        if dtype in multi_join_types:
            join_inputs[cancer]['omics'].append((source, dtype))

    return join_inputs


@pytest.mark.parametrize("cancer_name", get_join_inputs())
def test_join_omics_to_omics(cancer_name):
    print(cancer_name)
    _run_combos(cancer_name, 'omics')
    


def _run_combos(cancer_name, cat1, cat2=None):
    cancer_obj = getattr(cptac, cancer_name.title())()
    cancer_combos = get_join_inputs()[cancer_name]
    if cat2 is None:
        test_units = it.combinations(cancer_combos[cat1], 2)
        join_function_name = f"join_{cat1}_to_{cat1}"
    else:
        test_units = it.product(cancer_combos[cat1], cancer_combos[cat2], 2)
        join_function_name = f"join_{cat1}_to_{cat2}"
    join_function = getattr(cancer_obj, join_function_name)

    for source1, source2 in test_units:
        #TEMP? or Keep exclusion of awg data?
        if source1[0] == 'awg' or source2[0] == 'awg':
            continue
        if source1[0] == 'broad' or source2[0] == 'broad':
            print(f"   one of the datasets in ({source1}, {source2}) is too large")
            continue
        print(f"   testing join of {source1} to {source2}")
        try:
            joined_data = join_function(
                    df1_name = source1[1],
                    df2_name = source2[1],
                    df1_source = source1[0],
                    df2_source = source2[0]
            )
            #Include more robust testing?
            # df1 = cancer_obj.get_dataframe(source1[1], source1[0])
            # df2 = cancer_obj.get_dataframe(source2[1], source2[0])
            #assert df1.shape[1] + df2.shape[1] == joined_data.shape[1]
        except (Exception) as e:
            print("      test failed!")
            pytest.fail(f"Unable to perform join on {source1[1]} and {source2[1]} for {cancer_name}.\n{e}")



#def test_join_omics_to_mutations(self, get_cancer_test_units):
#    # loop through cancers
#    for cancer in get_cancer_test_units:
#        # generate omics-mutation genes combos
#        combos = self._combinations(cancer.omics, [cancer.mutation_genes])
#        # test each combo
#        self._run_combos(cancer, combos, cancer.cancer_object.join_omics_to_mutations)
#
#def test_join_metadata_to_metadata(self, get_cancer_test_units):
#    # loop through cancers
#    for cancer in get_cancer_test_units:
#        # generate metadata combos per cancer
#        combos = self._combinations(cancer.metadata)
#        # test each combo
#        self._run_combos(cancer, combos, cancer.cancer_object.join_metadata_to_metadata)
#
#def test_join_metadata_to_omics(self, get_cancer_test_units):
#    # loop through cancers
#    for cancer in get_cancer_test_units:
#        # generate metadata combos per cancer
#        combos = self._combinations(cancer.metadata, cancer.omics)
#        # test each combo
#        self._run_combos(cancer, combos, cancer.cancer_object.join_metadata_to_omics)
#
#def test_join_metadata_to_mutations(self, get_cancer_test_units):
#    # loop through cancers
#    for cancer in get_cancer_test_units:
#        # generate metadata-mutation genes combos
#        combos = self._combinations(cancer.metadata, [cancer.mutation_genes])
#        # test each combo
#        self._run_combos(cancer, combos, cancer.cancer_object.join_metadata_to_mutations)
#
#def test_multi_join(self, get_cancer_test_units):
#    for cancer in get_cancer_test_units:
#        try:
#            df = cancer.cancer_object.multi_join(cancer.multi_joinables)
#        except:
#            pytest.fail(f"Unable to perform multijoin on {cancer.multi_joinables} for {cancer.cancer_type}")
#
#
#def _combinations(self, list1, list2=None):
#    '''
#    Returns:
#        unique list of combo tuples
#    '''
#    if not list2:
#        combo_list = [ (a, b) for a, b in itertools.combinations(list1, 2) ] 
#    else:
#        combo_list = [ (a, b) for a, b in itertools.product(list1, list2) ]
#    
#    return combo_list
#
#def _run_combos(self, cancer, combos, func):
#    '''
#    Test each combo
#
#    Args:
#        cancer: the cancer whose data will be tested
#        combos: the combos to be run
#        func: the function the combos will be tested on
#    '''
#    for c in combos:
#            # ds is for dataset
#            ds1 = c[0]
#            # ds1_df = cancer.get_dataset(ds1)
#            ds2 = c[1]
#            # ds2_df = cancer.get_dataset(ds2)
#            # expected_columns = ds1_df.shape[1] + ds2_df.shape[1]
#            try:
#                df = func(ds1, ds2)
#            except (Exception) as e:
#                pytest.fail(f"Unable to perform join on {ds1} and {ds2} for {cancer.cancer_type}.\n{e}")
#            # verify the join worked based on column counts
#            #assert df.shape[1] == expected_columns






### OLD STUFF ##
            
## TODO: Check use cases for standard usage and then try to mess that up
#class TestJoin:
#    '''
#    Tests to verify that each join produces a working dataframe
#    '''
#   
#    def test_join_omics_to_omics(self, get_cancer_test_units):
#        # loop through cancers
#        for cancer in get_cancer_test_units:
#            # generate omics combos per cancer
#            combos = self._combinations(cancer.omics)
#            # test each combo
#            self._run_combos(cancer, combos, cancer.cancer_object.join_omics_to_omics)
#
#    def test_join_omics_to_mutations(self, get_cancer_test_units):
#        # loop through cancers
#        for cancer in get_cancer_test_units:
#            # generate omics-mutation genes combos
#            combos = self._combinations(cancer.omics, [cancer.mutation_genes])
#            # test each combo
#            self._run_combos(cancer, combos, cancer.cancer_object.join_omics_to_mutations)
#
#    def test_join_metadata_to_metadata(self, get_cancer_test_units):
#        # loop through cancers
#        for cancer in get_cancer_test_units:
#            # generate metadata combos per cancer
#            combos = self._combinations(cancer.metadata)
#            # test each combo
#            self._run_combos(cancer, combos, cancer.cancer_object.join_metadata_to_metadata)
#
#    def test_join_metadata_to_omics(self, get_cancer_test_units):
#        # loop through cancers
#        for cancer in get_cancer_test_units:
#            # generate metadata combos per cancer
#            combos = self._combinations(cancer.metadata, cancer.omics)
#            # test each combo
#            self._run_combos(cancer, combos, cancer.cancer_object.join_metadata_to_omics)
#
#    def test_join_metadata_to_mutations(self, get_cancer_test_units):
#        # loop through cancers
#        for cancer in get_cancer_test_units:
#            # generate metadata-mutation genes combos
#            combos = self._combinations(cancer.metadata, [cancer.mutation_genes])
#            # test each combo
#            self._run_combos(cancer, combos, cancer.cancer_object.join_metadata_to_mutations)
#
#    def test_multi_join(self, get_cancer_test_units):
#        for cancer in get_cancer_test_units:
#            try:
#                df = cancer.cancer_object.multi_join(cancer.multi_joinables)
#            except:
#                pytest.fail(f"Unable to perform multijoin on {cancer.multi_joinables} for {cancer.cancer_type}")
#
#
#    def _combinations(self, list1, list2=None):
#        '''
#        Returns:
#            unique list of combo tuples
#        '''
#        if not list2:
#            combo_list = [ (a, b) for a, b in itertools.combinations(list1, 2) ] 
#        else:
#            combo_list = [ (a, b) for a, b in itertools.product(list1, list2) ]
#        
#        return combo_list
#
#    def _run_combos(self, cancer, combos, func):
#        '''
#        Test each combo
#
#        Args:
#            cancer: the cancer whose data will be tested
#            combos: the combos to be run
#            func: the function the combos will be tested on
#        '''
#        for c in combos:
#                # ds is for dataset
#                ds1 = c[0]
#                # ds1_df = cancer.get_dataset(ds1)
#                ds2 = c[1]
#                # ds2_df = cancer.get_dataset(ds2)
#                # expected_columns = ds1_df.shape[1] + ds2_df.shape[1]
#                try:
#                    df = func(ds1, ds2)
#                except (Exception) as e:
#                    pytest.fail(f"Unable to perform join on {ds1} and {ds2} for {cancer.cancer_type}.\n{e}")
#                # verify the join worked based on column counts
#                #assert df.shape[1] == expected_columns
#                
