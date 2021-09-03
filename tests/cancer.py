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


#   The purpose of this class is to organize a cancer object's datasets by
#   type. dataset.py in the cptac package defines a lot of methods and members
#   but there is no built-in way to call them in batches by type for testing.

from _typeshed import Self


class Cancer:

    metadata_types = [
            'clinical',
            'derived_molecular',
            'experimental_design',
            # See dataset.py for why these aren't included:
                #'medical_history',
                #'treatment',
                #'followup'
        ]

    def __init__(self, cancer_type, cancer_object):
        """
        Initialize a Cancer object.

        Cancer class is used as a wrapper for cptac.[Cancer] objects that will be tested.

        Parameters:
        cancer_type (string): name of the cancer
        cancer_object (cptac.[Cancer]): Instance of the cptac.[Cancer] class

        """
        self.cancer_type = cancer_type
        self.cancer_object = cancer_object

        self.metadata = list()
        self.omics = list()
        # self.mutations = list()

        self.valid_getters = dict()
        self.invalid_getters = dict()

        self._sort_datasets()
        self._sort_getters()


    def _sort_datasets(self):
        # categorize datasets for join tests
        datasets = self.cancer_object.get_data_list().items()
        for (dataset, dimensions) in datasets:
            if dataset in Cancer.metadata_types:
                self.metadata.append(dataset)
            elif dataset.__contains__('omics'):
                self.omics.append(dataset)

    def _sort_getters(self):
        # collect all possible getters
        all_getters = set()
        for attribute in dir(self.cancer_object):
            if attribute.startswith("get_"):
                all_getters.add(attribute)

        ### sift valid and invalid getters
        datasets = self.cancer_object.get_data_list().keys()

        # valid getters
        for d in datasets:
            getter_name = "get_" + d
            valid_getter = getattr(self.cancer_object, getter_name)
            self.valid_getters[getter_name] = valid_getter

        # invalid getters
        for getter in all_getters:
            if getter_name not in self.valid_getters.keys():
                if getter_name.startswith("CNV") and self.cancer_type == "Ucecconf":
                    pass
                else:
                    g = getattr(self.cancer_object, getter_name)
                    self.invalid_getters[getter_name] = g

    
