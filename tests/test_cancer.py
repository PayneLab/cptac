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

import pytest

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

    valid_omics_dfs = [
            'acetylproteomics',
            'circular_RNA',
            'CNV',
            'lincRNA',
            'lipidomics',
            'metabolomics',
            'miRNA',
            'phosphoproteomics',
            'phosphoproteomics_gene',
            'proteomics',
            'somatic_mutation_binary',
            'transcriptomics', 
            'CNV_log2ratio',
            'CNV_gistic'
             ]

    important_mutation_genes = ["TP53", "KRAS", "ARID1A", "PTEN", "EGFR"]

    multi_join_types = [
        "acetylproteomics", 
        "CNV",
        "CNV_gistic",
        "CNV_log2ratio",
        "phosphoproteomics", 
        "phosphoproteomics_gene", 
        "proteomics", 
        "somatic_mutation_binary", 
        "somatic_mutation", 
        "transcriptomics", 
        "clinical", 
        "derived_molecular", 
        "experimental_design"
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
        self.multi_joinables = dict()

        self._sort_datasets()
        self._sort_getters()
        self._gather_mutation_genes()


    def _sort_datasets(self):
        # categorize datasets for join tests
        # omics, metadata, 
        
        datasets = self.cancer_object.get_data_list().items()
        for (dataset, dimensions) in datasets:
            if dataset in Cancer.metadata_types:
                self.metadata.append(dataset)
            elif dataset in Cancer.valid_omics_dfs:
                self.omics.append(dataset)
            if dataset in ["clinical", "transcriptomics", "proteomics"]:
                self.multi_joinables[dataset] = list()

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
            try: 
                if d.startswith("CNV") and self.cancer_type == "Ucecconf":
                    getter_name = "get_CNV"
                else:
                    getter_name = "get_" + d
                    valid_getter = getattr(self.cancer_object, getter_name)
                    self.valid_getters[getter_name] = valid_getter
            except:
                pytest.fail(f"unable to add get {d} attribute")

        # invalid getters
        for getter in all_getters:
            if getter_name not in self.valid_getters.keys():
                g = getattr(self.cancer_object, getter_name)
                self.invalid_getters[getter_name] = g

    def _gather_mutation_genes(self):
        self.mutation_genes = list()
        if "somatic_mutation" in self.cancer_object.get_data_list():
            recorded_genes = self.cancer_object.get_somatic_mutation()["Gene"].tolist()
            for g in self.important_mutation_genes:
               if g in recorded_genes:
                   self.mutation_genes.append(g)

    def get_dataset(self, dataset, CNV_type="log2ratio"):
        '''
        Args:
            dataset: the desired dataset
            CNV_type: if the desired dataset is CNV and the cancer type is Ucecconf,
                you can specify which version of the dataset is returned.

        Returns: 
            adataframe for the dataset desired
        '''
        if dataset == "CNV" and self.cancer_type == "Ucecconf":
            return self.valid_getters["get_CNV"](CNV_type)
        return self.valid_getters["get_" + dataset]()

    def get_omics(self):
        return self.omics

    def get_metadata(self):
        return self.metadata
    
    def get_mutation_genes(self):
        return self.mutation_genes

    
