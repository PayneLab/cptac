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
        # TODO: add docstring here
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

        # sift valid and invalid getters
        datasets = self.cancer_object.get_data_list().items()
        for (dataset, dimensions) in datasets:
            getter_name = "get_" + dataset
            if getter_name not in all_getters:
                g = getattr(self.cancer_object, getter_name)
                self.invalid_getters[getter_name] = g

            elif getter_name in all_getters:
                g = getattr(self.cancer_object, getter_name)
                self.valid_getters[getter_name] = g

            else:
                raise Exception(f"{getter_name} unable to be sorted.")

    
