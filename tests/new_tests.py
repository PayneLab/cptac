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
from collections import namedtuple

NO_INTERNET = True

class TestGetters:

    def test_all_getters(self):
        """Test all getters for a dataset."""

        getter_names = [name for name in dir(cptac.dataset.Dataset)
            if name.startswith("get_")
            and name not in ("get_cancer_type", "get_genotype_all_vars")]

        dss = self._get_dataset_tuples()

        for ds in dss:
            self._check_single_getter(ds, getter_names)

    # Helper functions
    def _check_single_getter(self, ds_tuple, getter_names):
        """Test a single getter from a dataset."""

        # Instantiate the dataset
        ds = ds_tuple.function(ds_tuple.version)

        for getter_name in getter_names:

            # Get the getter function
            getter = getattr(ds, getter_name)

            # Call the getter to get the dataframe
            try:
                df = getter()

            except cptac.exceptions.DataFrameNotIncludedError:
                warnings.warn(f"The {getter_name[4:]} dataframe was not found in the {ds_tuple.name} dataset.")

            else:
                # Check index and column names
                assert df.index.name == "Patient_ID"
                assert df.columns.name == "Name"

                # Check no duplicate indices/headers
                if getter_name[4:] not in ["somatic_mutation", "treatment", "medical_history", "gene_fusion", "followup"]:
                    assert df.index.duplicated().sum() == 0
                    assert df.columns.duplicated().sum() == 0
                    assert df.duplicated().sum() == 0

                # Check no null sample statuses
                if getter_name == "get_clinical":
                    assert df.Sample_Tumor_Normal.isna().sum() == 0

    def _get_dataset_tuples(self):
        """Load the datasets and return a list of them."""

        DatasetTuple = namedtuple("DatasetTuple", ["name", "version", "function"])

        dss = [
            DatasetTuple("brca", "3.1", cptac.Brca),
            DatasetTuple("brca", "3.1.1", cptac.Brca),

            DatasetTuple("ccrcc", "0.0", cptac.Ccrcc),
            DatasetTuple("ccrcc", "0.1", cptac.Ccrcc),

            DatasetTuple("ccrcc", "0.1.1", cptac.Ccrcc),
            DatasetTuple("colon", "0.0", cptac.Colon),
            DatasetTuple("colon", "0.0.1", cptac.Colon),

            DatasetTuple("endometrial", "2.1", cptac.Endometrial),
            DatasetTuple("endometrial", "2.1.1", cptac.Endometrial),

            DatasetTuple("gbm", "1.0", cptac.Gbm),
            DatasetTuple("gbm", "2.0", cptac.Gbm),
            DatasetTuple("gbm", "2.1", cptac.Gbm),
            DatasetTuple("gbm", "3.0", cptac.Gbm),

            DatasetTuple("hnscc", "0.1", cptac.Hnscc),
            DatasetTuple("hnscc", "2.0", cptac.Hnscc),

            DatasetTuple("lscc", "1.0", cptac.Lscc),

            DatasetTuple("luad", "2.0", cptac.Luad),
            DatasetTuple("luad", "3.1", cptac.Luad),
            DatasetTuple("luad", "3.1.1", cptac.Luad),

            DatasetTuple("ovarian", "0.0", cptac.Ovarian),
            DatasetTuple("ovarian", "0.0.1", cptac.Ovarian),
        ]

        return dss
