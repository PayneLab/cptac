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
import os
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
        """Parses index files to generate a list of named tuples of all exisiting dataset names and version numbers, and references to their loader functions."""

        # Download latest versions of all datasets, so we have the latest versions 
        # of all the index files and can get complete lists of all versions
        cptac.download(dataset="all", version="latest", redownload=False)

        # Construct the path to the cptac directory
        path_here = os.path.abspath(os.path.dirname(__file__))
        dirs = path_here.split(os.sep)
        dirs = dirs[:-1]
        dirs.append("cptac")
        cptac_path = os.sep.join(dirs)
        
        # Get a list of all data directories in the cptac directory, and
        # turn it into a list of the paths to all the index files
        all_dirs = [f.path for f in os.scandir(cptac_path) if f.is_dir()]
        data_dirs = [dir for dir in all_dirs if dir.split(os.sep)[-1].startswith("data_")]
        index_paths = [os.sep.join([dir, "index.txt"]) for dir in data_dirs]

        # Read in the index file from each data directory, and parse out all version numbers
        versions = {}
        for path in index_paths:

            # Get the dataset name
            dataset = path.split(os.sep)[-2].split("_", maxsplit=1)[1]

            # Read in the index file and parse out the version numbers
            with open(path, "r") as fp:
                versions_list = [line.strip()[1:] for line in fp.readlines() if line.startswith("#")]

            # Save the versions for that dataset
            versions[dataset] = versions_list

        # Create a named tuple for each version of each dataset
        DatasetTuple = namedtuple("DatasetTuple", ["name", "version", "function"])
        dataset_tuples = []

        for dataset in versions.keys():
            data_versions = versions[dataset]
            for version in data_versions:
                version_tuple = DatasetTuple(name=dataset, version=version, function=dataset)
                dataset_tuples.append(version_tuple)

        # TODO: Get functions. Can we scrape them from a namespace by checking inheritance from Dataset class?
        # Other option is to further standardize naming, but that would be more susceptible to bugs. If you do,
        # leave more troubleshooting notes detailing that issues might arise from incorrectly named datasets.

        # Also, move the downloading step out of this function, and into the class init equivalent. Have that also call
        # this, and download all previous data versions.

        return dataset_tuples
