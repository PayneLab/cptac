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
import cptac.exceptions as ex
import pandas as pd

class PancanDataset:

    def __init__(self, cancer_type, versions, no_internet):
        """Initialize variables for a PancanDataset object."""
        # Assumed instance variables:
        # self._datasets = {}

    def get_clinical(self, source, tissue_type="both"):
        """Get the clinical dataframe from the specified data source."""
        return self._get_dataframe("clinical", source, tissue_type)

    def get_proteomics(self, source, tissue_type="both"):
        """Get the proteomics dataframe from the specified data source."""
        return self._get_dataframe("clinical", source, tissue_type)

    # "Private" methods
    def _get_dataframe(self, name, source, tissue_type="both"):
        """Check that a given dataframe from a given source exists, and return a copy if it does."""

        if source in self._datasets.keys():
            return self._datasets[source]._get_dataframe(name, tissue_type)
        else:
            raise DataTypeNotInSourceError(f"{name} dataframe not available from source {source}.")
