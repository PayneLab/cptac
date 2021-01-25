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

        self._cancer_type = cancer_type
        self._datasets = {} # Child class __init__ needs to fill this

    # Data getters
    def get_clinical(self, source, tissue_type="both", imputed=False):
        """Get the clinical dataframe from the specified data source."""
        return self._get_dataframe("clinical", source, tissue_type, imputed=imputed)
    
    def get_demographic(self, source, tissue_type="both", imputed=False):
        """Get the demographic dataframe from the specified data source."""
        return self._get_dataframe("demographic", source, tissue_type, imputed=imputed)
    
    def get_general_medical(self, source, tissue_type="both", imputed=False):
        """Get the general_medical dataframe from the specified data source."""
        return self._get_dataframe("general_medical", source, tissue_type, imputed=imputed)
    
    def get_cancer_diagnosis(self, source, tissue_type="both", imputed=False):
        """Get the cancer_diagnosis dataframe from the specified data source."""
        return self._get_dataframe("cancer_diagnosis", source, tissue_type, imputed=imputed)
    
    def get_followup(self, source, tissue_type="both", imputed=False):
        """Get the followup dataframe from the specified data source."""
        return self._get_dataframe("followup", source, tissue_type, imputed=imputed)

    def get_proteomics(self, source, tissue_type="both", imputed=False):
        """Get the proteomics dataframe from the specified data source."""
        return self._get_dataframe("proteomics", source, tissue_type, imputed=imputed)

    def get_transcriptomics(self, source, tissue_type="both", imputed=False):
        """Get the proteomics dataframe from the specified data source."""
        return self._get_dataframe("transcriptomics", source, tissue_type, imputed=imputed)

    # Help functions
    def get_cancer_type(self):
        return self._cancer_type

    # "Private" methods
    def _get_dataframe(self, name, source, tissue_type, imputed):
        """Check that a given dataframe from a given source exists, and return a copy if it does."""

        if imputed:
            name = name + "_imputed"

        if source in self._datasets.keys():
            return self._datasets[source]._get_dataframe(name, tissue_type)
        else:
            raise ex.DataSourceNotFoundError(f"Data source {source} not found for the {self._cancer_type} dataset.")
