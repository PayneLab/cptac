#   Copyright 2018 Samuel Payne sam_payne@byu.edu
#   Licensed under the Apache License, Version 2.0 (the "License");
#   You may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from cptac.cancers.cancer import Cancer
from cptac.cancers.bcm.bcmucec import BcmUcec
from cptac.cancers.broad.broaducec import BroadUcec
from cptac.cancers.umich.umichucec import UmichUcec
from cptac.cancers.washu.washuucec import WashuUcec
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized

class Ucec(Cancer):
    """
    The Ucec class is a child class of the base Cancer class, intended for handling 
    Uterine Corpus Endometrial Carcinoma (UCEC) data from different sources.

    Attributes:
        _sources (dict): A dictionary of source objects that load and hold data 
                         pertaining to UCEC cancer.
    """

    def __init__(self, no_internet=False):
        """
        Initializes the Ucec object, loading the relevant data sources for UCEC data. 

        Parameters:
            no_internet (bool): If set to True, the data source objects are initialized 
                                without downloading data. Default is False.
        """
        super().__init__(cancer_type="ucec")

        # Initialize data source objects for various sources of UCEC data
        self._sources["bcm"] = BcmUcec(no_internet=no_internet)
        self._sources["broad"] = BroadUcec(no_internet=no_internet)
        self._sources["mssm"] = Mssm(filter_type='ucec', no_internet=no_internet)
        self._sources["umich"] = UmichUcec(no_internet=no_internet)
        self._sources["washu"] = WashuUcec(no_internet=no_internet)
        self._sources["harmonized"] = Harmonized(filter_type='ucec', no_internet=no_internet)

