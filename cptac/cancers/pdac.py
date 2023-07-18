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

# Importing the base Cancer class
from cptac.cancers.cancer import Cancer

# Importing data sources for PDAC (Pancreatic Ductal Adenocarcinoma)
from cptac.cancers.bcm.bcmpdac import BcmPdac
from cptac.cancers.broad.broadpdac import BroadPdac
from cptac.cancers.umich.umichpdac import UmichPdac
from cptac.cancers.washu.washupdac import WashuPdac
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized

class Pdac(Cancer):
    """
    The Pdac class inherits from the Cancer class and provides access to various data sources 
    for PDAC (Pancreatic Ductal Adenocarcinoma).

    Attributes:
        _sources (dict): A dictionary that stores instances of data source objects.
    """

    def __init__(self, no_internet=False):
        """
        The constructor for the Pdac class. It initializes instances of various data source 
        classes and adds them to the _sources attribute.

        Parameters:
            no_internet (bool): If set to True, the data source objects are initialized 
                                without downloading data. Default is False.
        """
        super().__init__(cancer_type="pdac")

        # Initialize data sources and add them to the _sources dictionary
        self._sources["bcm"] = BcmPdac(no_internet=no_internet)
        self._sources["broad"] = BroadPdac(no_internet=no_internet)
        self._sources["mssm"] = Mssm(filter_type='pdac', no_internet=no_internet)
        self._sources["umich"] = UmichPdac(no_internet=no_internet)
        self._sources["washu"] = WashuPdac(no_internet=no_internet)
        self._sources["harmonized"] = Harmonized(filter_type='pdac', no_internet=no_internet)

