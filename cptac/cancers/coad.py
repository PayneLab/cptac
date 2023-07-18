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

# Imports from the cptac.cancers package
from cptac.cancers.cancer import Cancer

# Imports from the cptac.cancers sub-packages
from cptac.cancers.bcm.bcmcoad import BcmCoad
from cptac.cancers.broad.broadcoad import BroadCoad
from cptac.cancers.umich.umichcoad import UmichCoad
from cptac.cancers.washu.washucoad import WashuCoad
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized


class Coad(Cancer):
    """
    The Coad class is a subclass of the Cancer class and provides an interface to
    various data sources with COAD (Colon adenocarcinoma) data.

    Attributes:
        _sources (dict): A dictionary that holds instances of data source objects.
    """

    def __init__(self, no_internet=False):
        """
        Constructor for the Coad class. It initializes instances of various data source
        classes and adds them to the _sources attribute.

        Parameters:
            no_internet (bool): A flag indicating whether to download data or not. If True,
                                the data source objects are initialized without downloading
                                data. Defaults to False.
        """
        super().__init__(cancer_type="coad")

        # Initialize data sources and add them to the _sources dictionary
        self._sources["bcm"] = BcmCoad(no_internet=no_internet)
        self._sources["broad"] = BroadCoad(no_internet=no_internet)
        self._sources["mssm"] = Mssm(filter_type='coad', no_internet=no_internet)
        self._sources["umich"] = UmichCoad(no_internet=no_internet)
        self._sources["washu"] = WashuCoad(no_internet=no_internet)
        self._sources["harmonized"] = Harmonized(filter_type='coad', no_internet=no_internet)

