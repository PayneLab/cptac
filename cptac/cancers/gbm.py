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

from cptac.cancers.cancer import Cancer

from cptac.cancers.awg.awggbm import AwgGbm
from cptac.cancers.awgconf.awgconfgbm import AwgConfGbm
from cptac.cancers.bcm.bcmgbm import BcmGbm
from cptac.cancers.broad.broadgbm import BroadGbm
from cptac.cancers.pdc.pdcgbm import PdcGbm
from cptac.cancers.umich.umichgbm import UmichGbm
from cptac.cancers.washu.washugbm import WashuGbm
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized


class Gbm(Cancer):

    def __init__(self, no_internet=False):
        """Load all the data sources with GBM data and provide an interface to them."""

        super().__init__(cancer_type="gbm")

        self._sources["awg"] = AwgGbm(version="latest", no_internet=no_internet)
        self._sources["awgconf"] = AwgConfGbm(version="latest", no_internet=no_internet)
        self._sources["bcm"] = BcmGbm(version="latest", no_internet=no_internet)
        self._sources["broad"] = BroadGbm(version="latest", no_internet=no_internet)
        self._sources["mssm"] = Mssm(filter_type='gbm', version="latest", no_internet=no_internet)
        self._sources["pdc"] = PdcGbm(version="latest", no_internet=no_internet)
        self._sources["umich"] = UmichGbm(version="latest", no_internet=no_internet)
        self._sources["washu"] = WashuGbm(version="latest", no_internet=no_internet)
        self._sources["harmonized"] = Harmonized(filter_type='gbm', version="latest", no_internet=no_internet)
