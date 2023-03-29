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

from cptac.cancers.awg.awgcoad import AwgCoad
from cptac.cancers.bcm.bcmcoad import BcmCoad
from cptac.cancers.broad.broadcoad import BroadCoad
from cptac.cancers.pdc.pdccoad import PdcCoad
from cptac.cancers.umich.umichcoad import UmichCoad
from cptac.cancers.washu.washucoad import WashuCoad
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized


class Coad(Cancer):

    def __init__(self, no_internet=False):
        """Load all the data sources with COAD data and provide an interface to them."""

        super().__init__(cancer_type="coad")

        self._sources["awg"] = AwgCoad(version="latest", no_internet=no_internet)
        self._sources["bcm"] = BcmCoad(version="latest", no_internet=no_internet)
        self._sources["broad"] = BroadCoad(version="latest", no_internet=no_internet)
        self._sources["mssm"] = Mssm(filter_type='coad', version="latest", no_internet=no_internet)
        self._sources["pdc"] = PdcCoad(version="latest", no_internet=no_internet)
        self._sources["umich"] = UmichCoad(version="latest", no_internet=no_internet)
        self._sources["washu"] = WashuCoad(version="latest", no_internet=no_internet)
        self._sources["harmonized"] = Harmonized(filter_type='coad', version="latest", no_internet=no_internet)
