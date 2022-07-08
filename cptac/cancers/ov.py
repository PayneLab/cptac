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

from cptac.cancers.awg.awgov import AwgOv
from cptac.cancers.bcm.bcmov import BcmOv
from cptac.cancers.broad.broadov import BroadOv
from cptac.cancers.pdc.pdcov import PdcOv
from cptac.cancers.umich.umichov import UmichOv
from cptac.cancers.washu.washuov import WashuOv
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized


class Ov(Cancer):

    def __init__(self, no_internet=False):
        """Load all the data sources with OV data and provide an interface to them."""

        super().__init__(cancer_type="ov")
        
        self._sources["awg"] = AwgOv(version="latest", no_internet=no_internet)
        self._sources["bcm"] = BcmOv(version="latest", no_internet=no_internet)
        self._sources["broad"] = BroadOv(version="latest", no_internet=no_internet)
        self._sources["mssm"] = Mssm(filter_type='ov', version="latest", no_internet=no_internet)
        self._sources["pdc"] = PdcOv(version="latest", no_internet=no_internet)
        self._sources["umich"] = UmichOv(version="latest", no_internet=no_internet)
        self._sources["washu"] = WashuOv(version="latest", no_internet=no_internet)
        self._sources["harmonized"] = Harmonized(filter_type='ov', version="latest", no_internet=no_internet)
