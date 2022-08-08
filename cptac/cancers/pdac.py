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

from cptac.cancers.awg.awgpdac import AwgPdac
from cptac.cancers.bcm.bcmpdac import BcmPdac
from cptac.cancers.broad.broadpdac import BroadPdac
from cptac.cancers.pdc.pdcpdac import PdcPdac
from cptac.cancers.umich.umichpdac import UmichPdac
from cptac.cancers.washu.washupdac import WashuPdac
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized


class Pdac(Cancer):

    def __init__(self, no_internet=False):
        """Load all the data sources with PDAC data and provide an interface to them."""

        super().__init__(cancer_type="pdac")

        self._sources["awg"] = AwgPdac(version="latest", no_internet=no_internet)
        self._sources["bcm"] = BcmPdac(version="latest", no_internet=no_internet)
        self._sources["broad"] = BroadPdac(version="latest", no_internet=no_internet)
        self._sources["mssm"] = Mssm(filter_type='pdac', version="latest", no_internet=no_internet)
        self._sources["pdc"] = PdcPdac(version="latest", no_internet=no_internet)
        self._sources["umich"] = UmichPdac(version="latest", no_internet=no_internet)
        self._sources["washu"] = WashuPdac(version="latest", no_internet=no_internet)
        self._sources["harmonized"] = Harmonized(filter_type='pdac', version="latest", no_internet=no_internet)
