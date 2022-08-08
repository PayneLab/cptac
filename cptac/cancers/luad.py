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

from cptac.cancers.awg.awgluad import AwgLuad
from cptac.cancers.bcm.bcmluad import BcmLuad
from cptac.cancers.broad.broadluad import BroadLuad
from cptac.cancers.pdc.pdcluad import PdcLuad
from cptac.cancers.umich.umichluad import UmichLuad
from cptac.cancers.washu.washuluad import WashuLuad
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized


class Luad(Cancer):

    def __init__(self, no_internet=False):
        """Load all the data sources with LUAD data and provide an interface to them."""

        super().__init__(cancer_type="luad")
        
        self._sources["awg"] = AwgLuad(version="latest", no_internet=no_internet)
        self._sources["bcm"] = BcmLuad(version="latest", no_internet=no_internet)
        self._sources["broad"] = BroadLuad(version="latest", no_internet=no_internet)
        self._sources["mssm"] = Mssm(filter_type='luad', version="latest", no_internet=no_internet)
        self._sources["pdc"] = PdcLuad(version="latest", no_internet=no_internet)
        self._sources["umich"] = UmichLuad(version="latest", no_internet=no_internet)
        self._sources["washu"] = WashuLuad(version="latest", no_internet=no_internet)
        self._sources["harmonized"] = Harmonized(filter_type='luad', version="latest", no_internet=no_internet)
