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

from cptac.cancers.awg.awgbrca import AwgBrca
from cptac.cancers.bcm.bcmbrca import BcmBrca
from cptac.cancers.broad.broadbrca import BroadBrca
from cptac.cancers.pdc.pdcbrca import PdcBrca
from cptac.cancers.umich.umichbrca import UmichBrca
from cptac.cancers.washu.washubrca import WashuBrca
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized


class Brca(Cancer):

    def __init__(self, no_internet=False):
        """Load all the data sources with BRCA data and provide an interface to them."""

        super().__init__(cancer_type="brca")

        self._sources["awg"] = AwgBrca(version="latest", no_internet=no_internet)
        self._sources["bcm"] = BcmBrca(version="latest", no_internet=no_internet)
        self._sources["broad"] = BroadBrca(version="latest", no_internet=no_internet)
        self._sources["mssm"] = Mssm(filter_type='brca', version="latest", no_internet=no_internet)
        self._sources["pdc"] = PdcBrca(version="latest", no_internet=no_internet)
        self._sources["umich"] = UmichBrca(version="latest", no_internet=no_internet)
        self._sources["washu"] = WashuBrca(version="latest", no_internet=no_internet)
        self._sources["harmonized"] = Harmonized(filter_type='brca', version="latest", no_internet=no_internet)
