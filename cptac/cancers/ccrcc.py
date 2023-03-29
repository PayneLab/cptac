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

from cptac.cancers.awg.awgccrcc import AwgCcrcc
from cptac.cancers.bcm.bcmccrcc import BcmCcrcc
from cptac.cancers.broad.broadccrcc import BroadCcrcc
from cptac.cancers.pdc.pdcccrcc import PdcCcrcc
from cptac.cancers.umich.umichccrcc import UmichCcrcc
from cptac.cancers.washu.washuccrcc import WashuCcrcc
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized


class Ccrcc(Cancer):

    def __init__(self, no_internet=False):
        """Load all the data sources with ccRCC data and provide an interface to them."""

        super().__init__(cancer_type="ccrcc")
        
        self._sources["awg"] = AwgCcrcc(version="latest", no_internet=no_internet)
        self._sources["bcm"] = BcmCcrcc(version="latest", no_internet=no_internet)
        self._sources["broad"] = BroadCcrcc(version="latest", no_internet=no_internet)
        self._sources["mssm"] = Mssm(filter_type='ccrcc', version="latest", no_internet=no_internet)
        self._sources["pdc"] = PdcCcrcc(version="latest", no_internet=no_internet)
        self._sources["umich"] = UmichCcrcc(version="latest", no_internet=no_internet)
        self._sources["washu"] = WashuCcrcc(version="latest", no_internet=no_internet)
        self._sources["harmonized"] = Harmonized(filter_type='ccrcc', version="latest", no_internet=no_internet)
