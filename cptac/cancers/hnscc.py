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

from cptac.cancers.awg.awghnscc import AwgHnscc
from cptac.cancers.bcm.bcmhnscc import BcmHnscc
from cptac.cancers.broad.broadhnscc import BroadHnscc
from cptac.cancers.pdc.pdchnscc import PdcHnscc
from cptac.cancers.umich.umichhnscc import UmichHnscc
from cptac.cancers.washu.washuhnscc import WashuHnscc
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized


class Hnscc(Cancer):

    def __init__(self, no_internet=False):
        """Load all the data sources with HNSCC data and provide an interface to them."""

        super().__init__(cancer_type="hnscc")
        
        self._sources["awg"] = AwgHnscc(version="latest", no_internet=no_internet)
        self._sources["bcm"] = BcmHnscc(version="latest", no_internet=no_internet)
        self._sources["broad"] = BroadHnscc(version="latest", no_internet=no_internet)
        self._sources["mssm"] = Mssm(filter_type='hnscc', version="latest", no_internet=no_internet)
        self._sources["pdc"] = PdcHnscc(version="latest", no_internet=no_internet)
        self._sources["umich"] = UmichHnscc(version="latest", no_internet=no_internet)
        self._sources["washu"] = WashuHnscc(version="latest", no_internet=no_internet)
        self._sources["harmonized"] = Harmonized(filter_type='hnscc', version="latest", no_internet=no_internet)
