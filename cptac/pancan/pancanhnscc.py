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

from .pancandataset import PancanDataset

from .bcmhnscc import BcmHnscc
from .mssmclinical import MssmClinical
from .pdchnscc import PdcHnscc
from .umichhnscc import UmichHnscc
from .washuhnscc import WashuHnscc
from .broadhnscc import BroadHnscc
from .harmonized import Harmonized
from .joiningdataset import JoiningDataset

#List sources to be downloaded
SOURCES = [
    "bcmhnscc",
    "mssmclinical",
    "pdchnscc",
    "umichhnscc",
    "washuhnscc",
    "broadhnscc",
    "harmonized"
]

class PancanHnscc(PancanDataset):

    def __init__(self, version="latest", no_internet=False):
        """Load all the data sources with HNSCC data and provide an interface to them."""

        super().__init__(cancer_type="pancanhnscc", version=version, no_internet=no_internet)

        self._datasets["bcm"] = BcmHnscc(no_internet=no_internet, version=self._get_version("bcm"))
        self._datasets["broad"] = BroadHnscc(no_internet=no_internet, version=self._get_version("broad"))
        self._datasets["mssm"] = MssmClinical(no_internet=no_internet, version=self._get_version("mssm"), filter_type='pancanhnscc')
        self._datasets["pdc"] = PdcHnscc(no_internet=no_internet, version=self._get_version("pdc"))
        self._datasets["umich"] = UmichHnscc(no_internet=no_internet, version=self._get_version("umich"))
        self._datasets["washu"] = WashuHnscc(no_internet=no_internet, version=self._get_version("washu"))
        self._datasets["harmonized"] = Harmonized(no_internet=no_internet, version=self._get_version("harmonized"), filter_type='pancanhnscc')
        
        join_dict = {k: v._data for k, v in self._datasets.items()}
        self._joining_dataset = JoiningDataset(join_dict)
        
