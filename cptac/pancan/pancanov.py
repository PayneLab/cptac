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

from .mssmclinical import MssmClinical
from .pdcov import PdcOv
from .umichov import UmichOv
from .washuov import WashuOv
from .bcmov import BcmOv

SOURCES = [
    "bcmov",
    "mssmclinical",
    "pdcov",
    "umichov",
    "washuov",
]

class PancanOv(PancanDataset):

    def __init__(self, versions="latest", no_internet=False):
        """Load all the data sources with OV data and provide an interface to them."""

        super().__init__(cancer_type="pancanov", versions=versions, no_internet=no_internet)
        
        self._datasets["bcm"] = BcmOv(no_internet=no_internet, version=self._get_version("ov"))
        self._datasets["mssm"] = MssmClinical(no_internet=no_internet, version=self._get_version("mssm"), filter_type='pancanov')
        self._datasets["pdc"] = PdcOv(no_internet=no_internet, version=self._get_version("pdc"))
        self._datasets["umich"] = UmichOv(no_internet=no_internet, version=self._get_version("umich"))
        self._datasets["washu"] = WashuOv(no_internet=no_internet, version=self._get_version("washu"))
