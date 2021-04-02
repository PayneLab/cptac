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
#from .bcmhnscc import BcmHnscc
from .umichhnscc import UmichHnscc
from .washuhnscc import WashuHnscc

SOURCES = [
    "mssmclinical",
    "umichgbm",
    "washuhnscc",
]

class PancanHnscc(PancanDataset):

    def __init__(self, versions="latest", no_internet=False):
        """Load all the data sources with HNSCC data and provide an interface to them."""

        super().__init__(cancer_type="pancanhnscc", versions=versions, no_internet=no_internet)

        def get_version(source):
            if versions == "latest":
                return versions
            else:
                return versions[source]
            
        self._datasets["mssm"] = MssmClinical(no_internet=no_internet, version=get_version("mssm"), filter_type='pancanhnscc')
        #self._datasets["bcm"] = BcmHnscc(no_internet=no_internet, version=get_version("bcm"))
        self._datasets["umich"] = UmichHnscc(no_internet=no_internet, version=get_version("umich"))
        
        
