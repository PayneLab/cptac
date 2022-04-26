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
from cptac.tools.joiningdataset import JoiningDataset

from cptac.cancers.awg.awgucec import AwgUcec
from cptac.cancers.bcm.bcmucec import BcmUcec
from cptac.cancers.broad.broaducec import BroadUcec
from cptac.cancers.pdc.pdcucec import PdcUcec
from cptac.cancers.umich.umichucec import UmichUcec
from cptac.cancers.washu.washuucec import WashuUcec
from cptac.cancers.mssm.mssmclinical import MssmClinical
from cptac.cancers.harmonized.harmonized import Harmonized

#List Sources to be downloaded
SOURCES = [
    "bcmucec",
    "mssmclinical",
    "pdcucec",
    "umichucec",
    "washuucec",
    "broaducec",
    "harmonized"
    
]

class Ucec(Cancer):

    def __init__(self, version="latest", no_internet=False):
        """Load all the data sources with UCEC data and provide an interface to them."""

        super().__init__(cancer_type="pancanucec", version=version, no_internet=no_internet)

        self._datasets["awg"] = AwgUcec(no_internet=no_internet, version=self._get_version("awg"))
        self._datasets["bcm"] = BcmUcec(no_internet=no_internet, version=self._get_version("bcm"))
        self._datasets["broad"] = BroadUcec(no_internet=no_internet, version=self._get_version("broad"))
        self._datasets["mssm"] = MssmClinical(no_internet=no_internet, version=self._get_version("mssm"), filter_type='pancanucec')
        self._datasets["pdc"] = PdcUcec(no_internet=no_internet, version=self._get_version("pdc"))
        self._datasets["umich"] = UmichUcec(no_internet=no_internet, version=self._get_version("umich"))
        self._datasets["washu"] = WashuUcec(no_internet=no_internet, version=self._get_version("washu"))
        self._datasets["harmonized"] = Harmonized(no_internet=no_internet, version=self._get_version("harmonized"), filter_type='pancanucec')
        
        join_dict = {k: v._data for k, v in self._datasets.items()}
        self._joining_dataset = JoiningDataset(join_dict)
        
        self._pancan_unionize_indices()
