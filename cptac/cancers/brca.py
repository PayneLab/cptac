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

from cptac.cancers.awg.awgbrca import AwgBrca
from cptac.cancers.bcm.bcmbrca import BcmBrca
from cptac.cancers.broad.broadbrca import BroadBrca
from cptac.cancers.pdc.pdcbrca import PdcBrca
from cptac.cancers.umich.umichbrca import UmichBrca
from cptac.cancers.washu.washubrca import WashuBrca
from cptac.cancers.mssm.mssmclinical import MssmClinical
from cptac.cancers.harmonized.harmonized import Harmonized


class Brca(Cancer):

    def __init__(self, no_internet=False):
        """Load all the data sources with BRCA data and provide an interface to them."""

        super().__init__(cancer_type="brca")

        self._sources["awg"] = AwgBrca(no_internet=no_internet, version="latest")
        self._sources["bcm"] = BcmBrca(no_internet=no_internet, version="latest")
        self._sources["broad"] = BroadBrca(no_internet=no_internet, version="latest")
        self._sources["mssm"] = MssmClinical(no_internet=no_internet, version="latest", filter_type='pancanbrca')
        self._sources["pdc"] = PdcBrca(no_internet=no_internet, version="latest")
        self._sources["umich"] = UmichBrca(no_internet=no_internet, version="latest")
        self._sources["washu"] = WashuBrca(no_internet=no_internet, version="latest")
        self._sources["harmonized"] = Harmonized(no_internet=no_internet, version="latest", filter_type= 'pancanbrca')
        
        join_dict = {k: v._data for k, v in self._sources.items()}
        self._joining_dataset = JoiningDataset(join_dict)
        
        self._pancan_unionize_indices()
