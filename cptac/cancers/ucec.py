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

from cptac.cancers.awg.awgucec import AwgUcec
from cptac.cancers.awgconf.awgconfucec import AwgConfUcec
from cptac.cancers.bcm.bcmucec import BcmUcec
from cptac.cancers.broad.broaducec import BroadUcec
from cptac.cancers.pdc.pdcucec import PdcUcec
from cptac.cancers.umich.umichucec import UmichUcec
from cptac.cancers.washu.washuucec import WashuUcec
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized


class Ucec(Cancer):

    def __init__(self, no_internet=False):
        """Load all the data sources with UCEC data and provide an interface to them."""

        super().__init__(cancer_type="ucec")

        self._sources["awg"] = AwgUcec(version="latest", no_internet=no_internet)
        self._sources["awgconf"] = AwgConfUcec(version="latest", no_internet=no_internet)
        self._sources["bcm"] = BcmUcec(version="latest", no_internet=no_internet)
        self._sources["broad"] = BroadUcec(version="latest", no_internet=no_internet)
        self._sources["mssm"] = Mssm(filter_type='ucec', version="latest", no_internet=no_internet)
        self._sources["pdc"] = PdcUcec(version="latest", no_internet=no_internet)
        self._sources["umich"] = UmichUcec(version="latest", no_internet=no_internet)
        self._sources["washu"] = WashuUcec(version="latest", no_internet=no_internet)
        self._sources["harmonized"] = Harmonized(filter_type='ucec', version="latest", no_internet=no_internet)
