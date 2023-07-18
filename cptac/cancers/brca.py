#   Copyright 2023 Samuel Payne sam_payne@byu.edu
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

from cptac.cancers.bcm.bcmbrca import BcmBrca
from cptac.cancers.broad.broadbrca import BroadBrca
from cptac.cancers.umich.umichbrca import UmichBrca
from cptac.cancers.washu.washubrca import WashuBrca
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized

class Brca(Cancer):
    """Class to manage BRCA data from various sources.
    
    Extends the Cancer class and initializes the BRCA data from different sources.
    """

    def __init__(self, no_internet=False):
        """Initializes the Brca object.
        
        Args:
            no_internet (bool): If True, the object will not attempt to download data from the internet. Default is False.
        """
        super().__init__(cancer_type="brca")

        classes = {'bcm': BcmBrca, 'broad': BroadBrca, 'mssm': Mssm, 
                   'umich': UmichBrca, 'washu': WashuBrca, 'harmonized': Harmonized}
        
        self._sources = {key: value(no_internet=no_internet) if key not in ['mssm', 'harmonized'] 
                         else value(filter_type='brca', no_internet=no_internet) 
                         for key, value in classes.items()}
