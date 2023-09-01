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

# Import various BRCA data sources
from cptac.cancers.bcm.bcmbrca import BcmBrca
from cptac.cancers.broad.broadbrca import BroadBrca
from cptac.cancers.umich.umichbrca import UmichBrca
from cptac.cancers.washu.washubrca import WashuBrca
from cptac.cancers.mssm.mssm import Mssm
from cptac.cancers.harmonized.harmonized import Harmonized

class Brca(Cancer):
    """Manages BRCA (Breast Cancer) data from various sources.
    
    This class extends the base Cancer class and initializes the BRCA data from 
    a variety of sources including BCM, Broad Institute, MSSM, University of Michigan, 
    Washington University, and a Harmonized dataset. 

    Attributes:
        _sources (dict): A dictionary holding data from different sources.
    """

    def __init__(self, no_internet=False):
        """Initializes the Brca object.
        
        Args:
            no_internet (bool): If True, the object will not attempt to download data from the internet. 
                                Default is False.

        Raises:
            ValueError: If the 'no_internet' argument is not of boolean type.
        """
        if not isinstance(no_internet, bool):
            raise ValueError("The 'no_internet' argument must be of boolean type.")

        super().__init__(cancer_type="brca")

        classes = {'bcm': BcmBrca, 'broad': BroadBrca, 'mssm': Mssm, 
                   'umich': UmichBrca, 'washu': WashuBrca, 'harmonized': Harmonized}
        
        # Populate the _sources attribute with data from different sources
        self._sources = {key: value(no_internet=no_internet) if key not in ['mssm', 'harmonized'] 
                         else value(filter_type='brca', no_internet=no_internet) 
                         for key, value in classes.items()}
