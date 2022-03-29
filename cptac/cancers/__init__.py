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

from cptac.file_download import get_box_token

from .file_download import download, download_pdc_id, list_pdc_datasets

from .pdcbrca import PdcBrca
from .pdcccrcc import PdcCcrcc
from .pdccoad import PdcCoad
from .pdcgbm import PdcGbm
from .pdchnscc import PdcHnscc
from .pdclscc import PdcLscc
from .pdcluad import PdcLuad
from .pdcov import PdcOv
from .pdcpdac import PdcPdac
from .pdcucec import PdcUcec

from .brca import Brca
from .ccrcc import Ccrcc
from .coad import Coad
from .gbm import Gbm
from .hnscc import Hnscc
from .lscc import Lscc
from .luad import Luad
from .ov import Ov
from .ucec import Ucec
from .pdac import Pdac

def list_datasets(print_list=True):
    """Print available datasets in the cptac.pancan module.
    Parameters:
    print_list (bool, optional): Whether to print the list. Default is True. Otherwise, it's returned as a string.
    """

    datasets = [
        "PancanBrca",
        "Ccrcc",
        "Coad",
        "Gbm",
        "Hnscc",
        "Lscc",
        "Luad",
        "Ov",
        "Ucec",
        "Pdac"
    ]

    str_result = "\n".join(datasets)

    if print_list:
        print(str_result)
    else:
        return str_result