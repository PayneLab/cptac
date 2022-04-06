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

from cptac.cancers.brca import Brca
from cptac.cancers.ccrcc import Ccrcc
from cptac.cancers.coad import Coad
from cptac.cancers.gbm import Gbm
from cptac.cancers.hnscc import Hnscc
from cptac.cancers.lscc import Lscc
from cptac.cancers.luad import Luad
from cptac.cancers.ov import Ov
from cptac.cancers.pdac import Pdac
from cptac.cancers.ucec import Ucec


def list_datasets(print_list=True):
    """Print available datasets in the cptac.pancan module.
    Parameters:
    print_list (bool, optional): Whether to print the list. Default is True. Otherwise, it's returned as a string.
    """

    datasets = [
        "Brca",
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