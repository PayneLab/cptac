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

import cptac
from cptac.exceptions import InvalidParameterError

def download(dataset, version="latest", redownload=False):

    if dataset == "pancanbrca":
        cptac.download("bcmbrca", version=version, redownload=redownload)
        cptac.download("broadbrca", version=version, redownload=redownload)
        cptac.download("mssmclinical", version=version, redownload=redownload)
        cptac.download("umichbrca", version=version, redownload=redownload)
        cptac.download("washubrca", version=version, redownload=redownload)
    else:
        raise InvalidParameterError(f"{dataset} is not a valid dataset.")
