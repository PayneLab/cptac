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

import time
import cptac
import logging
import pytest
from cptac.exceptions import InvalidParameterError

from .conftest import get_cancer_inputs


# I would include other attributes to each test. Ex, checking file location, 
#   negative test for getting data before downloading, etc.
@pytest.mark.parametrize("cancer, source, datatype", get_cancer_inputs())
def test_redownload(cancer, source, datatype):
    print(f"{cancer}, {source}, {datatype}")
    time.sleep(1)
    assert cptac.download(sources={source: datatype}, cancers=cancer, redownload=True)

@pytest.mark.parametrize("cancer, source, datatype", get_cancer_inputs())
def test_redownload(cancer, source, datatype):
    time.sleep(1)
    assert cptac.download(sources={source: datatype}, cancers=cancer, redownload=True)
    
