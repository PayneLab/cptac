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

from cptac.exceptions import DataTypeNotInSourceError, InvalidDataVersionError

class Source:

    def __init__(self, source, cancer, version, valid_versions):
        
        self.source = source
        self.cance.type = cancer
        self._data = list()
        self.load_functions = dict()
        self.version = self.set_version(version)
        self.valid_versions = valid_versions
        
    def _get_df(self, version, df_type):

        # if that df hasn't been loaded yet, load it
        if df_type not in self._data:

            # check to see if the df type requested is availabe from this source
            if df_type not in self.load_functions:
                raise DataTypeNotInSourceError(f"The {self.source} source does not have {df_type} data for {self.cancer_type} cancer.")
            else:
                # call the load function for that df type
                self.load_functions[df_type](version)

        return self._data[df_type]

    def _set_version(self, version):

        # check if version is valid
        if version not in self.valid_versions:
            raise InvalidDataVersionError(f"{version} is not a valid version. These are the valid versions: {self.valid_versions}")

    def _version(self):
        return self.version