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

class Error(Exception):
    """Base class for cptac exceptions."""

class InvalidVersionError(Error):
    """For when an invalid version is requested."""

    def __init__(self, invalid_version, valid_versions):
        """Initialize an InvalidVersionError with a message about valid options.

        Parameters:
        invalid_version (str): The invalid version the user passed.
        valid_versions (list or iterable of str): The valid options for versions.

        Returns:
        InvalidVersionError: Initialized with a message telling valid version options.
        """
        message = "'{}' is an invalid version for this dataset. Valid versions:".format(invalid_version)
        for version in valid_versions:
            message = message + "\n\t'{}'".format(version)

        super().__init__(message)

