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

# Exceptions
class CptacError(Exception):
    """Base class for all exceptions we'll raise."""
    pass

class NoInternetError(CptacError):
    """No internet."""
    pass

class InvalidParameterError(CptacError):
    """Invalid parameter."""
    pass

class AmbiguousLatestError(InvalidParameterError):
    """They pass "latest" for a version parameter, but index latest does not match latest version locally installed."""
    pass

class FileError(CptacError):
    """Base class for data-related errors."""
    pass

class DatasetNotInstalledError(FileError):
    """They requested a dataset they haven't installed."""
    pass

class DataVersionNotInstalledError(FileError):
    """They requested a version they haven't installed of a dataset."""
    pass

class MissingFileError(FileError):
    """A data file was missing."""
    pass

class DataError(CptacError):
    """Something was wrong with the data."""
    pass

class ReindexMapError(DataError):
    """Error in mapping keys when reindexing a dataframe."""
    pass

# Warnings
class CptacWarning(UserWarning):
    """Base class for all warnings we'll generate."""
    pass

class FailedReindexWarning(CptacWarning):
    """Error reindexing a dataframe."""
    pass

class NotApplicableFilterWarning(CptacWarning):
    """Filter value for multiple mutations existed for the dataset, but not for that particular gene."""
    pass

class NonexistentOmicsKeyWarning(CptacWarning):
    """Key for selecting from an omics dataframe didn't exist. Column created, but filled with NaN."""
    pass

class OldDataVersionWarning(CptacWarning):
    """They're using an old data version."""
    pass

class OldPackageVersionWarning(CptacWarning):
    """They're using an old version of the package."""
    pass

class DownloadingNewLatestWarning(CptacWarning):
    """Downloading a new latest data version. If they want to use an old version, they'll have to manually specify it."""
    pass
