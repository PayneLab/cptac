#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

class DataSet:

    def __init__(self):
        self.data = {}

    # Get metadata dataframes
    def get_clinical(self):
        """Get clinical dataframe."""
        return self._get_dataframe("clinical")

    def get_derived_molecular(self):
        """Get derived_molecular dataframe."""
        return self._get_dataframe("derived_molecular")

    def get_experimental_setup(self):
        """Get experimental_setup dataframe."""
        return self._get_dataframe("experimental_setup")

    def get_treatment(self):
        """Get treatment dataframe."""
        return self._get_dataframe("treatment")

    # Get omics dataframes
    def get_acetylproteomics(self):
        """Get acetylproteomics dataframe."""
        return self._get_dataframe("acetylproteomics")

    def get_proteomics(self):
        """Get proteomics dataframe."""
        return self._get_dataframe("proteomics")

    def get_transcriptomics(self):
        """Gets transcriptomics dataframe."""
        return self._get_dataframe("transcriptomics")

    def get_circular_RNA(self):
        """Gets circular_RNA dataframe."""
        return self._get_dataframe("circular_RNA")

    def get_miRNA(self):
        """Gets miRNA dataframe."""
        return self._get_dataframe("miRNA")

    def get_CNA(self):
        """Get the CNA dataframe."""
        return self._get_dataframe("CNA")

    def get_phosphoproteomics(self):
        """Gets the phosphoproteomics dataframe."""
        return self._get_dataframe("phosphoproteomics")

    def get_phosphoproteomics_gene(self):
        """Gets the phosphoproteomics_gene dataframe. The gene level phosphorylation measurement is an aggregate metric which potentially averages to
    gether individual measurements of different sites. Use get_phosphoproteomics() to view the data for individual sites."""
        return self._get_dataframe("phosphoproteomics_gene")

    # Get mutations dataframes
    def get_mutations(self):
        """Get the somatic_mutation dataframe."""
        return self._get_dataframe("somatic_mutation")

    def get_mutations_binary(self):
        """Gets the somatic_mutation_binary dataframe, which has a binary value indicating, for each location on each gene, whether there was a mutat
    ion in that gene at that location, for each sample."""
        return self._get_dataframe("somatic_mutation_binary")

    # Get map of Sample_ID to Sample_Status
    def get_sample_status_map(self):
        """Get a pandas Series from the clinical dataframe, with sample ids as the index, and each sample's status (tumor or normal) as the values."""
        clinical = get_clinical()
        raw_map = clinical["Proteomics_Tumor_Normal"]
        parsed_map = raw_map.where(raw_map == "Tumor", other="Normal") # Replace various types of normal (Adjacent_normal, Myometrium_normal, etc.) with just "Normal"
        parsed_map.name = "Sample_Status"
        return parsed_map

    # Utilities functions

    # Help functions
    def list_api(self):
        """Print docstrings for all accessible functions."""
        help(__name__)

    def list_data(self):
        """Print list of loaded data frames and dimensions."""
        print("Below are the available endometrial data frames contained in this package:")
        for dataframe in data.values():
            print("\t", dataframe.name)
            print("\t", "\t", "Dimensions:", dataframe.shape)

    # "Private" functions
    def _get_dataframe(self, name):
        """Check if a dataframe with the given name exists, and return it if it does.

        Parameters:
        name (str): The name of the dataframe to get.

        Returns:
        pandas DataFrame: The desired dataframe, if it exists in this dataset.
        """
        if name in self.data.keys():
            return self.data[name]
        else:
            print("{} dataframe not included in this dataset.".format(name))
            return
