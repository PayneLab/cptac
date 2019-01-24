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

class MolecularData:
    def __init__(self, proteomics, transcriptome, cna, phosphoproteomics):
        self.proteomics = proteomics
        self.transcriptome = transcriptome
        self.cna = cna
        self.phosphoproteomics = phosphoproteomics
    def get_proteomics(self):
        return self.proteomics
    def get_transcriptome(self):
        return self.transcriptome
    def get_cna(self):
        return self.cna
    def get_phosphoproteomics(self):
        return self.phosphoproteomics
