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
