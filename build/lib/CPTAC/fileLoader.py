class FileLoader:
    def __init__(self, fileName):
        self.fileName = fileName
    def getFileName(self):
        return self.fileName
    def setFileName(self, fileName):
        self.fileName = fileName
    def readFile(self):
        return open(self.fileName, 'r')
    def writeFile(self):
        return open(self.fileName, 'w')
