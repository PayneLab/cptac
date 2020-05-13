###############################
# This program takes a single argument from the command line
# sys.argv[1]: path to the wikipathways directory to obtain the lastest information.
# To search for the latest directory see https://www.wikipathways.org/index.php/Download_Pathways


# sys.argv[1]: path to uniprot file (used to fix parsing errors in the WikiPathways relsease xml) -- Called Uniprot_Proteome.tsv
# sys.argv[2]: path to wikipathways directory (https://www.wikipathways.org/index.php/Download_Pathways)
# sys.argv[3]: path to tabbed dataframe file
################################

# TODO: Make the output file the name we expect it to be in get_interacting_proteins
# TODO: Update the docstring

import pandas as pd
import sys
import os
from xml.dom import minidom


class ParseWikiPathways:
    def __init__(self):
        pass

    def ParseWikiPathwaysToTabFile(self):
        ### 1. get commandline input
        pathToUniprot = 'Uniprot_Proteome.tsv'
        pathToWikiPathwaysRelease = sys.argv[1]
        pathToDataframeFile = 'WikipathwaysDataframe.tsv'

        ### 2. Call some accessory functions just to set up data structures
        genesPerPathway = self.getGenesPerPathwayDict(pathToWikiPathwaysRelease,
                                                 pathToUniprot)  # this is returning a dictionary
        allPathways = self.getAllPathways(genesPerPathway)
        allGenes = self.getAllGenes(genesPerPathway)

        ### 3. Now we make the data frame
        # df = pd.DataFrame(genesPerPathway, columns = allPathways, index = allGenes)
        df = pd.DataFrame(columns=allPathways, index=allGenes)
        df = df.fillna(False)

        for pathway in df.columns:
            for gene in df.index:
                if gene in genesPerPathway[pathway]:
                    # df.replace(to_replace = False, value = True)
                    df.at[gene, pathway] = True
        ### 4. Now we need to save out the dataframe to a tabbed file
        df.to_csv(pathToDataframeFile, sep="\t")


    def fixNewLines(self, s):
        """
        @Param s:
            A protein name that is checked for new lines.

        @Return:
            A list of protein names.

        This function us used to fix parsing errors from the xml file of the WikiPathways release.
        If the parameter s contains a new line, it returns the protein names resulting from splitting s with a new line.
        If the parameter s doesn't contain a new line, it returns a list containing the only the parameter s.
        """
        if "\n" in s:
            nList = s.split("\n")
        else:
            nList = []
            nList.append(s)
        return nList

    def fixTabs(self, s):
        """
        @Param s:
            A protein name that is checked for tabs.

        @Return:
            A list of protein names.

        This function us used to fix parsing errors from the xml file of the WikiPathways release.
        If the parameter s contains a tab, it returns the protein names resulting from splitting s with a tab.
        If the parameter s doesn't contain a tab, it returns a list containing the only the parameter s.
        """

        if "\t" in s:
            tList = s.split("\t")
        else:
            tList = []
            tList.append(s)
        return tList

    def fixParen(self, s):
        """
        @Param s:
            A protein name that is checked for parenthesis.

        @Return:
            A list of protein names.

        This function us used to fix parsing errors from the xml file of the WikiPathways release.
        If the parameter s contains a parenthesis, it returns the protein names resulting from splitting s with a parenthesis.
        If the parameter s doesn't contain a parenthesis, it returns a list containing the only the parameter s.
        """

        fixedList = []
        list1 = s.split(")")
        for i in list1:
            i = i.strip()
            list2 = i.split("(")
            for j in list2:
                j = j.strip()
                if j != "":
                    fixedList.append(j)
        return fixedList


    def fixParsingErrors(self, interactingList):
        """
        @Param interactingList:
            A list of the names of interacting proteins obtained form each pathway found in the file containing the WikiPathways Release

        @Return:
            A list of the names of interacting proteins from the interactingList parameter, after parsing errors have been fixed

        This function fixes parsing errors from the WikiPathways Release file. This function will eliminate new lines, tabs, and parenthesis using the fixNewLines, fixTabs, and fixParen functions above.
        Where needed, this function will also separate protein names that were combinded with newlines, tabs, or parenthesis
        """


        fixedList = []
        for i in interactingList:
            i = i.strip()
            noSpaces = i.split(" ")
            for item in noSpaces:
                newlineFreeList = self.fixNewLines(item)
                for n in newlineFreeList:
                    tabFreeList = self.fixTabs(n)
                    for t in tabFreeList:
                        noParen = self.fixParen(t)
                        fixedList.extend(noParen)
        return list(set(fixedList))


    def getUniprotGeneList(self, pathToUniprot):
        """
        @Param pathToUniprot:
            A path to the Uniprot Proteome file that is obtained from the commandLine. This is a global variable declared at the top of this file.

        @Return:
            A list of all the Uniprot genes in the Uniprot Proteome file.

        This function takes the path to the Uniprot Proteome file and returns a list of all the genes in that file.
        This list is later intersected with the list of WikiPathways proteins to eliminate non-genes from the list of WikiPathways proteins.
        """


        filePath = pathToUniprot
        bioplex_interactions = pd.read_csv(filePath, sep='\t')
        genes = bioplex_interactions['Gene names'].tolist()
        uniprotList = set()
        for gene in genes:
            if str(gene) != 'nan':
                gene_arr = gene.split()
                for gene in gene_arr:
                    uniprotList.add(gene)
        return list(uniprotList)


    def intersectWithUniprot(self, interactingGenes, uniprotGenes):
        """
        @Param interactingGenes:
            A list of the protein names in a certain pathway.

        @Return:
            A list of protein names, after the intersection with Uniprot.

        This function takes a list of interacting proteins and intersects them with a list of Uniprot genes. This function also eliminates any possible repeats
        """


        geneSet = set(interactingGenes)
        uniprotSet = set(uniprotGenes)
        intersectSet = geneSet.intersection(uniprotSet)
        return list(intersectSet)


    def getGenesPerPathwayDict(self, pathToWikiPathwaysRelease, pathToUniprot):
        """
        @Param pathToWikiPathwaysRelease:
            A path to the WikiPathways release file. This is a global variable, declared at the top of this file.

        @Return:
            A dictionary of pathway names and the genes in each pathway. Key = pathway name. Value = list of proteins in the pathway.

        This function makes a dictionary of the WikiPathways release xml file. In this dictionary, the keys are the pathways names and the values are lists of proteins in that pathway.
        """


        path = pathToWikiPathwaysRelease
        genesPerPathway = {}
        uniprotGenes = self.getUniprotGeneList(pathToUniprot)
        for fileName in os.listdir(path):
            if fileName == ".DS_Store":
                continue
            # filePath = path + "/" + fileName
            filePath = f'{path}/{fileName}'
            pathwayGenes = []

            xmldoc = minidom.parse(filePath)
            pathway = xmldoc.getElementsByTagName("Pathway")[0]
            pathwayName = pathway.getAttribute("Name")

            dataNodes = pathway.getElementsByTagName("DataNode")
            for node in dataNodes:
                if node.getAttribute("Type") == "GeneProduct":
                    geneName = node.getAttribute("TextLabel")
                    pathwayGenes.append(geneName)

            pathwayGenes = self.fixParsingErrors(pathwayGenes)
            pathwayGenes = self.intersectWithUniprot(pathwayGenes, uniprotGenes)
            # pathwayGenes = intersectWithUniprot(pathwayGenes)
            genesPerPathway[pathwayName] = pathwayGenes
        return genesPerPathway


    def getAllPathways(self, genesPerPathway):
        """
        @Param genesPerPathway:
            A dictionary of pathway names and lists of proteins. Key = pathway name. Value = list of proteins in that pathway.

        @Return:
            A list of all the pathways (keys) in the genesPerPathway dicitonary.

        This function takes a dictionary of pathways and proteins in each pathway and return a list of all the pathways in the dictionary. This list is later used to make the column names in the dataframe.
        """


        allPathways = genesPerPathway.keys()
        allPathways = list(allPathways)
        return allPathways


    def getAllGenes(self, genesPerPathway):
        """
        @Param genesPerPathway:
            A dictionary of pathway names and lists of proteins. Key = pathway name. Value = list of proteins in that pathway.

        @Return:
            A list of all the proteins (values) in the genesPerPathway dicitonary.

        This function takes a dictionary of pathways and proteins in each pathway and return a list of all the proteins in the dictionary. This list is later used to make the row names in the dataframe.
        """


        allGenes = []
        for geneList in genesPerPathway.values():
            allGenes.extend(geneList)
        allGenes = list(set(allGenes))
        allGenes = self.fixParsingErrors(allGenes)
        return allGenes

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("USAGE: ParseWikiPathways.py <path_to_wikipathways_directory>")
        print(" This program takes a single argument from the command line \n \
        sys.argv[1]: path to the wikipathways directory to obtain the lastest information.\n \
        To search for the latest directory see https://www.wikipathways.org/index.php/Download_Pathways")
        sys.exit(0)
    p = ParseWikiPathways()
    p.ParseWikiPathwaysToTabFile()
