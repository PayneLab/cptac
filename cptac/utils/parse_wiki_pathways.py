import pandas as pd
import sys
import os
from time import time
from xml.dom import minidom

#### FORMAT GENES PER PATHWAY DICTIONARY SO IT IS A COMPLETE TABLE AND CAN BE CONVERTED DIRECTLY
#### INTO A DATAFRAME (df = pd.DataFrame.from_dict(dictionary))



def fixNewLines(s):
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


def fixTabs(s):
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


def fixParen(s):
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


def fixParsingErrors(interactingList):
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
            newlineFreeList = fixNewLines(item)
            for n in newlineFreeList:
                tabFreeList = fixTabs(n)
                for t in tabFreeList:
                    noParen = fixParen(t)
                    fixedList.extend(noParen)
    return list(set(fixedList))


def getUniprotGeneList(pathToUniprot):
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


def intersectWithUniprot(interactingGenes, uniprotGenes):
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



def getGenesPerPathwayDict(pathToWikiPathwaysRelease, pathToUniprot):
	"""
	@Param pathToWikiPathwaysRelease:
		A path to the WikiPathways release file. This is a global variable, declared at the top of this file.

	@Return:
		A dictionary of pathway names and the genes in each pathway. Key = pathway name. Value = list of proteins in the pathway.

	This function makes a dictionary of the WikiPathways release xml file. In this dictionary, the keys are the pathways names and the values are lists of proteins in that pathway.
	"""
    path = pathToWikiPathwaysRelease
    genesPerPathway = {}
    uniprotGenes = getUniprotGeneList(pathToUniprot)
    for fileName in os.listdir(path):
        if fileName == ".DS_Store":
            continue
        #filePath = path + "/" + fileName
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

        pathwayGenes = fixParsingErrors(pathwayGenes)
        pathwayGenes = intersectWithUniprot(pathwayGenes, uniprotGenes)
        # pathwayGenes = intersectWithUniprot(pathwayGenes)
        genesPerPathway[pathwayName] = pathwayGenes
    return genesPerPathway



def getAllPathways(genesPerPathway):
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


def getAllGenes(genesPerPathway):
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
    allGenes = fixParsingErrors(allGenes)
    return allGenes


def ParseWikiPathwaysToTabFile():
	"""
	This function creates the dataframe of pathways and proteins. The columns are pathway names and the rows are protein names. At any given index in the dataframe, if the corresponding protein (row) is in the
	corresponding pathway (column), the index has a value of True. Otherwise, the index value is false. This function then saves the dataframe as a tab separated file in the location provided by the pathToDataframeFile
	global variable declared at the top of this file.
	"""
    print("ParseWikiPathwaysToTabFile")
    start = time()
    ### 1. get commandline input
    # This was done at the top of the file with sys.argv[x].
    pathToUniprot = sys.argv[1]
    pathToWikiPathwaysRelease = sys.argv[2]
    pathToDataframeFile = sys.argv[3]

    ### 2. Call some accessory functions just to set up data structures
    genesPerPathway = getGenesPerPathwayDict(pathToWikiPathwaysRelease, pathToUniprot)  # this is returning a dictionary
    allPathways = getAllPathways(genesPerPathway)
    allGenes = getAllGenes(genesPerPathway)

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
    df.to_csv(pathToDataframeFile, sep=",")
    end = time()
    print(end - start)


###############################
# This program gets three variables from the command line
# global variables from command line
# sys.argv[1]: path to uniprot file (used to fix parsing errors in the WikiPathways relsease xml)
# sys.argv[2]: path to wikipathways directory .gpml
# sys.argv[3]: path to tabbed dataframe file


################################
# calls the function above to save the dataframe

if len(sys.argv) < 3:
    print("USAGE: ParseWikiPathways.py <path_to_uniprot_file> <path_to_wikipathways> <path_to_tabbed_dataframe>")
    sys.exit(0)
ParseWikiPathwaysToTabFile()
