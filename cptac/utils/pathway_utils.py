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

import os
import pandas as pd
import numpy as np
import sys
import urllib3
import json
import requests
import webbrowser
import warnings

from cptac.exceptions import HttpResponseError, InvalidParameterError, ParameterWarning

def get_interacting_proteins_string(protein, num_results=25):
    """
    @Param protein:
        The name of the protein that you want to generate a list of interacting proteins for.

    @Param num_results (default=25):
        The number of interacting proteins that you want to get.

    @Return:
        A pandas.Series of proteins known by the String api to be interacting partners with the specified protein.
        This list will always also contain the protein you were looking for interactors for.

    This method takes as a parameter the name of a protein. It then accesses the STRING database, through
    a call to their public API, and generates a list of proteins known to be interacting partners with the specified
    protein. Optional second parameter is num_results (which by default is 25), which specifies in the API call how many
    interacting partners to retrieve from the database. The list of interacting proteins is returned to the caller
    as a pandas.Series.
    """

    # Send query to the STRING API
    query_url = "https://string-db.org/api/json/network"

    params = {
        "identifiers": protein,
        "species": "9606", # Homo sapiens
        "limit": num_results,
    }

    query_resp = requests.get(query_url, params=params)

    # Check that the response came back good
    if query_resp.status_code != requests.codes.ok:
        raise HttpResponseError(f"Submitting your query to the STRING API returned an HTTP status {query_resp.status_code}. The content returned from the request may be helpful:\n{query_resp.content.decode('utf-8')}")    

    # Put the response data in a dataframe
    resp_df = pd.DataFrame(query_resp.json())

    # Get the unique values of columns we want, as a list
    interactors = resp_df["preferredName_A"].\
        append(resp_df["preferredName_B"]).\
        unique()

    # Make sure the protein they searched is in the output list
    if protein not in interactors:
        interactors = np.insert(interactors, 0, protein)

    # Sort and convert to series
    interactors = np.sort(interactors)
    interactors = pd.Series(interactors)

    return interactors

def get_interacting_proteins_biogrid(protein):
    """Queries the BioGRID API to get interacting proteins for a given protein, based on curated literature references.

    Parameters:
    protein: The name of the protein that you want to generate a list of interacting proteins for.

    Returns:
    pandas.DataFrame: The interacting proteins, ranked by the number of literature references supporting them.
    """

    # Post query to the BioGRID API
    query_url = "https://webservice.thebiogrid.org/interactions/"
    
    params = {
        "searchNames": "true",
        "geneList": protein,
        "includeInteractors": "true",
        "includeInteractorInteractions": "false",
        "interSpeciesExcluded": "true",
        "format": "json",
        "taxId": "9606",
        "start": "0",
        "accesskey": "0ff59dcf3511928e78aad499688381c9"
    }

    query_resp = requests.get(query_url, params=params)

    # Check that the response came back good
    if query_resp.status_code != requests.codes.ok:
        raise HttpResponseError(f"Submitting your query to the STRING API returned an HTTP status {query_resp.status_code}. The content returned from the request may be helpful:\n{query_resp.content.decode('utf-8')}")    

    elif len(query_resp.json()) == 0:
        raise InvalidParameterError(f"No interactors found for '{protein}'. Are you sure you entered the identifier correctly?")

    # Put the response data in a dataframe
    resp_df = pd.DataFrame(query_resp.json()).transpose()

    # Get a list of all the interactors, and rank them by how many references each has
    interactors = resp_df["OFFICIAL_SYMBOL_A"].\
    where(resp_df["OFFICIAL_SYMBOL_A"] != protein, other=resp_df["OFFICIAL_SYMBOL_B"]).\
    value_counts().\
    to_frame("num_references")

    interactors.index.name = "protein"

    return interactors

def get_interacting_proteins_bioplex(protein, secondary_interactions=False):
    """
    @Param protein:
        The name of the protein that you want to generate a list of interacting proteins for.

    @Return:
        A list of proteins which are interacting partners with the specified protein, according to the bioplex data table.
        Returns None if specified protein isn't found, or no interacting partners are found.

    This method takes as a parameter the name of a protein. It then accesses the bioplex data table and returns a list of any protein found to be an interacting partner to the given protein.

    The Bioplex data table is the "BioPlex 3.0 Interactions (293T Cells)" file for the HEK293T cell line, downloaded from <https://bioplex.hms.harvard.edu/interactions.php>. The direct download link is <https://bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv>. IMPORTANT: After downloading the file, you need to compress it using gzip.
    """

    path_here = os.path.abspath(os.path.dirname(__file__))
    file_name = os.path.join("data", "BioPlex_293T_Network_10K_Dec_2019.tsv.gz")
    file_path = os.path.join(path_here, file_name)

    # Read in the file, then sort to prioritize the interactions with the
    # highest pInt (probability of interacting) and lowest pNI (probability of
    # no interaction). Then, reset the index so that index number will
    # correspond to rank under this sorting scheme.
    bioplex_interactions = pd.read_csv(file_path, sep='\t').\
        sort_values(by=["pInt", "pNI"], ascending=[False, True]).\
        reset_index()

    # Combine the two series, sorting by index so that they end up in order of
    # priority, in case we filter in the future.
    all_interactions = bioplex_interactions["SymbolA"].\
        append(bioplex_interactions["SymbolB"]).\
        sort_index().\
        tolist()

    if secondary_interactions:
        secondary_interactions_list = []
        for interaction in all_interactions:
            secondary = get_interacting_proteins_bioplex(interaction, False)
            for si in secondary:
                secondary_interactions_list.append(si)

        for asi in secondary_interactions_list:
            if asi not in all_interactions:
                all_interactions.append(asi)

    if len(all_interactions) > 0:
        return all_interactions
    else:
        return None
    
def get_interacting_proteins_wikipathways(protein):
    """
    @param protein:
    	String. The name of the protein
    @Return:
    	A list of proteins known by the most recent WikiPathways download to be interacting parters with the specified protein.
    	Returns None if specified protein is not found in the WikiPathways dataframe (which was intersected with Uniprot).

    This function takes a path to WikiPathways Dataframe file and protein name and returns a list of all the proteins that interact with it, using the pathways from the WikiPathways relsease file.
    This function loads the WikiPathways dataframe, and iterates through the row labelled with that protein name, return every protein in a pathway that also contains that protein.
    """

    path_here = os.path.abspath(os.path.dirname(__file__))
    data_dir_name = "data"
    file_name = "WikiPathwaysDataframe.tsv.gz"
    file_path = os.path.join(path_here, data_dir_name, file_name)
    proteinName = protein

    df = pd.read_csv(file_path, sep="\t", index_col=0)

    if (proteinName in df.index):
    	row = df.loc[proteinName]
    	filtered_df = df.loc[:, row.values.tolist()]

    	def has_true(values):
    		for val in values:
    			if val == True:
    				return True
    		return False

    	filtered_df_final = filtered_df.loc[filtered_df.apply(lambda row: has_true(row.values.tolist()), axis=1), :]
    	return filtered_df_final.index.tolist()

    return list()  # The protein was not found.

def list_pathways_wikipathways():
    """
    @ Return:
    	A list of all the WikiPathways pathways

    Uses the WikipathwaysDataFrame to return a list of all the possible pathways found.
    """
    path_here = os.path.abspath(os.path.dirname(__file__))
    data_dir_name = "data"
    file_name = "WikiPathwaysDataframe.tsv.gz"
    file_path = os.path.join(path_here, data_dir_name, file_name)

    df = pd.read_csv(file_path, sep="\t", index_col=0)
    return list(df.columns)

def get_pathways_with_proteins(proteins, database, reactome_resource="UniProt", quiet=False):
    """Query either the Reactome REST API or the WikiPathways downloaded dataframe to find pathways containing a particular gene or protein.

    Parameters:
    proteins (str or list of str): The protein(s) to look for matches to.
    database (str): The database to use; either 'reactome' or 'wikipathways'.
    reactome_resource (str, optional): If using Reactome, this is the resource the identifier(s) come from. Default is UniProt. Other options include HGNC, Ensembl, and GO. For more options, consult <https://reactome.org/content/schema/objects/ReferenceDatabase>. This parameter is meaningless if using WikiPathways.
    quiet (bool, optional): Whether to suppress warnings issued when identifiers are not found. Default False.

    Returns:
    pandas.DataFrame: A table of pathways containing the given genes or proteins, with pathway names and, if using Reactome, their Reactome identifiers (which are needed for the pathway_overlay function).
    """

    # Process string input
    if isinstance(proteins, str):
        proteins = [proteins]

    if database.lower() == "reactome":

        # Set headers and params
        headers = {"accept": "application/json"}
        params = {"species": "Homo sapiens"}

        # Loop over proteins and get the interacting pathways
        all_pathway_df = pd.DataFrame()
        for id in proteins:
            url = f"https://reactome.org/ContentService/data/mapping/{reactome_resource}/{id}/pathways"
            resp = requests.get(url, headers=headers, params=params)

            # Check that the response came back good
            if resp.status_code == 404:
                try:
                    msg = resp.json()["messages"]
                except (json.JSONDecodeError, KeyError):
                    raise HttpResponseError(f"Your query returned an HTTP status {resp.status_code}. The content returned from the request may be helpful:\n{resp.content.decode('utf-8')}") from None
                else:
                    if not quiet:
                        warnings.warn(f"The query for '{id}' returned HTTP 404 (not found). You may have mistyped the gene/protein ID or the reactome_resource name. The server gave the following message: {msg}", ParameterWarning, stacklevel=2)
                    continue
            elif resp.status_code != requests.codes.ok:
                raise HttpResponseError(f"Your query returned an HTTP status {resp.status_code}. The content returned from the request may be helpful:\n{resp.content.decode('utf-8')}")

            # Parse out pathway IDs and names
            pathway_dict = resp.json()
            names = []
            pathway_ids = []
            for pathway in pathway_dict:
                names.append(pathway["displayName"])
                pathway_ids.append(pathway["stId"])

            pathway_df = pd.DataFrame({"id": id, "pathway": names, "pathway_id": pathway_ids})
            pathway_df = pathway_df.sort_values(by="pathway_id")
            all_pathway_df = all_pathway_df.append(pathway_df)

    elif database.lower() == "wikipathways":

        path_here = os.path.abspath(os.path.dirname(__file__))
        data_dir_name = "data"
        file_name = "WikiPathwaysDataframe.tsv.gz"
        file_path = os.path.join(path_here, data_dir_name, file_name)
        df = pd.read_csv(file_path, sep="\t", index_col=0)
        all_pathway_df = pd.DataFrame()

        for protein in proteins:

            if protein in df.index:
                # Column headers are pathways; select pathways where the row for the protein has a
                # True for that pathway's column, indicating membership
                pathways = df.columns[df.loc[protein, :]].values

                prot_df = pd.DataFrame({"id": protein, "pathway": pathways})
                all_pathway_df = all_pathway_df.append(prot_df)

            else:
                if not quiet:
                    warnings.warn(f"The protein '{protein}' was not found in the WikiPathways data.", ParameterWarning, stacklevel=2)
    else:
        raise InvalidParameterError(f"Database '{database}' not recognized. Valid options: 'reactome', 'wikipathways'")

    all_pathway_df = all_pathway_df.reset_index(drop=True)
    return all_pathway_df

def get_proteins_in_pathways(pathways, database, quiet=False):
    """Query either the Reactome REST API or the downloaded WikiPathways dataframe to get a list of proteins contained in a particular pathway.

    Parameters:
    pathways (str or list of str): The pathway(s) to get the contained proteins for. If using Reactome, these must be pathway IDs (e.g. "R-HSA-140877").
    database (str): The database to use; either 'reactome' or 'wikipathways'.
    quiet (bool, optional): Whether to suppress warnings issued when identifiers are not found. Default False.

    Returns:
    pandas.DataFrame: The proteins contained in the pathways.
    """

    # Process string input
    if isinstance(pathways, str):
        pathways = [pathways]

    if database.lower() == "reactome":
        # Set headers and url
        headers = {"accept": "application/json"}

        # Loop over ids and get the interacting pathways
        all_protein_df = pd.DataFrame()
        for pathway_id in pathways:

            # Send the request
            url = f"https://reactome.org/ContentService/data/participants/{pathway_id}"
            resp = requests.get(url, headers=headers)

            if resp.status_code == 404 or (resp.status_code == requests.codes.ok and (len(resp.content.decode("utf-8")) == 0 or len(resp.json()) == 0)):
                if not quiet:
                    warnings.warn(f"The query for '{pathway_id}' found no results. You may have mistyped the pathway ID.", ParameterWarning, stacklevel=2)
                continue
            elif resp.status_code != requests.codes.ok:
                raise HttpResponseError(f"Your query returned an HTTP status {resp.status_code}. The content returned from the request may be helpful:\n{resp.content.decode('utf-8')}")

            # Parse all the proteins/genes out of the response
            members_df = pd.json_normalize(resp.json(), record_path=["refEntities"])
            prot_df = members_df[members_df["displayName"].str.startswith("UniProt:")]
            
            prot_names = prot_df["displayName"].str.rsplit(" ", n=1, expand=True)[1].\
                drop_duplicates(keep="first").\
                sort_values().\
                reset_index(drop=True)
            
            pathway_df = pd.DataFrame({"pathway": pathway_id, "member": prot_names})
            all_protein_df = all_protein_df.append(pathway_df)

        all_protein_df = all_protein_df.drop_duplicates(keep="first")

    elif database.lower() == "wikipathways":

        path_here = os.path.abspath(os.path.dirname(__file__))
        data_dir_name = "data"
        file_name = "WikiPathwaysDataframe.tsv.gz"
        file_path = os.path.join(path_here, data_dir_name, file_name)
        df = pd.read_csv(file_path, sep="\t", index_col=0)
        all_protein_df = pd.DataFrame()

        for pathway in pathways:

            if pathway in df.columns:
                prot_names = df.index[df[pathway]].values
                pathway_df = pd.DataFrame({"pathway": pathway, "member": prot_names})
                all_protein_df = all_protein_df.append(pathway_df)
                
            else:
                if not quiet:
                    warnings.warn(f"The pathway '{pathway}' was not found in the WikiPathways data.", ParameterWarning, stacklevel=2)

    else:
        raise InvalidParameterError(f"Database '{database}' not recognized. Valid options: 'reactome', 'wikipathways'")

    all_protein_df = all_protein_df.reset_index(drop=True)
    return all_protein_df

def reactome_pathway_overlay(pathway, df=None, analysis_token=None, open_browser=True, export_path=None, image_format="png", display_col_idx=0, diagram_colors="Modern", overlay_colors="Standard", quality=7):
    """Visualize numerical data (e.g. protein expression) on a Reactome pathway diagram, with each node's color corresponding to the expression value provided for that molecule.

    Parameters:
    pathway (str): The Reactome ID for the pathway you want to overlay the data on, e.g. "R-HSA-73929".
    df (pandas.DataFrame or pandas.Series, optional): If you haven't previously analyzed your data with Reactome, give this parameter the data you want to overlay. Each row corresponds to a particular gene/protein/etc, and each column is expression or other data for a sample or aggregate. Index must be unique identifiers. Multiple data columns allowed. All dtypes must be numeric. Default None assumes you are instead passing a token for previously analyzed data to the "analysis_token" parameter.
    analysis_token (str, optional): If the data you want to visualize has been recently analyzed in Reactome already, pass the token for that analysis to this parameter to overlay it on the specified pathway. This will allow this function to reaccess the archived results, thus avoiding wasting time by repeating the work of submitting and analyzing the data. Default of None assumes you are instead passing data to the "df" parameter.
    open_browser (bool, optional): Whether to automatically open the diagram in a new web browser tab. Default True.
    export_path (str, optional): A string providing a path to export the diagram to. Must end in a file name with the same extension as the "image_format" parameter. Default None causes no figure to be exported.
    image_format (str, optional): If export_path is not none, this specifies the format to export the diagram to. Options are "png", "gif", "svg", "jpg", "jpeg", or "pptx". Must match the file extension in the export path. If you're creating a gif and you want more than one column's data to be included in the image, make sure to pass None to the display_col_idx parameter. Default "png".
    display_col_idx (int, optional): If export_path is not none, this specifies which column in the dataframe to overlay expression data from. Must be a valid column index for the given table, or None. None will cause the first column to be displayed, unless you're creating a gif, in which case it will cause all columns to be included in the gif. Default None.
    diagram_colors (str, optional): If export_path is not none, this specifies the Reactome color scheme to use for the underlying diagram. Options are "Modern" or "Standard". Default "Modern".
    overlay_colors (str, optional): If export_path is not none, this specifies the Reactome color scheme to use for the data overlay. Options are "Standard", "Strosobar", or "Copper Plus". Default "Standard".
    quality (int, optional): If export_path is not none, this specifies what relative quality to export the image at. Must be between 1 and 10 inclusive. Default 7.

    Returns:
    list of float: The mean of the data values for all proteins in the pathways, for each column in the data table, in the order of the columns in the data table. I.e. each value in this list is the average of the data from a particular column for all the proteins in the pathway.
    str: If export_path is None, returns URL to diagram with data overlaid in Reactome Pathway Browser. Otherwise returns the path the image was exported to.
    """
    # Parameter checking
    if df is None and analysis_token is None:
        raise InvalidParameterError("You passed None to both the 'df' and 'analysis_token' parameters. You must pass a value to one of them.")

    if df is not None and analysis_token is not None:
        raise InvalidParameterError("You passed values to both the 'df' and 'analysis_token' parameters. You may only pass a value to one of them.")

    if export_path is not None:

        if image_format not in ("png", "gif", "svg", "jpg", "jpeg", "pptx"):
            raise InvalidParameterError(f"Invalid value for 'image_format' parameter. Valid options are 'png', 'gif', 'svg', 'jpg', 'jpeg', or 'pptx'. You passed '{image_format}'.")

        if display_col_idx is None:
            display_col_idx = ""        
        elif df is not None and display_col_idx not in range(0, df.shape[1] if isinstance(df, pd.DataFrame) else 1):
            raise InvalidParameterError(f"Invalid value for 'display_col_idx' parameter. Must be either None, or an int between 0 and one less than the number of columns in df (which is {df.shape[1] - 1} for this df), inclusive. You passed {display_col_idx}.")

        if diagram_colors not in ("Modern", "Standard"):
            raise InvalidParameterError(f"Invalid value for 'diagram_colors' parameter. Valid options are 'Modern' or 'Standard'. You passed '{diagram_colors}'.")

        if overlay_colors not in ("Standard", "Strosobar", "Copper Plus"):
            raise InvalidParameterError(f"Invalid value for 'overlay_colors' parameter. Valid options are 'Standard', 'Strosobar', or 'Copper Plus'. You passed '{overlay_colors}'.")

        if quality not in range(1, 11):
            raise InvalidParameterError(f"Invalid value for 'quality' parameter. Must be an int between 1 and 10 inclusive. You passed {quality}.")

        if image_format != export_path.split('.')[-1]:
            raise InvalidParameterError(f"The file extension in the 'export_path' parameter must match the 'image_format' parameter. For the image_format parameter, you passed '{image_format}'. The extension at the end of your export path was '{export_path.split('.')[-1]}'.")

        if export_path[:2] == "~/":
            raise InvalidParameterError("The export path you provided appeared to start with a reference to the user home directory. To avoid confusion, this function will not expand that reference. Please provide a full path instead.")

    if df is not None:

        df = df.copy(deep=True)

        # If they gave us a series, make it a dataframe
        if isinstance(df, pd.Series):
            if df.name is None:
                df.name = "data"
            df = pd.DataFrame(df)

        # Check that the index is strings (gene/protein names)
        if df.index.dtype != np.dtype("object"):
            raise InvalidParameterError(f"The dataframe you passed does not have an index of strings. The dtype of your dataframe's index is {df.index.dtype}.")

        elif df.index.str.isnumeric().any():
            raise InvalidParameterError(f"The dataframe you passed has an index of strings, but some of the values are completely numbers as strings. Gene/protein identifier always have letters in them.")

        # The identifier series (the index) needs to have a name starting with "#"
        if df.index.name is None:
            df.index.name = "#identifier"
        elif not df.index.name.startswith("#"):
            df.index.name = "#" + df.index.name

        # Take care of NaNs
        df = df.astype(str) # This represents NaNs as 'nan', which Reactome is OK with

        # Get the df as a tab-separated string
        df_str = df.to_csv(sep='\t')

        # Post the data to the Reactome analysis service
        analysis_url = "https://reactome.org/AnalysisService/identifiers/projection"
        headers = {"Content-Type": "text/plain"}
        params = {"pageSize": "0", "page": "1"} # We only need the analysis token, so set pageSize to 0 so we don't worry about getting any of the data for individual pathways.

        view_resp = requests.post(analysis_url, headers=headers, params=params, data=df_str)

        # Check that the response came back good
        if view_resp.status_code != requests.codes.ok:
            raise HttpResponseError(f"Submitting your data for analysis returned an HTTP status {view_resp.status_code}. The content returned from the request may be helpful:\n{view_resp.content.decode('utf-8')}")    

        # Get the token for accessing the analysis results
        token = view_resp.json()["summary"]["token"]

    else:
        token = analysis_token

    # Get the mean data values
    expr_url = f"https://reactome.org/AnalysisService/token/{token}/filter/pathways?resource=TOTAL&pValue=1"

    headers = {
        "accept": "application/json",
        "content-type": "text/plain",
    }

    expr_resp = requests.post(expr_url, headers=headers, data=pathway)

    # Check that the response came back good
    if expr_resp.status_code != requests.codes.ok:
        raise HttpResponseError(f"Submitting your data for analysis returned an HTTP status {expr_resp.status_code}. The content returned from the request may be helpful:\n{expr_resp.content.decode('utf-8')}")    

    # Get the expression list
    expr_json = expr_resp.json()

    if len(expr_json) > 0:
        expr_list = expr_resp.json()[0]["entities"]["exp"]
    else:
        expr_list = []

    # Use the token and the pathway ID to open the pathway diagram with the data overlaid in the Reactome Pathway Browser
    viewer_url = f"https://reactome.org/PathwayBrowser/#/{pathway}&DTAB=AN&ANALYSIS={token}"
    if open_browser:
        webbrowser.open(viewer_url)

    if export_path is not None:

        # Get the diagram
        export_url = f"https://reactome.org/ContentService/exporter/diagram/{pathway}.{image_format}?token={token}&resource=TOTAL&diagramProfile={diagram_colors}&analysisProfile={overlay_colors}&expColumn={display_col_idx}&quality={quality}"
        export_resp = requests.get(export_url)

        # Check that the response came back good
        if export_resp.status_code != requests.codes.ok:
            raise HttpResponseError(f"Submitting your data for analysis returned an HTTP status {export_resp.status_code}. The content returned from the request may be helpful:\n{export_resp.content.decode('utf-8')}")    

        # Save the image
        with open(export_path, 'wb') as dest:
            dest.write(export_resp.content)

    if export_path is None:
        return expr_list, viewer_url
    else:
        return expr_list, export_path

def reactome_enrichment_analysis(analysis_type, data, sort_by, ascending, include_high_level_diagrams=True, disease_pathways=True, include_interactors=False):
    """Use the Reactome Analysis Service API to do a gene set enrichment analysis.

    Parameters:
    analysis_type (str): The type of enrichment analysis you want to perform. Either "ranked" or "nonranked".
    data (pandas.DataFrame or pandas.Series, or list or array-like): The data you want to overlay. Format depends on the analysis type.
        If ranked analysis:
            - data is a DataFrame or Series where the index is unique gene/protein identifiers and column(s) are ranking values (e.g. expression values for genes).
            - Multiple data columns allowed and are analyzed as separate ranked enrichment analyses.
            - All dtypes must be numeric.
        If unranked analysis:
            - data is a list or pandas.Index of identifiers to test pathways for enrichment with.
    sort_by (str): The metric by which to sort the pathways when selecting the top ones. You can pass "p_value" to sort by the P value (hypergeometric distribution), or pass one of the metrics directly supported by the Reactome API, listed below. (Yes, our function just maps "p_value" to "ENTITIES_PVALUE" .)
        "NAME",
        "TOTAL_ENTITIES",
        "TOTAL_INTERACTORS",
        "TOTAL_REACTIONS",
        "FOUND_ENTITIES",
        "FOUND_INTERACTORS",
        "FOUND_REACTIONS",
        "ENTITIES_RATIO", # Note: This value is the ratio of the total number of proteins in the pathway, to the total number of proteins in all of Reactome.
        "ENTITIES_PVALUE",
        "ENTITIES_FDR",
        "REACTIONS_RATIO",
    ascending (bool): When sorting pathways by the specified metric, whether to put smaller values first.
    include_high_level_diagrams (bool, optional): Whether to include pathway diagrams in the output that are just summaries of lower level pathways and don't show actual reactions. If False, this will exclude any Reactome pathways that have EHLD (Enhanced Higher Level Diagrams). Default True includes these pathways in results.
    disease_pathways (bool, optional): Whether to include pathways that describe disease related function. Default True.
    include_interactors (bool, optional): Whether to include computationally inferred interactors when identifying pathways that are enriched with your submitted proteins/genes. Default False. You may want to set this to True if a large portion of the identifiers you submitted do not match a Reactome pathway when it is set to False.

    Returns:
    pandas.DataFrame: A dataframe with info on the top enriched pathways, sorted by the specified metric.
    """
    # Check the sort_by parameter
    if sort_by == "p_value":
        parsed_sort_by = "ENTITIES_PVALUE"
    else:
        parsed_sort_by = sort_by

    valid_sort_bys = [
        "NAME",
        "TOTAL_ENTITIES",
        "TOTAL_INTERACTORS",
        "TOTAL_REACTIONS",
        "FOUND_ENTITIES",
        "FOUND_INTERACTORS",
        "FOUND_REACTIONS",
        "ENTITIES_RATIO",
        "ENTITIES_PVALUE",
        "ENTITIES_FDR",
        "REACTIONS_RATIO"]
    
    if parsed_sort_by not in valid_sort_bys: 
        newline = "\n"
        single_qt = "'"
        raise InvalidParameterError(f"Invalid value for 'sort_by' parameter. You passed '{sort_by}'. Must be one of the following:\n'p_value'\n{newline.join([f'{single_qt}{x}{single_qt}' for x in valid_sort_bys])}")

    if analysis_type == "ranked":

        # Copy the data
        data = data.copy(deep=True)
        
        # If they gave us a series, make it a dataframe
        if isinstance(data, pd.Series):
            if data.name is None:
                data.name = "data"
            data = pd.DataFrame(data)

        # Check that the index is strings (gene/protein names)
        if data.index.dtype != np.dtype("object"):
            raise InvalidParameterError(f"The dataframe you passed does not have an index of strings. The dtype of your dataframe's index is {data.index.dtype}.")

        elif data.index.str.isnumeric().any():
            raise InvalidParameterError(f"The dataframe you passed has an index of strings, but some of the values are completely numbers as strings. Gene/protein identifier always have letters in them.")

        # The identifier series (the index) needs to have a name starting with "#"
        if data.index.name is None:
            data.index.name = "#identifier"
        elif not data.index.name.startswith("#"):
            data.index.name = "#" + data.index.name

        # Take care of NaNs and small numbers
        # This represents NaNs as 'nan', which Reactome is OK with
        # Also rounds all numbers without scientific notation
        data.iloc[:, 0] = data.iloc[:, 0].apply(lambda x: "{:.9f}".format(x)) 

        # Get the data as a tab-separated string
        data_str = data.to_csv(sep='\t')

    elif analysis_type == "unranked":

        # Format data
        data = pd.Index(data) # Convert it to an index if it wasn't
        data = data.dropna() # Drop NaNs
        data = data.astype(str) # Make it strings

        # Check that they're actual gene/protein names
        if data.str.isnumeric().any():
            raise InvalidParameterError(f"The data you passed has some values that are completely numeric. Gene/protein identifier always have letters in them.")

        # The first item needs to be a column header string starting with '#'
        if not data[0].startswith("#"):
            data = data.insert(0, "#identifier")

        # Get the list as a newline-separated string
        data_str = "\n".join(data)

    else:
        raise InvalidParameterError(f"Invalid value for 'analysis_type' parameter. You passed '{analysis_type}'. Must be 'ranked' or 'unranked'.")

    # Post the data to the Reactome analysis service
    analysis_url = "https://reactome.org/AnalysisService/identifiers/projection"
    headers = {"Content-Type": "text/plain"}
    params = {
        "interactors": include_interactors,
        "sortBy": parsed_sort_by,
        "order": "ASC" if ascending else "DESC",
        "includeDisease": disease_pathways,
    }

    resp = requests.post(analysis_url, headers=headers, params=params, data=data_str)

    # Check that the response came back good
    if resp.status_code != requests.codes.ok:
        raise HttpResponseError(f"Submitting your data for analysis returned an HTTP status {resp.status_code}. The content returned from the request may be helpful:\n{resp.content.decode('utf-8')}")

    warnings_list = resp.json()["warnings"]
    if len(warnings_list) != 0:
        newline = "\n"
        raise InvalidParameterError(f"Your analysis request returned the following warnings. You may have a data formatting problem. Check that your data matches the format specified in the docstring. Here's up to the first ten warnings:\n{newline.join(warnings_list[0:len(warnings_list)] if len(warnings_list) < 10 else warnings_list[0:10] + ['...'])}")

    # Process the JSON response
    resp_json = resp.json()
    analysis_token = resp_json["summary"]["token"]
    pathways_table = pd.json_normalize(resp_json["pathways"], sep="_")
    
    # Select the columns we want
    pathways_table = pathways_table[["stId", "name", "entities_ratio", "entities_pValue", "entities_fdr", "entities_found", "entities_total"]]

    # If requested, filter out EHLD pathways
    if not include_high_level_diagrams:

        # Download a list of all diagrams with EHLD diagrams
        ehld_url = "https://reactome.org/download/current/ehld/svgsummary.txt"
        ehld_resp = requests.get(ehld_url)

        # Check that the response came back good
        if ehld_resp.status_code != requests.codes.ok:
            raise HttpResponseError(f"Checking whether pathways are high level pathways returned an HTTP status {ehld_resp.status_code}. The content returned from the request may be helpful:\n{ehld_resp.content.decode('utf-8')}")

        # Parse the response
        ehld_list = ehld_resp.content.decode("utf-8")
        ehld_list = ehld_list.split("\n")
        ehld_list = [pathway_id for pathway_id in ehld_list if pathway_id.startswith("R-HSA-")]

        has_ehld = pathways_table["stId"].isin(ehld_list)
        pathways_table = pathways_table[~has_ehld]

        # Make the index look normal
        pathways_table = pathways_table.reset_index(drop=True)

    return analysis_token, pathways_table
