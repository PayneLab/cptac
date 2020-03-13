import sys
import pandas as pd



"""
@Return:
	A list of proteins known by the most recent WikiPathways download to be interacting parters with the specified protein.
	Returns None if specified protein is not found in the WikiPathways dataframe (which was intersected with Uniprot).

This function takes a path to WikiPathways Dataframe file and protein name and returns a list of all the proteins that interact with it, using the pathways from the WikiPathways relsease file.
This function loads the WikiPathways dataframe, and iterates through the row labelled with that protein name, return every protein in a pathway that also contains that protein.
"""

def get_interacting_proteins_wikipathways():
	#get variables from command line
		#sys.argv[1]: path to WikiPathways Dataframe .csv file
		#sys.argv[2]: protein name that the interacting proteins are based on

	WikiPathwaysDataframePath = sys.argv[1]
	proteinName = sys.argv[2]

	df = pd.read_csv(WikiPathwaysDataframePath, index_col=0) #NEED PATH TO file
	row = df.loc[proteinName]
	filtered_df = df.loc[:, row.values.tolist()]
	def has_true(values):
		for val in values:
			if val == True:
				return True
		return False
	filtered_df_final = filtered_df.loc[filtered_df.apply(lambda row: has_true(row.values.tolist()), axis=1), :]
	return filtered_df_final.index.tolist()


print(get_interacting_proteins_wikipathways())
