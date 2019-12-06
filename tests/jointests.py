import cptac
import pytest
import sys
import traceback

class JoinTest:
	def __init__(self):
		cptac.download(dataset="endometrial", version='latest')
		#cptac.download(dataset="brca", version='latest')
		#cptac.download(dataset="gbm", version='latest')
		#cptac.download(dataset="hsncc", version='latest')
		#cptac.download(dataset="luad", version='latest')
		cptac.download(dataset="ovarian", version='latest')
		cptac.download(dataset="ccrcc", version='latest')
		cptac.download(dataset="colon", version='latest')
		self.en = cptac.Endometrial()
		#self.brca = cptac.Brca()
		#self.gbm = cptac.Gbm()
		#self.hsncc = cptac.Hnscc()
		#self.luad= cptac.Luad()
		self.ovarian = cptac.Ovarian()
		self.ccrcc = cptac.Ccrcc()
		self.colon = cptac.Colon()
		#self.datasets = list(self.en,self.brca,self.gbm,self.hsncc,self.luad,self.ovarian,self.ccrcc)
		self.datasets = list([self.en,self.ovarian,self.ccrcc,self.colon])
	

	def testOmicsToOmics(self):
		ds1 = self.datasets
		error = False
		for dataset1 in ds1:
			valid_omics = set(dataset1._valid_omics_dfs).intersection(set(dataset1._data.keys()))
			for omic in valid_omics:
				for omic2 in valid_omics:
					if omic == omic2:
						continue
					else:
						print(f"\n\njoining {omic}, {omic2} in dataset {dataset1.get_cancer_type()}\n\n")
						try:
							cross = dataset1.join_omics_to_omics(df1_name=omic, df2_name=omic2)
						except Exception as e:
							print(f"Error: {e}\n\n")
							print(f"In datasets: {dataset1.get_cancer_type()}, {omic} did not successfully join with {omic2}. \n\n")
							traceback.print_exc()
							sys.exit(0)
		assert(not error)
	
	def testMetaDataToOmics(self):
		ds1 = self.datasets
		error = False
		for dataset1 in ds1:
			valid_omics = set(dataset1._valid_omics_dfs).intersection(set(dataset1._data.keys()))
			valid_metadata = set(dataset1._valid_metadata_dfs).intersection(set(dataset1._data.keys()))
			for omic in valid_omics:
				for metadata in valid_metadata:
					print(f"joining {omic}, {metadata} in dataset {dataset1.get_cancer_type()}\n\n")
					try:
						cross = dataset1.join_metadata_to_omics(metadata_df_name=metadata, omics_df_name=omic)
					except Exception as e:
						error = True
						print(f"Error: {e}\n\n")
						print(f"In datasets: {dataset1.get_cancer_type()}, {omic} did not successfully join with {metadata}. \n\n")
						traceback.print_exc()
						sys.exit(0)
		assert(not error)

	def testMetaDataToMetaData(self):
		ds1 = self.datasets
		error = False
		for dataset1 in ds1:
			valid_metadata = set(dataset1._valid_metadata_dfs).intersection(set(dataset1._data.keys()))
			for md1 in valid_metadata:
				for md2 in valid_metadata:
					if md1 == md2:
						continue
					print(f"joining {md1}, {md2} in dataset {dataset1.get_cancer_type()}\n\n")
					try:
						cross = dataset1.join_metadata_to_metadata(df1_name=md1, df2_name=md2)
					except Exception as e:
						error = True
						print(f"Error: {e}\n\n")
						print(f"In datasets: {dataset1.get_cancer_type()}, {md1} did not successfully join with {md2}. \n\n")
						traceback.print_exc()
						sys.exit(0)
		assert(not error)

	def testOmicsToMutations(self, specific_mutations = None):
		ds1 = self.datasets
		error = False
		for dataset1 in ds1:
			if dataset1.get_cancer_type() == 'brca' or dataset1.get_cancer_type() == 'luad':
				continue
			valid_omics = set(dataset1._valid_omics_dfs).intersection(set(dataset1._data.keys()))
			valid_mutations= set(dataset1.get_somatic_mutation())
			if specific_mutations:
				for mutation in specific_mutations:
					if mutation not in valid_mutations:
						print(self.usage())
						return
				valid_mutations = specific_mutations
			for omic in valid_omics:
				for mutation in valid_mutations:
					print(f"joining {omic}, {mutation} in dataset {dataset1.get_cancer_type()}\n\n")
					try:
						cross = dataset1.join_omics_to_mutations(omics_df_name=omic, mutations_genes=mutation)
					except Exception as e:
						error = True
						print(f"Error: {e}\n\n")
						print(f"In datasets: {dataset1.get_cancer_type()}, {omic} did not successfully join with {mutation}. \n\n")
						traceback.print_exc()
						sys.exit(0)
		assert(not error)
		if not error:
			print(f"All joins for OmicsToMutations were successful")

	def testMetaDataToMutations(self):
		ds1 = self.datasets
		error = False
		for dataset1 in ds1:
			if dataset1.get_cancer_type() == 'brca' or dataset1.get_cancer_type() == 'luad':
				continue
			valid_metadata = set(dataset1._valid_metadata_dfs).intersection(set(dataset1._data.keys()))
			valid_mutations = dataset1.get_somatic_mutation()
			for metadata in valid_metadata:
				for mutation in valid_mutations["Gene"].drop_duplicates():
					print(f"joining {metadata}, {mutation} in dataset {dataset1.get_cancer_type()}\n\n")
					try:
						cross = dataset1.join_metadata_to_mutations(metadata_df_name=metadata, mutations_genes=mutation)
					except Exception as e:
						error = True
						print(f"Error: {e}\n\n")
						print(f"In datasets: {dataset1.get_cancer_type()}, {metadata} did not successfully join with {mutation}. \n\n")
						traceback.print_exc()
						sys.exit(0)
		assert(not error)
		if not error:
			print(f"All joins for Metadata to Mutations were successful")

	def testLostData(self):
		#test for lost data in the columns
		return

	def usage(self):
		return f"\nUSAGE: python jointests.py (-metadata, -omics, -mutations [genes])\n"
	
	def handleMutations(args):
		if args[0][0] == "-":
			if args[0] == "-metadata":
				self.joinMetadataToMutations()
			elif args[0] == "-omics":
				self.testOmicsToMutations()
			else:
				print(self.usage())
		else:
			if args[-1] not in ["-omics", "-metadata"]:
				print(f"\nERROR: last paramter must be any of -omics or -metadata" + self.usage())
			else:
				mutations = args[:-1]


	def main(args):
		args = []
		key_args = ["-mutations", "-metadata", "-omics"]
		t = joinTest()
		if !args:
			t.testOmicsToOmics()
			t.testMetaDataToMutations()
			t.testOmicstoMutations()
		else:
			if args[0] not in key_args:
				print(t.usage())
			else:
				if args[0] == "-mutations":
					t.handleMutations(args[1:])

							
					


if __name__ == "__main__":
	main(sys.argv[1:])
