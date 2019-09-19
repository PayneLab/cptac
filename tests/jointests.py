import cptac
import pytest

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
		ds2 = self.datasets
		for dataset1 in ds1:
			for dataset2 in ds2:
				valid_omics = set(dataset1._valid_omics_dfs + dataset2._valid_omics_dfs)
				for omic in valid_omics:
					ds1.join_omics_to_omics=

				#if omic1 in dataset._data.keys() and omic2 in dataset.data.keys():
					#do the join


t = JoinTest()
t.testOmicsToOmics()
