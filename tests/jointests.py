from typing import Union

import cptac
import pytest
import sys
import traceback
import random

from cptac import Ovarian, Ccrcc, Colon, Endometrial


class InvalidArgumentsError(Exception):
    def __init__(self, message):
        self.message = message


class JoinTest:
    def __init__(self):
#        cptac.download(dataset="brca", version='latest')
#        cptac.download(dataset="ccrcc", version='latest')
#        cptac.download(dataset="colon", version='latest')
#        cptac.download(dataset="endometrial", version='latest')
#        cptac.download(dataset="gbm", version='latest')
#        cptac.download(dataset="hsncc", version='latest')
#        cptac.download(dataset="luad", version='latest')
#        cptac.download(dataset="ovarian", version='latest')

        self.brca = cptac.Brca("3.1.1")
        self.ccrcc = cptac.Ccrcc("0.1.1")
        self.colon = cptac.Colon("0.0.1")
        self.en = cptac.Endometrial("2.1.1")
        self.gbm = cptac.Gbm("3.0")
        self.hsncc = cptac.Hnscc("2.0")
        self.luad= cptac.Luad("3.1")
        self.ovarian = cptac.Ovarian("0.0.1")

        self.datasets = [
            self.brca,
            self.ccrcc,
            self.colon,
            self.en,
            self.gbm,
            self.hsncc,
            self.luad,
            self.ovarian,
            ]

    @staticmethod
    def RepresentsInt(s: str) -> bool:
        """tests to see if the given string can be parsed into an integer"""
        try:
            int(s)
            return True
        except ValueError:
            return False

    @staticmethod
    def usage() -> object:
        """Gets a string that describes how to use this class"""
        return "USAGE: ['-omics, -metadata, -mutations'] mutation names (if joining mutations)"

    def parseArguments(self, system_arguments: list):
        """
        Parses the given system arguements taken from the main function and then executes the appropriate function
        @param system_arguments any of the system arguments in -mutations, -metadata, -omics where mutations is
        followed by a list of mutations
        @raises Exception if there are incorrect parameters
        """
        if len(system_arguments) == 0:
            raise InvalidArgumentsError(self.usage())
        try:
            if system_arguments[0] == "-mutations":  # wrong order, mutations must be attached to metadata or omics
                raise Exception("mutations must come after either metadata or omics")

            elif system_arguments[0] == "-metadata":  # join metadata to something
                if system_arguments[1] == "-mutations":  # join metadata to mutations
                    if self.RepresentsInt(system_arguments[2]):  # int x amount of random mutations
                        amt_random_mutations = int(system_arguments[2])
                        self.testMetaDataToMutations(True, amt_random_mutations)
                    else:  # specific mutations listed after the -mutations argument
                        mutations = system_arguments[2:]
                        self.testMetaDataToMutations(False, mutations)
                elif system_arguments[1] == "-omics":  # join metadata to omics
                    self.testMetaDataToOmics()
                elif system_arguments[1] == "-metadata":  # join metadata to metadata
                    self.testMetaDataToMetaData()
                else:  # invalid arguments
                    raise InvalidArgumentsError(self.usage())

            elif system_arguments[0] == "-omics":  # join omics to something
                if system_arguments[1] == "-mutations":  # join omics to mutations
                    if self.RepresentsInt(system_arguments[2]):  # try int x amount of random mutations
                        amt_random_mutations = int(system_arguments[2])
                        is_random_join = True
                        self.testOmicsToMutations(is_random_join, amt_random_mutations)
                    else:  # specific mutations listed after the -mutations argument
                        mutations = system_arguments[2:]
                        is_random_join = False
                        self.testOmicsToMutations(is_random_join, mutations)
                elif system_arguments[1] == "-metadata":  # join omics to metadata
                    self.testMetaDataToOmics()
                elif system_arguments[1] == "-omics":  # join omics to omics
                    self.testOmicsToOmics()
                else:  # invalid arguments
                    raise InvalidArgumentsError(self.usage())
            else:
                raise InvalidArgumentsError(self.usage())
        except InvalidArgumentsError as i:
            print(i.message)
            raise Exception(self.usage())

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
                            print(
                                f"In datasets: {dataset1.get_cancer_type()}, {omic} did not successfully join with {omic2}. \n\n")
                            traceback.print_exc()
                            sys.exit(0)
        assert (not error)

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
                        print(
                            f"In datasets: {dataset1.get_cancer_type()}, {omic} did not successfully join with {metadata}. \n\n")
                        traceback.print_exc()
                        sys.exit(0)
        assert (not error)

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
                        print(
                            f"In datasets: {dataset1.get_cancer_type()}, {md1} did not successfully join with {md2}. \n\n")
                        traceback.print_exc()
                        sys.exit(0)
        assert (not error)

    def testOmicsToMutations(self, is_random_join: bool, *mutations):
        ds1 = self.datasets
        error = False
        if not is_random_join and mutations is None:
            for dataset1 in ds1:
                if dataset1.get_cancer_type() == 'brca' or dataset1.get_cancer_type() == 'luad':
                    continue
                valid_omics = set(dataset1._valid_omics_dfs).intersection(set(dataset1._data.keys()))
                valid_mutations = set(dataset1.get_somatic_mutation())
                for omic in valid_omics:
                    for mutation in valid_mutations:
                        print(f"joining {omic}, {mutation} in dataset {dataset1.get_cancer_type()}\n\n")
                        try:
                            cross = dataset1.join_omics_to_mutations(omic, mutation)
                        except Exception as e:
                            error = True
                            print(f"Error: {e}\n\n")
                            print(
                                f"In datasets: {dataset1.get_cancer_type()}, {omic} did not successfully join with {mutation}. \n\n")
                            traceback.print_exc()

        elif is_random_join and isinstance(mutations, int):
            amt_random_mutations: int = mutations
            dataset1: Union[Endometrial, Ovarian, Ccrcc, Colon]
            for dataset1 in ds1:
                if dataset1.get_cancer_type() == 'brca' or dataset1.get_cancer_type() == 'luad':
                    continue
                valid_omics = set(dataset1._valid_omics_dfs).intersection(dataset1._data.keys())
                valid_mutations = set(dataset1.get_somatic_mutation())
                random_mutations = random.choices(valid_mutations, k=amt_random_mutations)
                for omic in valid_omics:
                    for mutation in random_mutations:
                        print(f"joining {omic}, {mutation} in dataset {dataset1.get_cancer_type()}\n\n")
                        try:
                            cross = dataset1.join_omics_to_mutations(omic, mutation)
                        except Exception as e:
                            error = True
                            print(f"Error: {e}\n\n")
                            print(
                                f"In datasets: {dataset1.get_cancer_type()}, {omic} did not successfully join with {mutation}. \n\n")
                            traceback.print_exc()

        elif not is_random_join and isinstance(mutations, list):
            mutations_list: list = mutations
            for dataset1 in ds1:
                if dataset1.get_cancer_type() == 'brca' or dataset1.get_cancer_type() == 'luad':
                    continue
                valid_omics = set(dataset1._valid_omics_dfs).intersection(dataset1._data.keys())
                valid_mutations = set(dataset1.get_somatic_mutation())
                for omic in valid_omics:
                    for mutation in mutations_list:
                        print(f"joining {omic}, {mutation} in dataset {dataset1.get_cancer_type()}\n\n")
                        try:
                            cross = dataset1.join_omics_to_mutations(omic, mutation)
                        except Exception as e:
                            error = True
                            print(f"Error: {e}\n\n")
                            print(
                                f"In datasets: {dataset1.get_cancer_type()}, {omic} did not successfully join with {mutation}. \n\n")
                            traceback.print_exc()

        else:
            error = True
            print("mutations must be an instance of list or followed by an integer")

        assert (not error)
        if not error:
            print(f"All joins for OmicsToMutations were successful")

    def testMetaDataToMutations(self, is_random_join: bool, *mutations):
        ds1 = self.datasets
        error = False

        if not is_random_join and mutations is None:
            for dataset1 in ds1:
                if dataset1.get_cancer_type() == 'brca' or dataset1.get_cancer_type() == 'luad':
                    continue
                valid_metadata = set(dataset1._valid_metadata_dfs).intersection(set(dataset1._data.keys()))
                valid_mutations = dataset1.get_somatic_mutation()
                for metadata in valid_metadata:
                    for mutation in valid_mutations["Gene"].drop_duplicates():
                        print(f"joining {metadata}, {mutation} in dataset {dataset1.get_cancer_type()}\n\n")
                        try:
                            cross = dataset1.join_metadata_to_mutations(metadata, mutation)
                        except Exception as e:
                            error = True
                            print(f"Error: {e}\n\n")
                            print(
                                f"In datasets: {dataset1.get_cancer_type()}, {metadata} did not successfully join with {mutation}. \n\n")
                        traceback.print_exc()
                        sys.exit(0)

        elif not is_random_join and isinstance(mutations, list):
            for dataset1 in ds1:
                if dataset1.get_cancer_type() == 'brca' or dataset1.get_cancer_type() == 'luad':
                    continue
                valid_metadata = set(dataset1._valid_metadata_dfs).intersection(set(dataset1._data.keys()))
                valid_mutations = mutations
                for metadata in valid_metadata:
                    for mutation in valid_mutations["Gene"].drop_duplicates():
                        print(f"joining {metadata}, {mutation} in dataset {dataset1.get_cancer_type()}\n\n")
                        try:
                            cross = dataset1.join_metadata_to_mutations(metadata, mutation)
                        except Exception as e:
                            error = True
                            print(f"Error: {e}\n\n")
                            print(
                                f"In datasets: {dataset1.get_cancer_type()}, {metadata} did not successfully join with {mutation}. \n\n")

        elif is_random_join and isinstance(mutations, int):
            amt_random_mutations = mutations
            for dataset1 in ds1:
                if dataset1.get_cancer_type() == 'brca' or dataset1.get_cancer_type() == 'luad':
                    continue
                valid_metadata = set(dataset1._valid_metadata_dfs).intersection(set(dataset1._data.keys()))
                valid_mutations = list(dataset1.get_somatic_mutation()["Gene"].drop_duplicates())
                random_mutations = random.choices(valid_mutations, k=amt_random_mutations)
                for metadata in valid_metadata:
                    for mutation in valid_mutations["Gene"].drop_duplicates():
                        print(f"joining {metadata}, {mutation} in dataset {dataset1.get_cancer_type()}\n\n")
                        try:
                            cross = dataset1.join_metadata_to_mutations(metadata, mutation)
                        except Exception as e:
                            error = True
                            print(f"Error: {e}\n\n")
                            print(
                                f"In datasets: {dataset1.get_cancer_type()}, {metadata} did not successfully join with {mutation}. \n\n")
        assert (not error)
        if not error:
            print(f"All joins for Metadata to Mutations were successful")



if __name__ == '__main__':
    arguments = sys.argv[1:]
    t = JoinTest()
    t.parseArguments(arguments)
