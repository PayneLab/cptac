import pytest
import cptac
import sys
from itertools import product

@pytest.fixture(scope="session")
def cancer_source_dtype_combos():
    # get options df from cptac
    options = cptac.get_options()
    options.drop(["Loadable datatypes"],  axis=1, inplace=True)
    options.set_index("Cancers")

    # convert Datatypes from type string to type list
    options.Datatypes = options.Datatypes.apply(lambda x: x.split(", "))

    # make a list of [(cancer, {source:datatype})] items
    combos = options.apply(lambda x: [(x[0], {source: datatype}) for (source, datatype) in list(product([x[1]], x[2]))], axis=1)
    
    # return a flattened list of (cancer, {source:datatype}) items
    return [item for sublist in combos.to_list() for item in sublist]