# import itertools
# import logging
# import pytest
# import cptac
# from collections import defaultdict
# from .conftest import get_cancer_inputs

# LOGGER = logging.getLogger(__name__)

# def get_join_inputs():
#     join_inputs = defaultdict(lambda: defaultdict(list))

#     for cancer, source, dtype in get_cancer_inputs():
#         if dtype.startswith("acetylproteomics"):
#             join_inputs[cancer]['omics'].append((source, dtype))
#         elif dtype.startswith("clinical"):
#             join_inputs[cancer]['metadata'].append((source, dtype))
#         else:
#             join_inputs[cancer]['multi'].append((source, dtype))

#     return join_inputs

# @pytest.mark.parametrize("cancer_name, cat1", get_join_inputs().items())
# def test_join_omics_to_omics(cancer_name, cat1):
#     _run_combos(cancer_name, cat1)

# @pytest.mark.parametrize("cancer_name, cat1, cat2", itertools.product(get_join_inputs().keys(), ['omics'], ['metadata', 'multi']))
# def test_join_omics_to_metadata(cancer_name, cat1, cat2):
#     _run_combos(cancer_name, cat1, cat2)

# def _run_combos(cancer_name, cat1, cat2=None):
#     # Arrange
#     cancer_obj = getattr(cptac, cancer_name.title())()
#     cancer_combos = get_join_inputs()[cancer_name]
    
#     if cat2 is None:
#         test_units = itertools.combinations(cancer_combos[cat1], 2)
#         join_func = cancer_obj.join_omics_to_omics
#     else:
#         test_units = itertools.product(cancer_combos[cat1], cancer_combos[cat2])
#         join_func = cancer_obj.join_omics_to_metadata

#     # Act & Assert
#     for source1, source2 in test_units:
#         if source1[0] == 'broad' or source2[0] == 'broad':
#             LOGGER.debug(f"one of the datasets in ({source1}, {source2}) is too large")
#             continue

#         LOGGER.debug(f"testing join of {source1} to {source2}")

#         try:
#             joined_data = join_func(
#                 df1_name=source1[1],
#                 df2_name=source2[1],
#                 df1_source=source1[0],
#                 df2_source=source2[0]
#             )
#             assert isinstance(joined_data, pd.DataFrame), f"Join operation on {source1[1]} and {source2[1]} did not return a DataFrame"

#         except Exception as e:
#             pytest.fail(f"Failed to perform join on {source1[1]} and {source2[1]} for {cancer_name}. Error: {str(e)}")