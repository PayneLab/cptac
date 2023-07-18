import cptac

def get_cancer_class(cancer_str):
    """
    Converts a string to a corresponding cancer class.

    Args:
    cancer_str (str): A string identifying the cancer type. This should match one of the following: 
        'brca', 'ccrcc', 'coad', 'gbm', 'hnscc', 'lscc', 'luad', 'ov', 'pdac', 'ucec', 'all_cancers'.
    
    Returns:
    A cancer class from the cptac module.

    Raises:
    ValueError: If the provided string does not correspond to a known cancer class.
    """
    # This dictionary should be updated as necessary
    mapping = {
        "brca": cptac.Brca,
        "ccrcc": cptac.Ccrcc,
        "coad": cptac.Coad,
        "gbm": cptac.Gbm,
        "hnscc": cptac.Hnscc,
        "lscc": cptac.Lscc,
        "luad": cptac.Luad,
        "ov": cptac.Ov,
        "pdac": cptac.Pdac,
        "ucec": cptac.Ucec,
        "all_cancers" : cptac.Ucec
    }

    try:
        return mapping[cancer_str.lower()]
    except KeyError:
        raise ValueError(f"'{cancer_str}' is not a known cancer class. Valid options are: {list(mapping.keys())}")

