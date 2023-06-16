import pytest
import cptac

def test_multi_join():
    """Test that multi_join correctly joins multiple dataframes"""
    # Instantiate a cancer object
    brca = cptac.Brca()

    # Prepare a join dictionary
    join_dict = {
        ("umich", "proteomics"): ["TP53", "EGFR"],
        ("washu", "transcriptomics"): ["TP53", "EGFR"],
        ("mssm", "clinical"): ["age", "sex"]
    }

    clinical_df = brca.get_dataframe('clinical', 'mssm')
    print(clinical_df.columns)

    # Use multi_join to join them
    result = brca.multi_join(join_dict, how='inner', flatten=True)
    print(result.columns)

    # Assert the resulting DataFrame is not empty
    assert not result.empty, "The join operation resulted in an empty DataFrame."

    # Assert the joined DataFrame has the expected columns
    expected_columns = ["TP53_umich_proteomics", "EGFR_umich_proteomics", "TP53_washu_transcriptomics", "EGFR_washu_transcriptomics", "age", "sex"]
    assert all(col in result.columns for col in expected_columns), "Joined DataFrame does not have the expected columns."

    # Check the number of rows and columns in the resulting DataFrame
    assert result.shape[0] > 0, "The joined DataFrame should have more than 0 rows."
    assert result.shape[1] == len(expected_columns), "The joined DataFrame does not have the expected number of columns."

    # Add more assertions to check the contents of the DataFrame if needed