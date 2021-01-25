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
import requests
import cptac.exceptions as ex

def _pdc_download_cancer_type(cancer_type, version, redownload):
    """Download data for the specified cancer type from the PDC."""

    # Steps
    # 1. Use GraphQL queries to download JSON
    # 2. Convert to tables
    # 3. Save to tsv.gz
    # 4. Create index files.
    #     - Leave hash and url fields blank, but do mark with \t characters. They're only used in download functions.

    study_submitter_id_map = {
        "brca": {
            "phosphoproteome_bi": "S039-2",
            "proteome_bi": "S039-1",
        },
        "ccrcc": {
            "phosphoproteome": "S044-2",
            "proteome": "S044-1",
        },
        "colon": {
            "phosphoproteome_pnnl": "S037-3",
            "proteome_pnnl": "S037-2",
        },
        "gbm": {
            "acetylome": "CPTAC GBM Discovery Study - Acetylome",
            "phosphoproteome": "CPTAC GBM Discovery Study - Phosphoproteome",
            "proteome": "CPTAC GBM Discovery Study - Proteome",
        },
        "hnscc": {
            "phosphoproteome": "CPTAC HNSCC Discovery Study - Phosphoproteome",
            "proteome": "CPTAC HNSCC Discovery Study - Proteome",
        },
        "lscc": {
            "acetylome": "CPTAC LSCC Discovery Study - Acetylome",
            "phosphoproteome": "CPTAC LSCC Discovery Study - Phosphoproteome",
            "proteome": "CPTAC LSCC Discovery Study - Proteome",
            "ubiquitylome": "CPTAC LSCC Discovery Study - Ubiquitylome",
        },
        "luad": {
            "acetylome": "CPTAC LUAD Discovery Study - Acetylome",
            "phosphoproteome": "S046-2",
            "proteome": "S046-1",
        },
        "ovarian": {
            "phosphoproteome_pnnl": "S038-3",
            "proteome_jhu": "S038-1",
            "proteome_pnnl": "S038-2",
        },
        "ucec": {
            "acetylome": "CPTAC UCEC Discovery Study - Acetylome",
            "phosphoproteome": "S043-2",
            "proteome": "S043-1",
        },
    }

    cancer_type = str.lower(cancer_type)
    if cancer_type != "all" and cancer_type not in study_submitter_id_map.keys():
        raise ex.InvalidParameterError(f"cancer_type must be either 'all' or one of the following:\n{list(study_submitter_id_map.keys())}\nYou passed '{cancer_type}'.")

    if cancer_type == "all":

        tables = {}
        for cancer, cancer_type_ids in study_submitter_id_map.items():
            tables[cancer_type] = _download_cancer_type(cancer_type_ids)

        return tables

    else:
        cancer_type_ids = study_submitter_id_map[cancer_type]
        return _download_cancer_type(cancer_type_ids)


def _download_cancer_type(cancer_type_ids):
    """Download the tables for all study_submitter_id for a particular cancer type."""

    tables = {}
    clin = {}

    for data_type, study_submitter_id in cancer_type_ids.items():

        clin_df = _pdc_download_study_clin(study_submitter_id)

        quant_df = _pdc_download_study_quant(study_submitter_id).\
        join(
            other=clin_df["case_submitter_id"],
            how="outer",
        ).\
        reset_index(drop=False).\
        set_index("case_submitter_id")

        clin_df = clin_df.\
        reset_index(drop=False).\
        set_index("case_submitter_id")

        clin[data_type] = clin_df
        tables[data_type] = quant_df

    master_clin = None

    for data_type in clin.keys():

        if master_clin is None:
            master_clin = clin[data_type].sort_index()
            print(f"\nMaster is {data_type}")
        else:
            new_clin = clin[data_type].sort_index()
            if not master_clin.equals(new_clin):
                print(f"Clin table was different: {data_type}")
            else:
                print(f"Clin table was same: {data_type}")

    tables["clin"] = master_clin

    return tables


def _pdc_download_study_clin(study_submitter_id):
    """Download PDC clinical data for a particular study."""

    clinical_query = '''
    query {
        clinicalPerStudy(study_submitter_id: "''' + study_submitter_id + '''") {
            case_id
            case_submitter_id
            age_at_diagnosis
            status
            ethnicity
            gender
            race
            morphology
            primary_diagnosis
            site_of_resection_or_biopsy
            tissue_or_organ_of_origin
            tumor_grade
            tumor_stage  
        }
    }
    '''

    result_json = _query_pdc(clinical_query)
    result_df = pd.\
    DataFrame(result_json["data"]["clinicalPerStudy"]).\
    set_index("case_id")

    return result_df

def _pdc_download_study_quant(study_submitter_id):
    """Download PDC quantitative data for a particular study."""

    proteome_query = '''
    query {
        quantDataMatrix(study_submitter_id: "''' + study_submitter_id + '''", data_type: "log2_ratio")
    }
    '''

    result_json = _query_pdc(proteome_query)
    result_df = pd.DataFrame(result_json["data"]["quantDataMatrix"])

    if result_df.shape[1] != 0:
        result_df = result_df.set_index(0).transpose()
        result_df = result_df.set_index(result_df.columns[0])
        result_df.index.name = "case_id"
    else:
        print(f"MISSING: {study_submitter_id} quant")

    return result_df

def _query_pdc(query):
    """Send a GraphQL query to the PDC and return the results."""

    url = 'https://pdc.cancer.gov/graphql'

    try:
        response = requests.post(url, json={'query': query})
        response.raise_for_status() # Raises a requests.HTTPError if the response code was unsuccessful

    except requests.RequestException: # Parent class for all exceptions in the requests module
        raise ex.NoInternetError("Insufficient internet. Check your internet connection.") from None

    return response.json()
