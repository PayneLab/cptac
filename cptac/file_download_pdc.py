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
import shutil
import warnings

from .exceptions import DatasetAlreadyInstalledWarning, InvalidParameterError, NoInternetError

def pdc_download(dataset, version, redownload):
    """Download data for the specified cancer type from the PDC."""

    study_ids_map = {
        "pdcbrca": {
            "acetylome": "PDC000239", # Prospective Breast BI Acetylome
            "phosphoproteome": "PDC000121", # Prospective BRCA Phosphoproteome S039-2
            "proteome": "PDC000120", # Prospective BRCA Proteome S039-1
        },
        "pdcccrcc": {
            "phosphoproteome": "PDC000128", # CPTAC CCRCC Discovery Study - Phosphoproteme S044-2
            "proteome": "PDC000127", # CPTAC CCRCC Discovery Study - Proteome S044-1
        },
        "pdccoad": {
            "phosphoproteome": "PDC000117", # Prospective COAD Phosphoproteome S037-3
            "proteome": "PDC000116", # Prospective COAD Proteome S037-2
        },
        "pdcgbm": {
            "acetylome": "PDC000245", # CPTAC GBM Discovery Study - Acetylome
            "phosphoproteome": "PDC000205", # CPTAC GBM Discovery Study - Phosphoproteome
            "proteome": "PDC000204", # CPTAC GBM Discovery Study - Proteome
        },
        "pdchnscc": {
            "phosphoproteome": "PDC000222", # CPTAC HNSCC Discovery Study - Phosphoproteome
            "proteome": "PDC000221", # CPTAC HNSCC Discovery Study - Proteome
        },
        "pdclscc": {
            "acetylome": "PDC000233", # CPTAC LSCC Discovery Study - Acetylome
            "phosphoproteome": "PDC000232", # CPTAC LSCC Discovery Study - Phosphoproteome
            "proteome": "PDC000234", # CPTAC LSCC Discovery Study - Proteome
            "ubiquitylome": "PDC000237", # CPTAC LSCC Discovery Study - Ubiquitylome
        },
        "pdcluad": {
            "acetylome": "PDC000224", # CPTAC LUAD Discovery Study - Acetylome
            "phosphoproteome": "PDC000149", # CPTAC LUAD Discovery Study - Phosphoproteome
            "proteome": "PDC000153", # CPTAC LUAD Discovery Study - Proteome
        },
        "pdcov": {
            "phosphoproteome": "PDC000119", # Prospective OV Phosphoproteome S038-3
            "proteome": "PDC000118", # Prospective OV Proteome S038-2
        },
        "pdcpda": {
            "proteome": "PDC000270", # CPTAC PDA Discovery Study - Proteome
            "phosphoproteome": "PDC000271", # CPTAC PDA Discovery Study - Phosphoproteome
        },
        "pdcucec": {
            "acetylome": "PDC000226", # CPTAC UCEC Discovery Study - Acetylome
            "phosphoproteome": "PDC000126", # UCEC Discovery - Phosphoproteome S043-2
            "proteome": "PDC000125", # UCEC Discovery - Proteome S043-1
        },
    }

    dataset = str.lower(dataset)

    if dataset == "pdcall":
        overall_result = True
        for dataset in study_ids_map.keys():
            if not pdc_download(dataset, version, redownload):
                overall_result = False

        return overall_result

    if not dataset.startswith("pdc"):
        raise InvalidParameterError(f"pdc_download function can only be used for PDC datasets, which start with the prefix 'pdc'. You tried to download '{dataset}'.")

    if dataset not in study_ids_map.keys():
        raise InvalidParameterError(f"PDC dataset must be one of the following:\n{list(study_ids_map.keys())}\nYou passed '{dataset}'.")

    dataset_ids = study_ids_map[dataset]

    # Get the directory to where to store the data, and see if it exists
    path_here = os.path.abspath(os.path.dirname(__file__))
    cancer_dir = os.path.join(path_here, "pancan", f"data_{dataset}")

    if os.path.isdir(cancer_dir):

        # Only redownload if they explicitly wanted that
        if redownload:
            shutil.rmtree(cancer_dir)
        else:
            warnings.warn(f"The {dataset} dataset has been downloaded previously. If you wish to erase the existing download and re-download it, pass 'redownload=True'.", DatasetAlreadyInstalledWarning, stacklevel=4)

            return True

    os.mkdir(cancer_dir)
    data_dir = os.path.join(cancer_dir, f"{dataset}_v0.0")
    os.mkdir(data_dir)

    # We'll combine all the clinical tables in case there are differences
    master_clin = pd.DataFrame()

    for data_type in dataset_ids.keys():

        print(data_type)

        # Append the clinical dataframe
        clin = _download_study_clin(dataset_ids[data_type])
        master_clin = master_clin.append(clin)

        # Get the quantDataMatrix table, then join in patient IDs from the biospecimenPerStudy table
        bio = _download_study_biospecimen(dataset_ids[data_type]).\
        set_index("aliquot_submitter_id").\
        sort_index()

        quant = _download_study_quant(dataset_ids[data_type])

        quant = quant.\
        assign(aliquot_submitter_id=quant.iloc[:, 0].str.split(":", n=1, expand=True)[1]).\
        drop(columns=quant.columns[0]).\
        set_index("aliquot_submitter_id").\
        sort_index()

        quant = bio.join(quant, how="right")

        # Save the quantitative table
        quant.to_csv(os.path.join(data_dir, f"{data_type}.tsv.gz"), sep="\t")

    # Drop any duplicated rows in combined clinical table, then save it too
    master_clin = master_clin.\
    set_index("case_submitter_id").\
    drop_duplicates(keep="first")

    master_clin.to_csv(os.path.join(data_dir, "clinical.tsv.gz"), sep="\t")

    # Write a dummy index with just version numbers
    index_path = os.path.join(cancer_dir, "index.txt")

    with open(index_path, "w") as index_file:
        index_file.write("#0.0\n")
        
    return True

def _download_study_clin(pdc_study_id):
    """Download PDC clinical data for a particular study."""

    clinical_query = '''
    query {
        clinicalPerStudy(pdc_study_id: "''' + pdc_study_id + '''", acceptDUA: true) {
            age_at_diagnosis, ajcc_clinical_m, ajcc_clinical_n, ajcc_clinical_stage, ajcc_clinical_t, ajcc_pathologic_m,
            ajcc_pathologic_n, ajcc_pathologic_stage, ajcc_pathologic_t, ann_arbor_b_symptoms, ann_arbor_clinical_stage,
            ann_arbor_extranodal_involvement, ann_arbor_pathologic_stage, best_overall_response, burkitt_lymphoma_clinical_variant,
            case_id, case_submitter_id, cause_of_death, circumferential_resection_margin, classification_of_tumor, colon_polyps_history,
            days_to_best_overall_response, days_to_birth, days_to_death, days_to_diagnosis, days_to_hiv_diagnosis, days_to_last_follow_up,
            days_to_last_known_disease_status, days_to_new_event, days_to_recurrence, demographic_id, demographic_submitter_id,
            diagnosis_id, diagnosis_submitter_id, disease_type, ethnicity, figo_stage, gender, hiv_positive, hpv_positive_type, hpv_status,
            icd_10_code, iss_stage, last_known_disease_status, laterality, ldh_level_at_diagnosis, ldh_normal_range_upper,
            lymphatic_invasion_present, lymph_nodes_positive, method_of_diagnosis, morphology, new_event_anatomic_site, new_event_type,
            overall_survival, perineural_invasion_present, primary_diagnosis, primary_site, prior_malignancy, prior_treatment,
            progression_free_survival, progression_free_survival_event, progression_or_recurrence, race, residual_disease,
            site_of_resection_or_biopsy, status, synchronous_malignancy, tissue_or_organ_of_origin, tumor_cell_content, tumor_grade,
            tumor_stage, vascular_invasion_present, vital_status, year_of_birth, year_of_death, year_of_diagnosis
        }
    }
    '''

    result_json = _query_pdc(clinical_query)
    result_df = pd.\
    DataFrame(result_json["data"]["clinicalPerStudy"])

    return result_df

def _download_study_biospecimen(pdc_study_id):
    """Download PDC biospecimen data for a particular study."""

    biospecimen_query = '''
    query {
        biospecimenPerStudy(pdc_study_id: "''' + pdc_study_id + '''", acceptDUA: true) {
            aliquot_submitter_id
            case_submitter_id
        }
    }
    '''

    result_json = _query_pdc(biospecimen_query)
    result_df = pd.\
    DataFrame(result_json["data"]["biospecimenPerStudy"])

    return result_df

def _download_study_quant(pdc_study_id):
    """Download PDC quantitative data for a particular study."""

    proteome_query = '''
    query {
        quantDataMatrix(pdc_study_id: "''' + pdc_study_id + '''", data_type: "log2_ratio", acceptDUA: true)
    }
    '''

    result_json = _query_pdc(proteome_query)
    result_df = pd.DataFrame(result_json["data"]["quantDataMatrix"])

    if result_df.shape[1] != 0:
        result_df = result_df.set_index(result_df.columns[0]).transpose()
    else:
        print(f"MISSING: {pdc_study_id} quant")

    return result_df

def _query_pdc(query):
    """Send a GraphQL query to the PDC and return the results."""

    url = 'https://pdc.cancer.gov/graphql'

    try:
        response = requests.post(url, json={'query': query})
        response.raise_for_status() # Raises a requests.HTTPError if the response code was unsuccessful

    except requests.RequestException: # Parent class for all exceptions in the requests module
        raise NoInternetError("Insufficient internet. Check your internet connection.") from None

    return response.json()

def _check_ids_match(ids_map):
    """Check that the ids in the download function's study_ids_map match up."""
    
    for cancer in ids_map.values():
        for data in cancer.values():
            pdc_study_id = data["pdc_study_id"]
            study_submitter_id = data["study_submitter_id"]

            query = '''
            query {
              study (pdc_study_id: "''' + pdc_study_id + '''" acceptDUA: true) {
                pdc_study_id,
                study_submitter_id
              }
            }
            '''

            idres = _query_pdc(query)

            server_psi = idres["data"]["study"][0]["pdc_study_id"]
            server_ssi = idres["data"]["study"][0]["study_submitter_id"]

            assert server_psi == pdc_study_id
            assert server_ssi == study_submitter_id

            print(f"{server_psi} == {pdc_study_id}")
            print(f"{server_ssi} == {study_submitter_id}")
            print()
