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

# This file contains the scripts necessary for downloading data from the pdc source
import os
import pandas as pd
import requests
import threading
from cptac import CPTAC_BASE_DIR

from cptac.tools.download_tools.box_download import validate_version, get_index
from cptac.exceptions import CptacDevError, PdcDownloadError, NoInternetError, PdcDownloadError

STUDY_IDS_MAP = {
    "pdcbrca": {
        "acetylproteomics": "PDC000239", # Prospective Breast BI Acetylproteomics
        "phosphoproteomics": "PDC000121", # Prospective BRCA Phosphoproteomics S039-2
        "proteomics": "PDC000120", # Prospective BRCA Proteomics S039-1
    },
    "pdcccrcc": {
        "phosphoproteomics": "PDC000128", # CPTAC CCRCC Discovery Study - Phosphoproteme S044-2
        "proteomics": "PDC000127", # CPTAC CCRCC Discovery Study - Proteomics S044-1
    },
    "pdccoad": {
        "phosphoproteomics": "PDC000117", # Prospective COAD Phosphoproteomics S037-3
        "proteomics": "PDC000116", # Prospective COAD Proteomics S037-2
    },
    "pdcgbm": {
        "acetylproteomics": "PDC000245", # CPTAC GBM Discovery Study - Acetylproteomics
        "phosphoproteomics": "PDC000205", # CPTAC GBM Discovery Study - Phosphoproteomics
        "proteomics": "PDC000204", # CPTAC GBM Discovery Study - Proteomics
    },
    "pdchnscc": {
        "phosphoproteomics": "PDC000222", # CPTAC HNSCC Discovery Study - Phosphoproteomics
        "proteomics": "PDC000221", # CPTAC HNSCC Discovery Study - Proteomics
    },
    "pdclscc": {
        "acetylproteomics": "PDC000233", # CPTAC LSCC Discovery Study - Acetylproteomics
        "phosphoproteomics": "PDC000232", # CPTAC LSCC Discovery Study - Phosphoproteomics
        "proteomics": "PDC000234", # CPTAC LSCC Discovery Study - Proteomics
        "ubiquitylome": "PDC000237", # CPTAC LSCC Discovery Study - Ubiquitylome
    },
    "pdcluad": {
        "acetylproteomics": "PDC000224", # CPTAC LUAD Discovery Study - Acetylproteomics
        "phosphoproteomics": "PDC000149", # CPTAC LUAD Discovery Study - Phosphoproteomics
        "proteomics": "PDC000153", # CPTAC LUAD Discovery Study - Proteomics
    },
    "pdcov": {
        "phosphoproteomics": "PDC000119", # Prospective OV Phosphoproteomics S038-3
        "proteomics": "PDC000118", # Prospective OV Proteomics S038-2
    },
    "pdcpdac": {
        "proteomics": "PDC000270", # CPTAC PDAC Discovery Study - Proteomics
        "phosphoproteomics": "PDC000271", # CPTAC PDAC Discovery Study - Phosphoproteomics
    },
    "pdcucec": {
        "acetylproteomics": "PDC000226", # CPTAC UCEC Discovery Study - Acetylproteomics
        "phosphoproteomics": "PDC000126", # UCEC Discovery - Phosphoproteomics S043-2
        "proteomics": "PDC000125", # UCEC Discovery - Proteomics S043-1
    },
}

def pdc_download(cancer, datatypes, version, redownload):
    """Download data for the specified cancer type and datatype from the PDC."""

    studyID = "pdc" + cancer
    dataset_ids = STUDY_IDS_MAP[studyID]

    # filter out the datasets not requested
    ids_to_remove = set(dataset_ids.keys()) - set(datatypes)
    [ dataset_ids.pop(key) for key in ids_to_remove ]

    data_dir = os.path.join(CPTAC_BASE_DIR, f"data/data_pdc_{cancer}")

    # Validate the version number, including parsing if it's "latest"
    dataset = "pdc_" + cancer
    version_number = validate_version(version, dataset, use_context="download")

    # Construct the path to the directory for this version
    version_path = os.path.join(data_dir, f"{dataset}_v{version_number}")

    # Check that the index file exists. If not, there was an uncaught error in the mapping file download.
    index_path = os.path.join(data_dir, "index.txt")
    if not os.path.isfile(index_path):
        raise CptacDevError(f"Index file not found at {index_path}. Mapping file download probably failed.")

    # If any of the files are missing, we're going to delete any remaining and redownload all, in case the missing files are a sign of a previous data problem
    data_files = [f"{data_type}.tsv.gz" for data_type in dataset_ids.keys()] + ["clinical.tsv.gz"]
    for data_file in data_files:
        data_file_path = os.path.join(version_path, data_file)
        if not os.path.isfile(data_file_path):
            redownload = True
            break

    if redownload:
        for data_file in data_files:
            data_file_path = os.path.join(version_path, data_file)
            if os.path.isfile(data_file_path):
                os.remove(data_file_path)
    else:
        return True # If all the files are there and the user didn't ask to redownload, we're done.

    # We'll combine all the clinical tables in case there are differences
    master_clin = pd.DataFrame()

    for data_type in dataset_ids.keys():

        # Print an update
        download_msg = f"Downloading {cancer} {data_type} files..."
        print(download_msg, end="\r")

        # Get the clinical and quantitative tables for the study ID
        clin, quant = download_pdc_id(dataset_ids[data_type], _download_msg=False)

        # Print a new update
        print(" " * len(download_msg), end="\r")
        save_msg = f"Saving {cancer} {data_type} files..."
        print(save_msg, end="\r")

        # Append the clinical dataframe
        master_clin = pd.concat([master_clin, clin])

        # Save the quantitative table
        quant.to_csv(os.path.join(version_path, f"{data_type}.tsv.gz"), sep="\t")

        # Erase update
        print(" " * len(save_msg), end="\r")

    # Print an update
    save_msg = f"Saving {cancer} clinical file..."
    print(save_msg, end="\r")

    # Drop any duplicated rows in combined clinical table, then save it too
    master_clin = master_clin.drop_duplicates(keep="first")

    master_clin.to_csv(os.path.join(version_path, "clinical.tsv.gz"), sep="\t")

    # Erase update
    print(" " * len(save_msg), end="\r")

    return True

def _download_mapping_files():
    pass

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

    proteomics_query = '''
    query {
        quantDataMatrix(pdc_study_id: "''' + pdc_study_id + '''", data_type: "log2_ratio", acceptDUA: true)
    }
    '''

    result_json = _query_pdc(proteomics_query)
    result_df = pd.DataFrame(result_json["data"]["quantDataMatrix"])

    if result_df.shape[1] != 0:
        result_df = result_df.set_index(result_df.columns[0]).transpose()
    else:
        raise PdcDownloadError(f"quantDataMatrix table returned for PDC study ID {pdc_study_id} was empty.")

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
    """Check that the ids in the download function's STUDY_IDS_MAP match up."""
    
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

def download_pdc_id(pdc_id, _download_msg=True):
    """Download a PDC dataset by its PDC study id.
    
    Returns:
    pandas.DataFrame: The clinical table for the study id.
    pandas.DataFrame: The quantitative table for the study id.
    """

    if _download_msg:
        clin_msg = f"Downloading clinical table for {pdc_id}..."
        print(clin_msg, end="\r")

    # Download the clinical table
    clin = _download_study_clin(pdc_id).\
    set_index("case_submitter_id").\
    sort_index()

    if _download_msg:
        print(" " * len(clin_msg), end="\r")
        bio_msg = f"Downloading biospecimenPerStudy table for {pdc_id}..."
        print(bio_msg, end="\r")

    # The the biospecimenPerStudy table, which has both patient IDs and aliquot IDs
    bio = _download_study_biospecimen(pdc_id).\
    set_index("aliquot_submitter_id").\
    sort_index()

    if _download_msg:
        print(" " * len(bio_msg), end="\r")
        quant_msg = f"Downloading quantitative table for {pdc_id}..."
        print(quant_msg, end="\r")

    # Get the quantitative data table
    quant = _download_study_quant(pdc_id)

    if _download_msg:
        print(" " * len(quant_msg), end="\r")
        format_msg = f"Formatting tables for {pdc_id}..."
        print(format_msg, end="\r")

    # Join the patient IDs from the biospecimenPerStudy table into the quant table
    quant = quant.\
    assign(aliquot_submitter_id=quant.iloc[:, 0].str.split(":", n=1, expand=True)[1]).\
    drop(columns=quant.columns[0]).\
    set_index("aliquot_submitter_id").\
    sort_index()

    quant = bio.\
    join(quant, how="inner").\
    reset_index().\
    set_index(["case_submitter_id", "aliquot_submitter_id"]).\
    sort_index()

    # Clear message
    if _download_msg:
        print(" " * len(format_msg), end="\r")

    return clin, quant

def list_pdc_datasets():
    for dataset in STUDY_IDS_MAP.keys():
        print(f"Pdc{dataset[3:].title()}:")
        for data_type in STUDY_IDS_MAP[dataset].keys():
            print(f"\t{data_type}: {STUDY_IDS_MAP[dataset][data_type]}")
