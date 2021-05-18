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

import cptac
from cptac.file_download import get_box_token
from cptac.exceptions import DatasetAlreadyInstalledWarning, InvalidParameterError, NoInternetError, PdcDownloadError

from .pancanbrca import SOURCES as BRCA_SOURCES
from .pancanccrcc import SOURCES as CCRCC_SOURCES
from .pancancoad import SOURCES as COAD_SOURCES
from .pancangbm import SOURCES as GBM_SOURCES
from .pancanhnscc import SOURCES as HNSCC_SOURCES
from .pancanlscc import SOURCES as LSCC_SOURCES
from .pancanluad import SOURCES as LUAD_SOURCES
from .pancanov import SOURCES as OV_SOURCES
from .pancanucec import SOURCES as UCEC_SOURCES

STUDY_IDS_MAP = {
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


def download(dataset, version="latest", redownload=False):

    dataset = dataset.lower()

    if dataset.startswith("pdc"):
        return _pdc_download(dataset, version=version, redownload=redownload)

    elif dataset.startswith("pancan") or dataset == "all":
        box_token = get_box_token()

        if dataset == "pancanbrca":
            sources = BRCA_SOURCES
        elif dataset == "pancanccrcc":
            sources = CCRCC_SOURCES
        elif dataset == "pancancoad":
            sources = COAD_SOURCES
        elif dataset == "pancangbm":
            sources = GBM_SOURCES
        elif dataset == "pancanhnscc":
            sources = HNSCC_SOURCES
        elif dataset == "pancanlscc":
            sources = LSCC_SOURCES
        elif dataset == "pancanluad":
            sources = LUAD_SOURCES
        elif dataset == "pancanov":
            sources = OV_SOURCES
        elif dataset == "pancanucec":
            sources = UCEC_SOURCES
        elif dataset == "all":
            sources = sorted(set(BRCA_SOURCES + CCRCC_SOURCES + COAD_SOURCES + GBM_SOURCES + HNSCC_SOURCES + LSCC_SOURCES + LUAD_SOURCES + OV_SOURCES + UCEC_SOURCES))
        else:
            raise InvalidParameterError(f"{dataset} is not a valid dataset.")

        overall_success = True
        for source in sources:

            if source.startswith("pdc"):
                single_success = download(source, version=version, redownload=redownload)
            else:
                single_success = cptac.download(source, version=version, redownload=redownload, _box_auth=True, _box_token=box_token)

            if not single_success:
                overall_success = False

        return overall_success

    else:
        return cptac.download(dataset, version=version, redownload=redownload, _box_auth=True)

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

# Helper functions
    
def _pdc_download(dataset, version, redownload):
    """Download data for the specified cancer type from the PDC."""

    dataset = str.lower(dataset)

    if dataset == "pdcall":
        overall_result = True
        for dataset in STUDY_IDS_MAP.keys():
            if not pdc_download(dataset, version, redownload):
                overall_result = False

        return overall_result

    if not dataset.startswith("pdc"):
        raise InvalidParameterError(f"pdc_download function can only be used for PDC datasets, which start with the prefix 'pdc'. You tried to download '{dataset}'.")

    if dataset not in STUDY_IDS_MAP.keys():
        raise InvalidParameterError(f"PDC dataset must be one of the following:\n{list(STUDY_IDS_MAP.keys())}\nYou passed '{dataset}'.")

    dataset_ids = STUDY_IDS_MAP[dataset]

    # Get the directory to where to store the data, and see if it exists
    path_here = os.path.abspath(os.path.dirname(__file__))
    cancer_dir = os.path.join(path_here, f"data_{dataset}")

    if os.path.isdir(cancer_dir):

        index_path = os.path.join(cancer_dir, "index.txt")

        # Check that they also have the index
        if not os.path.isfile(index_path):
            redownload = True
        else:
            # The PDC doesn't have a versioning scheme for the tables they serve, so originally we just called it version 0.0 but later decided it would be better to call it 1.0. So, check if theirs is called 0.0; if so, replace it with 1.0.

            with open(index_path, "r") as index_file:
                first_line = index_file.readline()

            if first_line.startswith("#0.0"):
                redownload=True

        if redownload:
            shutil.rmtree(cancer_dir)
        else:

            return True

    os.mkdir(cancer_dir)
    data_dir = os.path.join(cancer_dir, f"{dataset}_v1.0")
    os.mkdir(data_dir)

    # We'll combine all the clinical tables in case there are differences
    master_clin = pd.DataFrame()

    for data_type in dataset_ids.keys():

        # Print an update
        download_msg = f"Downloading {dataset} {data_type} files..."
        print(download_msg, end="\r")

        # Get the clinical and quantitative tables for the study ID
        clin, quant = download_pdc_id(dataset_ids[data_type], _download_msg=False)

        # Print a new update
        print(" " * len(download_msg), end="\r")
        save_msg = f"Saving {dataset} {data_type} files..."
        print(save_msg, end="\r")

        # Append the clinical dataframe
        master_clin = master_clin.append(clin)

        # Save the quantitative table
        quant.to_csv(os.path.join(data_dir, f"{data_type}.tsv.gz"), sep="\t")

        # Erase update
        print(" " * len(save_msg), end="\r")

    # Print an update
    save_msg = f"Saving {dataset} clinical file..."
    print(save_msg, end="\r")

    # Drop any duplicated rows in combined clinical table, then save it too
    master_clin = master_clin.drop_duplicates(keep="first")

    master_clin.to_csv(os.path.join(data_dir, "clinical.tsv.gz"), sep="\t")

    # Write a dummy index with just version numbers
    index_path = os.path.join(cancer_dir, "index.txt")

    with open(index_path, "w") as index_file:
        index_file.write("#1.0\n")
        
    # Erase update
    print(" " * len(save_msg), end="\r")

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
