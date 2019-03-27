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

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import CPTAC.Endometrial as en
from utilities import *

def check_dataframe(name, df, exp_dim, exp_headers, coordinates, values): # private
    """Test a dataframe's dimensions and headers, and three test values, then print whether it passed the test.

    Parameters
    name: string containing the name of the dataframe gotten by the getter we're testing
    df: the dataframe gotten by the getter we are testing
    exp_dim: a tuple containing the expected dimensions of the dataframe, in the format (rows, columns)
    exp_headers: if the dataframe has up to 20 columns, all of the headers for the dataframe, in order. If it has more than 20 columns, then a list containing the first ten and last ten headers, in order.
    coordinates: a tuple with three elements, each element being a tuple with two elements, the first element being the int index of the row of a test value, and the second element being the int index of the column of a test value
    values: a tuple with three elements, each element being the expected value of the test value corresponding to the coordinates at the same index in the coordinates parameter 
    """
    PASS = True

    # Check that we actually got a dataframe
    if df is None:
        print("Error loading {} dataframe. Getter returned None.".format(name))

    # Check dimensions
    act_dim = df.shape
    if exp_dim != act_dim:
        print("Error: {} dataframe dimensions did not match expected values.\n\tExpected: {}\n\tActual: {}\n".format(name, exp_dim, act_dim))
        PASS = False

    # Check headers
    act_headers_all = list(df.columns.values)
    if len(df.columns.values) <= 20:
        act_headers = act_headers_all
    else:
        act_headers = act_headers_all[:10] + act_headers_all[-10:]

    if len(exp_headers) != len(act_headers):
        print("Unexpected number of test headers in {} dataframe. Expected number of headers: {}. You passed {} headers.\n".format(name, len(act_headers), len(exp_headers)))
        PASS = False
    else:
        for i, header in enumerate(exp_headers):
            if header != act_headers[i]:
                print("Error: {} dataframe header did not match expected value.\n\tExpected: {}\n\tActual: {}\n".format(name, header, act_headers[i]))
                PASS = False

    # Check test values
    act_values = [
        df.iloc[coordinates[0][0], coordinates[0][1]],
        df.iloc[coordinates[1][0], coordinates[1][1]],
        df.iloc[coordinates[2][0], coordinates[2][1]]
    ]

    for i, value in enumerate(values):
        if act_values[i] != value:
            print("Error: {} dataframe value did not match expected value.\n\tColumn: {}\n\tIndex: {}\n\tExpected: {}\n\tActual: {}\n".format(name, df.columns.values[coordinates[i][1]], df.index.values[coordinates[i][0]], value, act_values[i]))
            PASS = False

    # Print whether the dataframe passed the test
    if PASS:
        print("PASS")
    else:
        print("FAIL\n")

def test_get_clinical_filtered():
    """Test get_clinical() with the default parameter unfiltered=False."""

    print('Testing get_clinical with the default parameter unfiltered=False...')

    clinical_name = "Clinical"
    clinical_df = en.get_clinical()
    clinical_dim = (144, 27)
    clinical_headers = ['Proteomics_Participant_ID', 'Case_excluded', 'Proteomics_Tumor_Normal', 'Country', 'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity', 'Path_Stage_Primary_Tumor-pT', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site', 'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm', 'Num_full_term_pregnancies']
    clinical_test_coord = ((79, 16), (15, 25), (88, 2))
    clinical_test_vals = (23.88, 3.2, 'Tumor')

    check_dataframe(clinical_name, clinical_df, clinical_dim, clinical_headers, clinical_test_coord, clinical_test_vals)

def test_get_clinical_unfiltered():
    """Test get_clinical with parameter unfiltered=True."""

    print('Testing get_clinical with parameter unfiltered=True...')

    clinical_excluded_name = "Clinical (with unfiltered samples)"
    clinical_excluded_df = en.get_clinical(unfiltered=True)
    clinical_excluded_dim = (153, 27)
    clinical_excluded_headers = ['Proteomics_Participant_ID', 'Case_excluded', 'Proteomics_Tumor_Normal', 'Country', 'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity', 'Path_Stage_Primary_Tumor-pT', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site', 'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm', 'Num_full_term_pregnancies']
    clinical_excluded_test_coord = ((23, 8), (151, 1), (32, 26))
    clinical_excluded_test_vals = ('Normal', 'No', '3')

    check_dataframe(clinical_excluded_name, clinical_excluded_df, clinical_excluded_dim, clinical_excluded_headers, clinical_excluded_test_coord, clinical_excluded_test_vals)
    print("The unfiltered data warning above was expected.") # To avoid confusion

def test_get_derived_molecular_filtered():
    """Test get_derived_molecular() with default parameter unfiltered=False."""

    print('Testing get_derived_molecular() with default parameter unfiltered=False...')
    print("UNDER CONSTRUCTION")

def test_get_derived_molecular_unfiltered():
    """Test get_derived_molecular with parameter unfiltered=True."""

    print('Testing get_derived_molecular() with parameter unfiltered=True...')
    print("UNDER CONSTRUCTION")

def test_get_acetylproteomics_filtered():
    """Test get_acetylproteomics() with default parameter unfiltered=False."""

    print('Test get_acetylproteomics() with default parameter unfiltered=False...')
    print("UNDER CONSTRUCTION")

def test_get_acetylproteomics_unfiltered():
    """Test get_acetylproteomics with parameter unfiltered=True."""

    print('Testing get_acetylproteomics with parameter unfiltered=True...')
    print("UNDER CONSTRUCTION")

def test_get_proteomics():
    """Test get_proteomics()."""

    print('Testing get_proteomics()...')

    proteomics_name = "Proteomics"
    proteomics_df = en.get_proteomics()
    proteomics_dim = (144, 10999)
    proteomics_headers = ['A1BG', 'A2M', 'A2ML1', 'A4GALT', 'AAAS', 'AACS', 'AADAT', 'AAED1', 'AAGAB', 'AAK1', 'ZSWIM8', 'ZSWIM9', 'ZW10', 'ZWILCH', 'ZWINT', 'ZXDC', 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3']
    proteomics_test_coord = ((34, 6003), (99, 9544), (143, 32))
    proteomics_test_vals = (0.0461, 1.68, 0.904)

    check_dataframe(proteomics_name, proteomics_df, proteomics_dim, proteomics_headers, proteomics_test_coord, proteomics_test_vals)

def test_get_transcriptomics_linear():
    """Test get_transcriptomics() with default parameter data_type="linear"."""

    print('Testing get_transcriptomics() with default parameter data_type="linear"...')

    transcriptomics_linear_name = "Transcriptomics (linear)"
    transcriptomics_linear_df = en.get_transcriptomics()
    transcriptomics_linear_dim = (109, 28057)
    transcriptomics_linear_headers = ['A1BG', 'A1BG-AS1', 'A1CF', 'A2M', 'A2M-AS1', 'A2ML1', 'A2MP1', 'A3GALT2', 'A4GALT', 'A4GNT', 'ZWILCH', 'ZWINT', 'ZXDA', 'ZXDB', 'ZXDC', 'ZYG11A', 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3']
    transcriptomics_linear_test_coord = ((22, 25483), (108, 23), (101, 17748))
    transcriptomics_linear_test_vals = (0.82, 12.0, 6.19)

    check_dataframe(transcriptomics_linear_name, transcriptomics_linear_df, transcriptomics_linear_dim, transcriptomics_linear_headers, transcriptomics_linear_test_coord, transcriptomics_linear_test_vals)

def test_get_transcriptomics_circular():
    """Test get_transcriptomics() with parameter data_type="circular"."""

    print('Testing get_transcriptomics() with parameter data_type="circular"...')

    transcriptomics_circular_name = "Transcriptomics (circular)"
    transcriptomics_circular_df = en.get_transcriptomics(data_type="circular")
    transcriptomics_circular_dim = (109, 4945)
    transcriptomics_circular_headers = ['circ_chr10_100260218_100262063_CWF19L1', 'circ_chr10_100923975_100926019_SLF2', 'circ_chr10_100923978_100926019_SLF2', 'circ_chr10_100937402_100944128_SLF2', 'circ_chr10_100937402_100950753_SLF2', 'circ_chr10_101584602_101586156_POLL', 'circ_chr10_101667886_101676436_FBXW4', 'circ_chr10_101672915_101676436_FBXW4', 'circ_chr10_101792839_101807901_OGA', 'circ_chr10_101792839_101810314_OGA', 'circ_chrX_80288906_80310233_CHMP1B2P', 'circ_chrX_80289664_80310233_CHMP1B2P', 'circ_chrX_80707427_80719656_BRWD3', 'circ_chrX_80791854_80793772_BRWD3', 'circ_chrX_84096194_84164387_RPS6KA6', 'circ_chrX_84134782_84164387_RPS6KA6', 'circ_chrX_85067127_85074391_APOOL', 'circ_chrX_85978767_85981809_CHM', 'circ_chrX_91414904_91418871_PABPC5-AS1', 'circ_chrX_9691579_9693419_TBL1X']
    transcriptomics_circular_test_coord = ((108, 1), (30, 4935), (73, 2003))
    transcriptomics_circular_test_vals = (9.08, 6.56, 0.0)

    check_dataframe(transcriptomics_circular_name, transcriptomics_circular_df, transcriptomics_circular_dim, transcriptomics_circular_headers, transcriptomics_circular_test_coord, transcriptomics_circular_test_vals)

def test_get_transcriptomics_miRNA():
    """Test get_transcriptomics() with parameter data_type="miRNA"."""

    print('Testing get_transcriptomics() with parameter data_type="miRNA"...')

    mirna_name = "Transcriptomics (miRNA)"
    mirna_df = en.get_transcriptomics(data_type="miRNA")
    mirna_dim = (99, 2337)
    mirna_headers = ['hsa-let-7a-2-3p', 'hsa-let-7a-3p', 'hsa-let-7a-5p', 'hsa-let-7b-3p', 'hsa-let-7b-5p', 'hsa-let-7c-3p', 'hsa-let-7c-5p', 'hsa-let-7d-3p', 'hsa-let-7d-5p', 'hsa-let-7e-3p', 'hsa-miR-9901', 'hsa-miR-9902', 'hsa-miR-9903', 'hsa-miR-9983-3p', 'hsa-miR-9985', 'hsa-miR-9986', 'hsa-miR-99a-3p', 'hsa-miR-99a-5p', 'hsa-miR-99b-3p', 'hsa-miR-99b-5p']
    mirna_test_coord = ((5, 0), (98, 1597), (54, 2231))
    mirna_test_vals = (1.79, 1.36, 0.26)
    
    check_dataframe(mirna_name, mirna_df, mirna_dim, mirna_headers, mirna_test_coord, mirna_test_vals)

def test_get_transcriptomics_with_invalid():
    """Test get_transcriptomics() with an invalid parameter, and make sure that it raises an exception."""

    print("Testing get_transcriptomics() with an invalid parameter, and make sure that it raises an exception...")

    try:
        en.get_transcriptomics("gobbledegook")
    except ValueError:
        print("PASS")
    else:
        print("Error: get_transcriptomics() did not raise ValueError as expected, when given invalid parameter.FAIL")

def test_get_CNA():
    """Test get_CNA()."""

    print('Testing get_CNA()...')

    cna_name = "CNA"
    cna_df = en.get_CNA()
    cna_dim = (95, 28057)
    cna_headers = ['MFSD14A', 'SASS6', 'TRMT13', 'LRRC39', 'DBT', 'RTCA-AS1', 'RTCA', 'MIR553', 'UBE4B', 'CDC14A', 'TSPY8', 'FAM197Y2', 'FAM197Y4', 'FAM197Y5', 'FAM197Y7', 'FAM197Y8', 'FAM197Y6', 'FAM197Y3', 'RBMY3AP', 'TTTY22']
    cna_test_coord = ((12, 27865), (60, 8), (94, 15439))
    cna_test_vals = (-0.07, 0.01, 0.03)

    check_dataframe(cna_name, cna_df, cna_dim, cna_headers, cna_test_coord, cna_test_vals)

def test_get_phosphoproteomics_site():
    """Test get_phosphoproteomics() with default parameter gene_level=False."""

    print('Testing get_phosphoproteomics() with default parameter gene_level=False...')

    phosphoproteomics_site_name = "Phosphoproteomics (site)"
    phosphoproteomics_site_df =  en.get_phosphoproteomics()
    phosphoproteomics_site_dim = (144, 73212)
    phosphoproteomics_site_headers = ['AAAS-S495', 'AAAS-S541', 'AAAS-Y485', 'AACS-S618', 'AAED1-S12', 'AAGAB-S310', 'AAGAB-S311', 'AAK1-S14', 'AAK1-S18', 'AAK1-S20', 'ZZZ3-S397', 'ZZZ3-S411', 'ZZZ3-S420', 'ZZZ3-S424', 'ZZZ3-S426', 'ZZZ3-S468', 'ZZZ3-S89', 'ZZZ3-T415', 'ZZZ3-T418', 'ZZZ3-Y399']
    phosphoproteomics_site_test_coord = ((36, 46), (12, 72436), (96, 45361))
    phosphoproteomics_site_test_vals = (0.579, 0.669, 0.156)

    check_dataframe(phosphoproteomics_site_name, phosphoproteomics_site_df, phosphoproteomics_site_dim, phosphoproteomics_site_headers, phosphoproteomics_site_test_coord, phosphoproteomics_site_test_vals)

def test_get_phosphoproteomics_gene():
    """Test get_phosphoproteomics() with parameter gene_level=True."""

    print('Testing get_phosphoproteomics() with parameter gene_level=True...')

    phosphoproteomics_gene_name = "Phosphoproteomics (gene)"
    phosphoproteomics_gene_df = en.get_phosphoproteomics(gene_level=True)
    phosphoproteomics_gene_dim = (144, 8466)
    phosphoproteomics_gene_headers = ['AAAS', 'AACS', 'AAED1', 'AAGAB', 'AAK1', 'AAMDC', 'AARS', 'AASDH', 'AATF', 'ABCA1', 'ZSCAN5C', 'ZSWIM3', 'ZSWIM8', 'ZUP1', 'ZW10', 'ZXDA', 'ZXDC', 'ZYX', 'ZZEF1', 'ZZZ3']
    phosphoproteomics_gene_test_coord =  ((2, 7999), (143, 1045), (71, 6543))
    phosphoproteomics_gene_test_vals = (-0.0879, 0.929, 0.153)

    check_dataframe(phosphoproteomics_gene_name, phosphoproteomics_gene_df, phosphoproteomics_gene_dim, phosphoproteomics_gene_headers, phosphoproteomics_gene_test_coord, phosphoproteomics_gene_test_vals)

def test_get_phosphosites():
    """Test get_phosphosites."""

    print('Testing get_phosphosites...')

    phosphosites_name = "Phosphosites for the AAK1-S14 gene"
    phosphosites_df = en.get_phosphosites('AAK1-S14')
    phosphosites_dim = (144, 1)
    phosphosites_headers = ['AAK1-S14']
    phosphosites_test_coord = ((27, 0), (76, 0), (128, 0))
    phosphosites_test_vals = (0.603, -0.272, 0.1395)

    check_dataframe(phosphosites_name, phosphosites_df, phosphosites_dim, phosphosites_headers, phosphosites_test_coord, phosphosites_test_vals)

def test_get_somatic_maf():
    """Test get_somatic() with default parameters binary=False, unparsed=False (this will return the Somatic Maf dataframe)."""

    print('Testing get_somatic() with default parameters binary=False, unparsed=False (this will return the Somatic Maf dataframe)...')

    somatic_name = "Somatic (maf)"
    somatic_df = en.get_somatic()
    somatic_dim = (52560, 5)
    somatic_headers = ['Clinical_Patient_Key', 'Patient_Id', 'Gene', 'Mutation', 'Location']
    somatic_test_coord = ((52000, 3), (12, 4), (34567, 0))
    somatic_test_vals = ('Missense_Mutation', 'p.T2121P', 'S059')

    check_dataframe(somatic_name, somatic_df, somatic_dim, somatic_headers, somatic_test_coord, somatic_test_vals)

def test_get_somatic_binary():
    """Test get_somatic with parameter binary=True and therefore unparsed=False."""

    print('Testing get_somatic with parameter binary=True and therefore unparsed=False...')

    somatic_binary_name = "Somatic (binary)"
    somatic_binary_df = en.get_somatic(binary=True)
    somatic_binary_dim = (95, 51559)
    somatic_binary_headers = ['A1BG_p.E298K', 'A1BG_p.S181N', 'A1CF_p.F487L', 'A1CF_p.S236Y', 'A2ML1_p.A8V', 'A2ML1_p.G1306D', 'A2ML1_p.L1347F', 'A2ML1_p.L82I', 'A2ML1_p.P712S', 'A2ML1_p.R443Q', 'ZYG11A_p.Q442H', 'ZYG11B_p.H315R', 'ZYG11B_p.R495M', 'ZYG11B_p.R728C', 'ZYX_p.C447Y', 'ZZEF1_p.A2723V', 'ZZEF1_p.D845Y', 'ZZEF1_p.K1251E', 'ZZEF1_p.K2387Sfs*40', 'ZZZ3_p.Y891C']
    somatic_binary_test_coord = ((94, 51558), (0, 0), (45, 25436))
    somatic_binary_test_vals = (0, 0, 0)

    check_dataframe(somatic_binary_name, somatic_binary_df, somatic_binary_dim, somatic_binary_headers, somatic_binary_test_coord, somatic_binary_test_vals)

def test_get_somatic_unparsed():
    """Test get_somatic with parameter unparsed=True and therefore binary=False."""

    print('Testing get_somatic with parameter unparsed=True and therefore binary=False...')

    somatic_unparsed_name = "Somatic (unparsed)"
    somatic_unparsed_df = en.get_somatic(unparsed=True)
    somatic_unparsed_dim = (53101, 124)
    somatic_unparsed_headers = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'ExAC_AC_AN_Adj', 'ExAC_AC_AN', 'ExAC_AC_AN_AFR', 'ExAC_AC_AN_AMR', 'ExAC_AC_AN_EAS', 'ExAC_AC_AN_FIN', 'ExAC_AC_AN_NFE', 'ExAC_AC_AN_OTH', 'ExAC_AC_AN_SAS', 'ExAC_FILTER']
    somatic_unparsed_test_coord = ((52265, 45), (12, 70), (27658, 1))
    somatic_unparsed_test_vals = ('strelkasnv-varssnv-mutectsnv', 'UPI0000167B91', 0)

    check_dataframe(somatic_unparsed_name, somatic_unparsed_df, somatic_unparsed_dim, somatic_unparsed_headers, somatic_unparsed_test_coord, somatic_unparsed_test_vals)

def evaluate_special_getters():
    print("Evaluating special getters...")
    results = []
    functions = {}
    results.append(en.get_clinical_cols()); functions[len(results)] = "clinical_cols"
    results.append(en.get_cohort_clinical(["Diabetes","BMI"])); functions[len(results)] = "cohort_meta"
    results.append(en.get_proteomics_quant(["S018","S100"])); functions[len(results)] = "proteomics_quant"
    results.append(en.get_proteomics_cols()); functions[len(results)] = "proteomics_cols"
    results.append(en.get_transcriptomics_cols()); functions[len(results)] = "transcriptomics_cols"
    results.append(en.get_cohort_proteomics(["A1BG","TP53"])); functions[len(results)] = "cohort_proteomics"
    results.append(en.get_cohort_transcriptomics(["A1BG","TP53"])); functions[len(results)] = "cohort_transcriptomics"
    results.append(en.get_cohort_cna(["SASS6","TTTY22"])); functions[len(results)] = "cohort_cna"
    results.append(en.get_cohort_phosphoproteomics(["TP53-S315","AAAS-S541"])); functions[len(results)] = "cohort_phosphoproteomics"
    results.append(en.get_patient_mutations("C3L-00157")); functions[len(results)] = "patient_mutations(Patient_Id)"
    results.append(en.get_patient_mutations("S013")); functions[len(results)] = "patient_mutations(Clinical_Patient_Key)"
    results.append(en.get_phosphosites("TP53")); functions[len(results)] = "phosphosites"
    PASS = True
    for x in range(0,len(results)):
        if results[x] is None:
            print("Error with get",functions[x+1], "function")
            PASS = False
    if PASS:
        print("PASS")
    else:
        print("FAIL")
def evaluate_utilities(): #compare_**** functions
    print("Evaluating utilities...")
    results = []
    functions = {}
    results.append(en.compare_gene(en.get_proteomics(), en.get_transcriptomics(), "A1BG")); functions[len(results)] = "compare_gene"
    results.append(en.compare_gene(en.get_proteomics(), en.get_transcriptomics(), ["A1BG","RPL11"])); functions[len(results)] = "compare_genes"
    results.append(en.compare_clinical(en.get_proteomics(), "BMI")); functions[len(results)] = "compare_clinical"
    results.append(en.compare_mutations(en.get_proteomics(),"TP53")); functions[len(results)] = "compare_mutations(Proteomics)"
    results.append(en.compare_mutations(en.get_proteomics(),"TP53","AURKA")); functions[len(results)] = "compare_mutations(Proteomics with Somatic)"
    results.append(en.compare_mutations(en.get_phosphoproteomics(), "IRS2")); functions[len(results)] = "compare_mutations(Phosphoproteomics)"
    results.append(en.compare_mutations(en.get_phosphoproteomics(), "IRS2","PIK3CA")); functions[len(results)] = "compare_mutations(Phosphoproteomics with Somatic)"
    results.append(en.compare_mutations_full(en.get_proteomics(),"TP53")); functions[len(results)] = "compare_mutations_full(Proteomics)"
    results.append(en.compare_mutations_full(en.get_proteomics(),"TP53","AURKA")); functions[len(results)] = "compare_mutations_full(Proteomics with Somatic)"
    results.append(en.compare_mutations_full(en.get_phosphoproteomics(), "IRS2")); functions[len(results)] = "compare_mutations_full(Phosphoproteomics)"
    results.append(en.compare_mutations_full(en.get_phosphoproteomics(), "IRS2","PIK3CA")); functions[len(results)] = "compare_mutations_full(Phosphoproteomics with Somatic)"
    results.append(en.compare_phosphosites("TP53")); functions[len(results)] = "compare_phosphosites"
    PASS = True
    for x in range(0,len(results)):
        if results[x] is None:
            print("Error with",functions[x+1],"function")
            PASS = False
    if PASS:
        print("PASS")
    else:
        print("FAIL")

def test_merged_column(original_df, merged_df, original_header, merged_header, merged_df_name): # private
    """
    Parameters
    original_df: the dataframe the column was taken from
    merged_df: the merged dataframe with the column
    original_header: the column's header in the original dataframe
    merged_header: the column's header in the merged dataframe
    merged_name: the name of the merged dataframe, in case we need to print an informative error message

    Returns
    bool indicating whether the column in the merged dataframe and the column in the original dataframe had the same values for each index
    """
    PASS = True

    for sample in merged_df.index.values:
        original_value = original_df.loc[sample, original_header]
        merged_value = merged_df.loc[sample, merged_header]
        if (merged_value != original_value) and (pd.notna(merged_value) or pd.notna(original_value)):
            print("Merged dataframe had incorrect values.\n\tDataframe: {}\n\tSample: {}\tColumn: {}\n\tExpected: {}\tActual: {}\n".format(merged_df_name, sample, merged_header, original_value, merged_value))
            PASS = False

    return PASS

def test_merged_column_from_row(source_df, merged_df, ID_column, filter_column, filter_value, source_column, dest_column, merged_df_name):
    """
    Parameters
    source_df: dataframe the data came from
    merged_df: dataframe the data was put in
    ID_column: string indiating the column in source_df that has the ID values that are the indices in merged_df
    filter_column: string indicating the column whose value was looked at to decide whether to take the value in the source column for a particular row, and put it in merged_df
    filter_value: the value in filter_column that indicates the data from source_column for that sample should go in merged_df
    source_column: string indicating the column in source_df from which data was taken
    dest_column: string indicating the column in merged_df where the data from source_df was put
    merged_df_name: string with the name of the merged dataframe, in case we need to print an informative error message.

    Returns
    bool indicating whether, for each sample and filter value, the data in merged_df matched the data in source_df for that sample ID and filter value.
    """
    PASS = True

    for sample in merged_df.index.values:
        sample_source_df = source_df.loc[source_df[ID_column] == sample] # Load a dataframe with all just the values from source_df for this sample
        source_filtered_df = sample_source_df.loc[sample_source_df[filter_column] == filter_value]
        original_values = source_filtered_df[source_column].values

        if len(original_values) == 0:
            if source_df.name.startswith('somatic'):
                if sample <= 'S104':
                    original_value = 'Wildtype_Tumor'
                else:
                    original_value = 'Wildtype_Normal'
            else:
                original_value == float('NaN')
        elif len(original_values) == 1:
            original_value = original_values[0]
        elif len(original_values) > 1 and source_df.name.startswith('somatic'):
            source_filtered_with_hierarchy = Utilities().add_mutation_hierarchy(source_filtered_df)
            source_filtered_with_hierarchy = source_filtered_with_hierarchy.sort_values(by = [ID_column, 'Mutation_Hierarchy'], ascending = [True,False]) #sorts by patient key, then by hierarchy so the duplicates will come with the lower number first
            original_value = source_filtered_with_hierarchy[source_column].iloc[0]
        else:
            raise ValueError('Unexpected duplicate entries in source dataframe for merged dataframe.\n\tSource dataframe: {}\n\tMerged dataframe: {}\n\tSample: {}\n\tColumn: {}\n\tValues found: {}\n'.format(source_df.name, merged_df_name, sample, source_column, original_values))

        merged_value = merged_df.loc[sample, dest_column]
        if (merged_value != original_value) and (pd.notna(merged_value) or pd.notna(original_value)):
            print("Merged dataframe had incorrect value.\n\tDataframe: {}\n\tSample: {}\tColumn: {}\n\tExpected: {}\tActual: {}\n".format(merged_df_name, sample, dest_column, original_value, merged_value))
            PASS = False

    return PASS

def evaluate_utilities_v2():
    # We will test all of the compare_**** functions, which either merge dataframes, or add columns to a dataframe
    # When dataframes are merged, we will make sure that the data in the merged dataframe was mapped to the proper identifier
    # When a column is added, we will make sure that data is not lost.
    # When values are imputed ('Wildtype' for NaN in somatic), we will make sure it is done only for the correct samples.

    print("Evaluating utilities v2...")

    PASS = True

    # Load our dataframes
    proteomics = en.get_proteomics()
    phosphoproteomics = en.get_phosphoproteomics()
    transcriptomics = en.get_transcriptomics()
    cna = en.get_CNA()
    somatic = en.get_somatic()

    # Test compare_gene, using the A1BG gene
    gene = 'A1BG'
    A1BG_compared = en.compare_gene(proteomics, transcriptomics, gene)
    A1BG_compared_name = 'A1BG_compared'

    ### Check the proteomics column
    if not test_merged_column(proteomics, A1BG_compared, gene, A1BG_compared.columns.values[0], A1BG_compared_name):
       PASS = False

    ### Check the transcriptomics column
    if not test_merged_column(transcriptomics, A1BG_compared, gene, A1BG_compared.columns.values[1], A1BG_compared_name):
       PASS = False

    # Test compare_gene, using a list of genes
    gene_list = ['A1BG', 'ZZEF1', 'SMURF1']
    sorted_gene_list = sorted(gene_list) # compare_gene should sort the genes
    list_compared = en.compare_gene(proteomics, transcriptomics, gene_list)
    list_compared_name = 'list_compared'

    ### Test the data from the first dataframe, which are in the first three columns of the merged dataframe
    for i in range(3):
        if not test_merged_column(proteomics, list_compared, sorted_gene_list[i], list_compared.columns.values[i], list_compared_name):
           PASS = False

    ### Test the data from the second dataframe, which are in the last three columns of the merged dataframe
    for i in range(3):
        if not test_merged_column(transcriptomics, list_compared, sorted_gene_list[i], list_compared.columns.values[i + 3], list_compared_name):
           PASS = False

    # Test compare_mutations, using functionality to compare a gene's omics data to its own somatic mutation data
    gene = 'TP53'
    TP53_mutation_compared = en.compare_mutations(proteomics, gene)
    TP53_mutation_compared_name = 'TP53_mutation_compared'

    ### Test data in 'TP53' column, which is the proteomics data for TP53
    if not test_merged_column(proteomics, TP53_mutation_compared, gene, gene, TP53_mutation_compared_name):
        PASS = False

    ### Test data in 'Mutation' column, which is from the somatic mutation data for TP53
    somatic_ID_column = 'Clinical_Patient_Key'
    somatic_filter_column = 'Gene'
    somatic_source_column = 'Mutation'
    somatic_dest_column = 'Mutation'
    if not test_merged_column_from_row(somatic, TP53_mutation_compared, somatic_ID_column, somatic_filter_column, gene, somatic_source_column, somatic_dest_column, TP53_mutation_compared_name):
        PASS = False

    ### Test data in 'Sample_Status' column, which should be 'Tumor' for all samples up to and including S100, and 'Normal' for the remaining ones
    TP53_sample_statuses =  TP53_mutation_compared['Sample_Status']
    for row in TP53_sample_statuses.iteritems():
        sample = row[0]
        value = row[1]
        if sample <= 'S104':
            if value != 'Tumor':
                print('Merged dataframe had incorrect value.\n\tDataframe: {}\n\tSample: {}\tColumn: {}\n\tExpected: {}\tActual: {}\n'.format(TP53_mutation_compared_name, sample, 'Sample_Status', 'Tumor', value))
                PASS = False
        else:
            if value != 'Normal':
                print('Merged dataframe had incorrect value.\n\tDataframe: {}\n\tSample: {}\tColumn: {}\n\tExpected: {}\tActual: {}\n'.format(TP53_mutation_compared_name, sample, 'Sample_Status', 'Normal', value))
                PASS = False

    # Test compare_mutations, using functionality to compare a gene's omcis data to the mutation data for another gene
    gene_2 = 'AURKA'
    multiple_mutation_compared = en.compare_mutations(proteomics, gene, gene_2)

    ### Test data in 'TP53' column, which is the proteomics data for TP53

    ### Test data in 'Mutation' column, which is from the somatic mutation data for AURKA

    ### Test data in 'Sample_Status' column, which is from the somatic mutation data for AURKA

    # Test compare_mutations (phosphosproteomics)

    # Test compare_mutations (phosphoproteomics with somatic)

    # Test compare_mutations_full (proteomics)

    # Test compare_mutations_full (proetomics with somatic)

    # Test compare_mutations_full (phosphosproteomics)

    # Test compare_mutations_full (phosphoproteomics with somatic)

    # Test compare_clinical

    # Test compare_phosphosites

    # Indicate whether the overall test passed
    if PASS:
        print("PASS")
    else:
        print("FAIL")

class Stats:
    def __init__(self):
        pass
    def evaluate(data, trait):
        data_trait = en.compare_clinical(data, trait)
        threshold = .05 / len(data.columns)
        tscutoff = .5
        significantTests = []
        significantGenes = []
        for num in range(1,len(data_trait.columns)):
            gene = data_trait.columns[num]
            oneGene = data_trait[[trait, gene]]
            oneGene = oneGene.dropna(axis=0)
            spearmanrTest = stats.spearmanr(oneGene[trait], oneGene[gene])
            if (abs(spearmanrTest[0]) >= tscutoff) and (spearmanrTest[1] <= threshold):
                significantTests.append(spearmanrTest)
                significantGenes.append(gene)
        if len(significantGenes) > 0:
            return "PASS"
        else:
            return "FAIL"
class Plotter:
    def __init__(self):
        pass
    def plot(data, column1, column2, method):
        if method == "scatterplot":
            plot = sns.relplot(x = column1, y = column2, data = data)
        elif method == "barplot":
            plot = sns.barplot(x = column1, y = column2, data = data)
        elif method == "boxplot":
            plot = sns.boxplot(x = column1, y = column2, data = data)
        else:
            message = method + " not a recognized method"
            print(message)
            return ""
        plt.show()

print("\nRunning tests:\n")

print("Testing getters...")
test_get_clinical_filtered()
test_get_clinical_unfiltered()
test_get_derived_molecular_filtered()
test_get_derived_molecular_unfiltered()
test_get_acetylproteomics_filtered()
test_get_acetylproteomics_unfiltered()
test_get_proteomics()
test_get_transcriptomics_linear()
test_get_transcriptomics_circular()
test_get_transcriptomics_miRNA()
test_get_transcriptomics_with_invalid()
test_get_CNA()
test_get_phosphoproteomics_site()
test_get_phosphoproteomics_gene()
test_get_phosphosites()
test_get_somatic_maf()
test_get_somatic_binary()
test_get_somatic_unparsed()

#evaluate_special_getters()
#evaluate_utilities()
#evaluate_utilities_v2()

# The below tests are not so necessary anymore, now that we have better tests above.

#print("Plotting...")
#Plotter().plot(en.get_proteomics(), "A1BG","PTEN","scatterplot")
#Plotter().plot(en.get_clinical(), "Diabetes","BMI","barplot")
#Plotter().plot(en.get_clinical(), "Diabetes","BMI","boxplot")
#print("PASS")

#print("Running statistics...")
#message = Stats().evaluate(en.get_proteomics(), "Tumor_Size_cm")
#print(message)

print("Version:",en.version())
