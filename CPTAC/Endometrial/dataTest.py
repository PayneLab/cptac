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

class Basic:
    def __init__(self):
        pass
    def evaluate_getters(self):
        print("Evaluating getters...")
        dataframes = []
        file_names = {}
        dataframes.append(en.get_clinical()); file_names[len(dataframes)] = "clinical"
        dataframes.append(en.get_proteomics()); file_names[len(dataframes)] = "proteomics"
        dataframes.append(en.get_phosphoproteomics()); file_names[len(dataframes)] = "site phosphoproteomics"
        dataframes.append(en.get_phosphoproteomics(gene_level = True)); file_names[len(dataframes)] = "gene phosphoproteomics"
        dataframes.append(en.get_transcriptomics()); file_names[len(dataframes)] = "linear transcriptomics"
        dataframes.append(en.get_transcriptomics(circular = True)); file_names[len(dataframes)] = "circular transcriptomics"
        dataframes.append(en.get_transcriptomics(miRNA = True)); file_names[len(dataframes)] = "miRNA"
        dataframes.append(en.get_CNA()); file_names[len(dataframes)] = "CNA"
        dataframes.append(en.get_somatic()); file_names[len(dataframes)] = "parsed somatic maf"
        dataframes.append(en.get_somatic(binary=True)); file_names[len(dataframes)] = "binary somatic"
        dataframes.append(en.get_somatic(unparsed=True)); file_names[len(dataframes)] = "unparsed somatic maf"
        PASS = True
        for x in range(0,len(dataframes)):
            if dataframes[x] is None:
                print("Error reading", file_names[x+1], "data")
                PASS = False
        if PASS:
            print("PASS")
        else:
            print("FAIL")

    def evaluate_getters_2(self):
        # ***UNDER CONSTRUCTION***

        # We will call each get function, and test the dimensions, headers, and some test values for it.
        print("Evaluating getters 2...")
        PASS = True

        # Test get_clinical() with default excluded=False
        clinical = en.get_clinical()

        ## Check dimensions
        exp_clinical_dim = (144, 170)
        act_clinical_dim = clinical.shape
        if exp_clinical_dim != act_clinical_dim:
            print("Error: Clinical dataframe dimensions did not match expected values.\nExpected: {}\nActual: {}".format(exp_clinical_dim, act_clinical_dim))
            PASS = False

        ## Check headers
        exp_clinical_headers = ['Proteomics_Participant_ID', 'Case_excluded', 'Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs', 'Proteomics_Aliquot_ID', 'Proteomics_Tumor_Normal', 'Proteomics_OCT', 'Country', 'Histologic_Grade_FIGO', 'Myometrial_invasion_Specify', 'Histologic_type', 'Treatment_naive', 'Tumor_purity', 'Path_Stage_Primary_Tumor-pT', 'Path_Stage_Reg_Lymph_Nodes-pN', 'Clin_Stage_Dist_Mets-cM', 'Path_Stage_Dist_Mets-pM', 'tumor_Stage-Pathological', 'FIGO_stage', 'LVSI', 'BMI', 'Age', 'Diabetes', 'Race', 'Ethnicity', 'Gender', 'Tumor_Site', 'Tumor_Site_Other', 'Tumor_Focality', 'Tumor_Size_cm', 'Estrogen_Receptor', 'Estrogen_Receptor_%', 'Progesterone_Receptor', 'Progesterone_Receptor_%', 'MLH1', 'MLH2', 'MSH6', 'PMS2', 'p53', 'Other_IHC_specify', 'MLH1_Promoter_Hypermethylation', 'Num_full_term_pregnancies', 'EPIC_Bcells', 'EPIC_CAFs', 'EPIC_CD4_Tcells', 'EPIC_CD8_Tcells', 'EPIC_Endothelial', 'EPIC_Macrophages', 'EPIC_NKcells', 'EPIC_otherCells', 'CIBERSORT_B _cells _naive', 'CIBERSORT_B _cells _memory', 'CIBERSORT_Plasma _cells', 'CIBERSORT_T _cells _CD8', 'CIBERSORT_T _cells _CD4 _naive', 'CIBERSORT_T _cells _CD4 _memory _resting', 'CIBERSORT_T _cells _CD4 _memory _activated', 'CIBERSORT_T _cells _follicular _helper', 'CIBERSORT_T _cells _regulatory _(Tregs)', 'CIBERSORT_T _cells _gamma _delta', 'CIBERSORT_NK _cells _resting', 'CIBERSORT_NK _cells _activated', 'CIBERSORT_Monocytes', 'CIBERSORT_Macrophages _M0', 'CIBERSORT_Macrophages _M1', 'CIBERSORT_Macrophages _M2', 'CIBERSORT_Dendritic _cells _resting', 'CIBERSORT_Dendritic _cells _activated', 'CIBERSORT_Mast _cells _resting', 'CIBERSORT_Mast _cells _activated', 'CIBERSORT_Eosinophils', 'CIBERSORT_Neutrophils', 'CIBERSORT_Absolute _score', 'ESTIMATE_StromalScore', 'ESTIMATE_ImmuneScore', 'ESTIMATE_ESTIMATEScore', 'Stemness_score', 'ER_ESR1', 'PR_PGR', 'Pathway_activity_EGFR', 'Pathway_activity_Hypoxia', 'Pathway_activity_JAK.STAT', 'Pathway_activity_MAPK', 'Pathway_activity_NFkB', 'Pathway_activity_PI3K', 'Pathway_activity_TGFb', 'Pathway_activity_TNFa', 'Pathway_activity_Trail', 'Pathway_activity_VEGF', 'Pathway_activity_p53', 'TP53_ATM', 'TP53_CHEK2', 'TP53_MDM4', 'TP53_RPS6KA3', 'TP53_TP53', 'TP53_pathway', 'PI3K_AKT1', 'PI3K_AKT2', 'PI3K_AKT3', 'PI3K_DEPDC5', 'PI3K_DEPTOR', 'PI3K_INPP4B', 'PI3K_MAPKAP1', 'PI3K_MLST8', 'PI3K_MTOR', 'PI3K_NPRL2', 'PI3K_NPRL3', 'PI3K_PDK1', 'PI3K_PIK3CA', 'PI3K_PIK3CB', 'PI3K_PIK3R1', 'PI3K_PIK3R2', 'PI3K_PPP2R1A', 'PI3K_PTEN', 'PI3K_RHEB', 'PI3K_RICTOR', 'PI3K_RPS6', 'PI3K_RPS6KB1', 'PI3K_RPTOR', 'PI3K_STK11', 'PI3K_TSC1', 'PI3K_TSC2', 'PI3K_pathway', 'HRD_BRCA1', 'HRD_BRCA2', 'HRD_BRCA1_or_BRCA2', 'CNV_clustering', 'CNV_1q_amplification', 'CNV_index', 'Purity_Immune', 'Purity_Cancer', 'Purity_Stroma', 'MSI_status', 'POLE_subtype', 'JAK1_MS_INDEL', 'JAK1_Mutation', 'Log2_variant_per_Mbp', 'Log2_SNP_per_Mbp', 'Log2_INDEL_per_Mbp', 'Log2_variant_total', 'Log2_SNP_total', 'Log2_INDEL_total', 'Mutation_signature_C>A', 'Mutation_signature_C>G', 'Mutation_signature_C>T', 'Mutation_signature_T>C', 'Mutation_signature_T>A', 'Mutation_signature_T>G', 'WXS_normal_sample_type', 'WXS_normal_filename', 'WXS_normal_UUID', 'WXS_tumor_sample_type', 'WXS_tumor_filename', 'WXS_tumor_UUID', 'WGS_normal_sample_type', 'WGS_normal_UUID', 'WGS_tumor_sample_type', 'WGS_tumor_UUID', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID', 'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality']

        act_clinical_headers = clinical.columns.values

        if len(exp_clinical_headers) != len(act_clinical_headers):
            # We shouldn't get to this point unless the dimension check failed.
            print("Error: Clinical dataframe had unexpected number of headers.\nExpected: {}\nActual: {}".format(len(exp_clinical_headers), len(act_clinical_headers)))
            PASS = False
        else:
            for i, header in enumerate(exp_clinical_headers):
                if header != act_clinical_headers[i]:
                    print("Error: Clinical dataframe header did not match expected value.\nExpected: {}\nActual: {}".format(header, act_clinical_headers[i]))
                    PASS = False

        ## Check some test values
        exp_clinical_values = [] # Actually make this a dict, with the keys the test values, and the values a tuple of the index and column
        act_clinical_values = []

        for i, value in exp_clinical_values:
            if value != act_clinical_values[i]:
# need to fix format                print("Error: Clinical dataframe value for did not match expected value.\nIndex: {}\nColumn: {}\nExpected: {}\nActual: {}".format(value, act_clinical_values[i]))
                PASS = False
                

        # Print whether we passed or failed
        if PASS:
            print("PASS")
        else:
            print("FAIL")
	
    def evaluate_special_getters(self):
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
    def evaluate_utilities(self): #compare_**** functions
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


class Stats:
    def __init__(self):
        pass
    def evaluate(self, data, trait):
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
    def plot(self, data, column1, column2, method):
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

Basic().evaluate_getters()
Basic().evaluate_special_getters()
Basic().evaluate_utilities()
# Basic().evaluate_getters_2()

print("Plotting...")
Plotter().plot(en.get_proteomics(), "A1BG","PTEN","scatterplot")
Plotter().plot(en.get_clinical(), "Diabetes","BMI","barplot")
Plotter().plot(en.get_clinical(), "Diabetes","BMI","boxplot")
print("PASS")

print("Running statistics...")
message = Stats().evaluate(en.get_proteomics(), "Pathway_activity_p53")
print(message)

print("Version:",en.version())
