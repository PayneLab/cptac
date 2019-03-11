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
        dataframes.append(en.get_transcriptomics(data_type="circular")); file_names[len(dataframes)] = "circular transcriptomics"
        dataframes.append(en.get_transcriptomics(data_type="miRNA")); file_names[len(dataframes)] = "miRNA"
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

    def check_dataframe(self, name, df, exp_dim, exp_headers, coordinates, values):
        """
        Parameters
        name: string containing the name of the dataframe gotten by the getter we're testing
        df: the dataframe gotten by the getter we are testing
        exp_dim: a tuple containing the expected dimensions of the dataframe, in the format (rows, columns)
        exp_headers: if the dataframe has up to 20 columns, all of the headers for the dataframe, in order. If it has more than 20 columns, then a list containing the first ten and last ten headers, in order.
        coordinates: a tuple with three elements, each element being a tuple with two elements, the first element being the int index of the row of a test value, and the second element being the int index of the column of a test value
        values: a tuple with three elements, each element being the expected value of the test value corresponding to the coordinates at the same index in the coordinates parameter 

        Returns
        Bool indicating if dataframe passed the test
        """
        PASS = True

        df = df

        ## Check dimensions
        act_dim = df.shape
        if exp_dim != act_dim:
            print("Error: {} dataframe dimensions did not match expected values.\n\tExpected: {}\n\tActual: {}".format(name, exp_dim, act_dim))
            PASS = False

        ## Check headers
        act_headers_all = list(df.columns.values)
        if len(df.columns.values) <= 20:
            act_headers = act_headers_all
        else:
            act_headers = act_headers_all[:10] + act_headers_all[-10:]

        if len(exp_headers) != len(act_headers):
            print("Unexpected number of test headers in {} dataframe. Expected number of headers: {}. You passed {} headers.".format(name, len(act_headers), len(exp_headers)))
            PASS = False
        else:
            for i, header in enumerate(exp_headers):
                if header != act_headers[i]:
                    print("Error: {} dataframe header did not match expected value.\n\tExpected: {}\n\tActual: {}".format(name, header, act_headers[i]))
                    PASS = False

        ## Check test values
        act_values = [
            df.iloc[coordinates[0][0], coordinates[0][1]],
            df.iloc[coordinates[1][0], coordinates[1][1]],
            df.iloc[coordinates[2][0], coordinates[2][1]]
        ]

        for i, value in enumerate(values):
            if act_values[i] != value:
                print("Error: {} dataframe value did not match expected value.\n\tColumn: {}\n\tIndex: {}\n\tExpected: {}\n\tActual: {}".format(name, df.columns.values[coordinates[i][1]], df.index.values[coordinates[i][0]], value, act_values[i]))
                PASS = False

        # The following line mostly to check everything's working during development
        print("\t\tPass: " + str(PASS))

        # Return whether the getter passed the test
        return PASS

    def evaluate_getters_v2(self):
        # ***UNDER CONSTRUCTION***

        print("Evaluating getters v2...")
        tester = Basic()
        PASS = True

        # Test get_clinical() with the default parameter excluded=False
        if not tester.check_dataframe(
            "Clinical",
            en.get_clinical(),
            (144, 171),
            ['Proteomics_Participant_ID', 'Case_excluded', 'Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs', 'Proteomics_Aliquot_ID', 'Proteomics_Tumor_Normal', 'Proteomics_OCT', 'Country', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID', 'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality'],
            ((79, 81), (15, 146), (88, 12)),
            (-1.03, 8.888888889, 'Serous')
        ):
            PASS = False

        # Test get_clinical(excluded=True)
        if not tester.check_dataframe(
            "Clinical (with excluded cases)",
            en.get_clinical(excluded=True),
            (153, 171),
            ['Proteomics_Participant_ID', 'Case_excluded', 'Proteomics_TMT_batch', 'Proteomics_TMT_plex', 'Proteomics_TMT_channel', 'Proteomics_Parent_Sample_IDs', 'Proteomics_Aliquot_ID', 'Proteomics_Tumor_Normal', 'Proteomics_OCT', 'Country', 'RNAseq_R1_sample_type', 'RNAseq_R1_filename', 'RNAseq_R1_UUID', 'RNAseq_R2_sample_type', 'RNAseq_R2_filename', 'RNAseq_R2_UUID', 'miRNAseq_sample_type', 'miRNAseq_UUID', 'Methylation_available', 'Methylation_quality'],
            ((23, 44), (151, 6), (32, 165)),
            (0.004118258, 'CPT0230400002,CPT0230400003,CPT0230400004,CPT0230410002,CPT0230410003,CPT0230410004,CPT0230420002,CPT0230420003,CPT0230420004', '171011_UNC31-K00269_0086_AHLJLCBBXX_CTTGTA_S7_L003_R2_001.fastq.gz')
        ):
            PASS = False

        # Test get_proteomics()
        if not  tester.check_dataframe(
            "Proteomics",
            en.get_proteomics(),
            (153, 10999),
            ['A1BG', 'A2M', 'A2ML1', 'A4GALT', 'AAAS', 'AACS', 'AADAT', 'AAED1', 'AAGAB', 'AAK1', 'ZSWIM8', 'ZSWIM9', 'ZW10', 'ZWILCH', 'ZWINT', 'ZXDC', 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3'],
            ((34, 6003), (99, 9544), (152, 32)),
            (-0.8170000000000001, -1.28, 0.904)
        ):
            PASS = False

        # Test get_transcriptomics() with default parameter data_type="linear"
        if not tester.check_dataframe(
            "Transcriptomics (linear)",
            en.get_transcriptomics(),
            (115, 28057),
            ['A1BG', 'A1BG-AS1', 'A1CF', 'A2M', 'A2M-AS1', 'A2ML1', 'A2MP1', 'A3GALT2', 'A4GALT', 'A4GNT', 'ZWILCH', 'ZWINT', 'ZXDA', 'ZXDB', 'ZXDC', 'ZYG11A', 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3'],
            ((22, 25483), (110, 23), (101, 17748)),
            (0.89, 11.83, 7.02)
        ):
            PASS = False

        # Test get_transcriptomics(data_type="circular")
        if not tester.check_dataframe(
            "Transcriptomics (circular)",
            en.get_transcriptomics(data_type="circular"),
            (115, 4945),
            ['circ_chr10_100260218_100262063_CWF19L1', 'circ_chr10_100923975_100926019_SLF2', 'circ_chr10_100923978_100926019_SLF2', 'circ_chr10_100937402_100944128_SLF2', 'circ_chr10_100937402_100950753_SLF2', 'circ_chr10_101584602_101586156_POLL', 'circ_chr10_101667886_101676436_FBXW4', 'circ_chr10_101672915_101676436_FBXW4', 'circ_chr10_101792839_101807901_OGA', 'circ_chr10_101792839_101810314_OGA', 'circ_chrX_80288906_80310233_CHMP1B2P', 'circ_chrX_80289664_80310233_CHMP1B2P', 'circ_chrX_80707427_80719656_BRWD3', 'circ_chrX_80791854_80793772_BRWD3', 'circ_chrX_84096194_84164387_RPS6KA6', 'circ_chrX_84134782_84164387_RPS6KA6', 'circ_chrX_85067127_85074391_APOOL', 'circ_chrX_85978767_85981809_CHM', 'circ_chrX_91414904_91418871_PABPC5-AS1', 'circ_chrX_9691579_9693419_TBL1X'],
            ((110, 1), (34, 4935), (73, 2003)),
            (8.85, 6.48, 0.0)
        ):
            PASS = False

        # Test get_transcriptomics(data_type="miRNA")
        if not tester.check_dataframe(
            "Transcriptomics (miRNA)",
            en.get_transcriptomics(data_type="miRNA"),
            (104, 2337),
            ['hsa-let-7a-2-3p', 'hsa-let-7a-3p', 'hsa-let-7a-5p', 'hsa-let-7b-3p', 'hsa-let-7b-5p', 'hsa-let-7c-3p', 'hsa-let-7c-5p', 'hsa-let-7d-3p', 'hsa-let-7d-5p', 'hsa-let-7e-3p', 'hsa-miR-9901', 'hsa-miR-9902', 'hsa-miR-9903', 'hsa-miR-9983-3p', 'hsa-miR-9985', 'hsa-miR-9986', 'hsa-miR-99a-3p', 'hsa-miR-99a-5p', 'hsa-miR-99b-3p', 'hsa-miR-99b-5p'],
            ((5, 0), (100, 1597), (54, 2231)),
            (1.79, 1.25, 1.86)
        ):
            PASS = False

        # Test get_transcriptomics() with an invalid parameter, and make sure that it raises an exception
        try:
            en.get_transcriptomics("lorem ipsum")
        except ValueError:
            pass
        else:
            print("Error: get_transcriptomics() did not raise ValueError as expected, when given invalid parameter.")
            PASS = False

        # Test get_CNA()
        if not tester.check_dataframe(
            "CNA",
            en.get_CNA(),
            (103, 28057),
            ['MFSD14A', 'SASS6', 'TRMT13', 'LRRC39', 'DBT', 'RTCA-AS1', 'RTCA', 'MIR553', 'UBE4B', 'CDC14A', 'TSPY8', 'FAM197Y2', 'FAM197Y4', 'FAM197Y5', 'FAM197Y7', 'FAM197Y8', 'FAM197Y6', 'FAM197Y3', 'RBMY3AP', 'TTTY22'],
            ((12, 27865), (67, 8), (102, 15439)),
            (-0.19, 0.01, 0.03)
        ):
            PASS = False

        # Test get_phosphoproteomics() with default parameter gene_level=False
        phosphoproteomics_site_name = "Phosphoproteomics (site)"
        phosphoproteomics_site_df =  en.get_phosphoproteomics()
        phosphoproteomics_site_dim = (153, 73212)
        phosphoproteomics_site_headers = ['AAAS-S495', 'AAAS-S541', 'AAAS-Y485', 'AACS-S618', 'AAED1-S12', 'AAGAB-S310', 'AAGAB-S311', 'AAK1-S14', 'AAK1-S18', 'AAK1-S20', 'ZZZ3-S397', 'ZZZ3-S411', 'ZZZ3-S420', 'ZZZ3-S424', 'ZZZ3-S426', 'ZZZ3-S468', 'ZZZ3-S89', 'ZZZ3-T415', 'ZZZ3-T418', 'ZZZ3-Y399']
        phosphoproteomics_site_test_coord = ((46, 45), (12, 72435), (96, 45362))
        phosphoproteomics_site_test_vals = (0.195, -0.27899999999999997, -0.13)

        if not tester.check_dataframe(phosphoproteomics_site_name, phosphoproteomics_site_df, phosphoproteomics_site_dim, phosphoproteomics_site_headers, phosphoproteomics_site_test_coord, phosphoproteomics_site_test_vals):
            PASS = False

        # Test get_phosphoproteomics(gene_level=True)
        phosphoproteomics_gene_name = "Phosphoproteomics (gene)"
        phosphoproteomics_gene_df = en.get_phosphoproteomics(gene_level=True)
        phosphoproteomics_gene_dim = (153, 8466)
        phosphoproteomics_gene_headers = ['AAAS', 'AACS', 'AAED1', 'AAGAB', 'AAK1', 'AAMDC', 'AARS', 'AASDH', 'AATF', 'ABCA1', 'ZSCAN5C', 'ZSWIM3', 'ZSWIM8', 'ZUP1', 'ZW10', 'ZXDA', 'ZXDC', 'ZYX', 'ZZEF1', 'ZZZ3']
        phosphoproteomics_gene_test_coord =  ((2, 7999), (148, 1045), (78, 6543))
        phosphoproteomics_gene_test_vals = (-0.0879, 1.33, 0.153)

        if not tester.check_dataframe(phosphoproteomics_gene_name, phosphoproteomics_gene_df, phosphoproteomics_gene_dim, phosphoproteomics_gene_headers, phosphoproteomics_gene_test_coord, phosphoproteomics_gene_test_vals):
            PASS = False

        # Test get_phosphosites
        phosphosites_name = "Phosphosites for the AAK1-S14 gene"
        phosphosites_df = en.get_phosphosites('AAK1-S14')
        phosphosites_dim = (153, 1)
        phosphosites_headers = ['AAK1-S14']
        phosphosites_test_coord = ((34, 0), (76, 0), (143, 0))
        phosphosites_test_vals = (1.46, -0.0511, 0.9940000000000001)

        if not tester.check_dataframe(phosphosites_name, phosphosites_df, phosphosites_dim, phosphosites_headers, phosphosites_test_coord, phosphosites_test_vals):
            PASS = False

        # Test get_somatic() with default parameters binary=False, unparsed=False (this will return the Somatic Maf dataframe)
        somatic_name = "Somatic (maf)"
        somatic_df = en.get_somatic()
        somatic_dim = (53101, 5)
        somatic_headers = ['Patient_Id', 'Gene', 'Mutation', 'Location', 'Clinical_Patient_Key']
        somatic_test_coord = ((53000, 3), (12, 4), (34567, 0))
        somatic_test_vals = ('p.M1259I', 'S001', 'C3N-00151')

        if not tester.check_dataframe(somatic_name, somatic_df, somatic_dim, somatic_headers, somatic_test_coord, somatic_test_vals):
            PASS = False

        # Test get_somatic(binary=True) and therefore unparsed=False
        somatic_binary_name = "Somatic (binary)"
        somatic_binary_df = en.get_somatic(binary=True)
        somatic_binary_dim = (103, 51559)
        somatic_binary_headers = ['A1BG_p.E298K', 'A1BG_p.S181N', 'A1CF_p.F487L', 'A1CF_p.S236Y', 'A2ML1_p.A8V', 'A2ML1_p.G1306D', 'A2ML1_p.L1347F', 'A2ML1_p.L82I', 'A2ML1_p.P712S', 'A2ML1_p.R443Q', 'ZYG11A_p.Q442H', 'ZYG11B_p.H315R', 'ZYG11B_p.R495M', 'ZYG11B_p.R728C', 'ZYX_p.C447Y', 'ZZEF1_p.A2723V', 'ZZEF1_p.D845Y', 'ZZEF1_p.K1251E', 'ZZEF1_p.K2387Sfs*40', 'ZZZ3_p.Y891C']
        somatic_binary_test_coord = ((102, 51558), (0, 0), (45, 25436))
        somatic_binary_test_vals = (0, 0, 0)

        if not tester.check_dataframe(somatic_binary_name, somatic_binary_df, somatic_binary_dim, somatic_binary_headers, somatic_binary_test_coord, somatic_binary_test_vals):
            PASS = False

        
        # Test get_somatic(unparsed=True) and therefore binary=False
        somatic_unparsed_name = "Somatic (unparsed)"
        somatic_unparsed_df = en.get_somatic(unparsed=True)
        somatic_unparsed_dim = (53101, 124)
        somatic_unparsed_headers = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'ExAC_AC_AN_Adj', 'ExAC_AC_AN', 'ExAC_AC_AN_AFR', 'ExAC_AC_AN_AMR', 'ExAC_AC_AN_EAS', 'ExAC_AC_AN_FIN', 'ExAC_AC_AN_NFE', 'ExAC_AC_AN_OTH', 'ExAC_AC_AN_SAS', 'ExAC_FILTER']
        somatic_unparsed_test_coord = ((52265, 45), (12, 70), (27658, 1))
        somatic_unparsed_test_vals = ('strelkasnv-varssnv-mutectsnv', 'UPI0000167B91', 0)

        if not tester.check_dataframe(somatic_unparsed_name, somatic_unparsed_df, somatic_unparsed_dim, somatic_unparsed_headers, somatic_unparsed_test_coord, somatic_unparsed_test_vals):
            PASS = False

        # Indicate whether the overall test passed
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

    def evaluate_utilities_v2(self):
        # We will test all of the compare_**** functions, which either merge dataframes, or add columns to a dataframe
        # When dataframes are merged, we will make sure that the data in the merged dataframe was mapped to the proper identifier
        # When a column is added, we will make sure that data is not lost.

        # Load our dataframes
        proteomics = en.get_proteomics()
        phosphoproteomics = en.get_phosphoproteomics()
        transcriptomics = en.get_transcriptomics()
        cna = en.get_CNA()

        # Test compare_gene, using the A1BG gene
        compared_A1BG = en.compare_gene(proteomics, transcriptomics, 'A1BG')
        # Loop through each row in the merged dataframe, and make sure that the values for each sample match the values for that sample in the original two dataframes.

        # Test compare_genes

        # Test compare_mutations (proteomics)

        # Test compare_mutations (proetomics with somatic)

        # Test compare_mutations (phosphosproteomics)

        # Test compare_mutations (phosphoproteomics with somatic)

        # Test compare_mutations_full (proteomics)

        # Test compare_mutations_full (proetomics with somatic)

        # Test compare_mutations_full (phosphosproteomics)

        # Test compare_mutations_full (phosphoproteomics with somatic)

        # Test compare_clinical

        # Test compare_phosphosites

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
Basic().evaluate_getters_v2()

print("Plotting...")
Plotter().plot(en.get_proteomics(), "A1BG","PTEN","scatterplot")
Plotter().plot(en.get_clinical(), "Diabetes","BMI","barplot")
Plotter().plot(en.get_clinical(), "Diabetes","BMI","boxplot")
print("PASS")

print("Running statistics...")
message = Stats().evaluate(en.get_proteomics(), "Pathway_activity_p53")
print(message)

print("Version:",en.version())
