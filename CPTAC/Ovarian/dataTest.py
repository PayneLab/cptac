import pandas as pd
import CPTAC.Ovarian as ov

class Basic:
    def __init__(self):
        pass
    def evaluate_getters(self):
        print("Evaluating getters...")
        dataframes = []
        file_names = {}
        dataframes.append(ov.get_clinical()); file_names[len(dataframes)] = "clinical"
        dataframes.append(ov.get_cnv()); file_names[len(dataframes)] = "cnv"
        dataframes.append(ov.get_phosphoproteomics()); file_names[len(dataframes)] = "phosphoproteomics"
        dataframes.append(ov.get_proteomics()); file_names[len(dataframes)] = "proteomics"
        dataframes.append(ov.get_somatic_mutations()); file_names[len(dataframes)] = "somatic_19"
        dataframes.append(ov.get_transcriptomics()); file_names[len(dataframes)] = "transcripmics"
        PASS = True
        for x in range(0,len(dataframes)):
            if dataframes[x] is None:
                print("Error reading", file_names[x+1], "data")
                PASS = False
        if PASS:
            print("PASS")
        else:
            print("FAIL")


print("\nRunning tests:\n")

Basic().evaluate_getters()
