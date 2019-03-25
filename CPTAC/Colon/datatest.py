import pandas as pd
import CPTAC.Colon as co

class Basic:
    def __init__(self):
        pass

    def evaluate_comparers(self):
        # Load our data
        prot = co.get_proteomics()
        phos = co.get_phosphoproteomics()

        # Set our variables
        gene_1 = 'TP53'
        gene_2 = 'AURKA'

        # Test compare_mutations with one gene
        comp_1 = co.compare_mutations(prot, gene_1)
        comp_2 = co.compare_mutations(phos, gene_1)

        # Test compare_mutations with two genes. Should return omics data for first gene, with somatic mutation data for second gene
        comp_3 = co.compare_mutations(prot, gene_1, gene_2)
        comp_4 = co.compare_mutations(phos, gene_1, gene_2)

        # Test compare_mutations_full with one gene
        comp_5 = co.compare_mutations_full(prot, gene_1)
        comp_6 = co.compare_mutations_full(phos, gene_1)

        # Test compare_mutations_full with two genes. Should return omics data for first gene, with somatic mutation data for second gene
        comp_7 = co.compare_mutations_full(prot, gene_1, gene_2)
        comp_8 = co.compare_mutations_full(phos, gene_1, gene_2)

        print(comp_1)
        print(comp_2)
        print(comp_3)
        print(comp_4)
        print(comp_5)
        print(comp_6)
        print(comp_7)
        print(comp_8)

        print('PASS')

print("\nRunning tests:\n")

Basic().evaluate_comparers()

