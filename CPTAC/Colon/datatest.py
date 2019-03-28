import CPTAC.Colon as co

prot = co.get_proteomics()
phos = co.get_proteomics()

gene = 'TP53'
gene2 = 'AURKA'

compareds = []

compareds.append(co.compare_mutations(prot, gene))
compareds.append(co.compare_mutations(prot, gene2))

compareds.append(co.compare_mutations(prot, gene, gene2))
compareds.append(co.compare_mutations(prot, gene2, gene))

compareds.append(co.compare_mutations_full(prot, gene))
compareds.append(co.compare_mutations_full(prot, gene2))

compareds.append(co.compare_mutations_full(prot, gene, gene2))
compareds.append(co.compare_mutations_full(prot, gene2, gene))

compareds.append(co.compare_mutations(phos, gene))
compareds.append(co.compare_mutations(phos, gene2))

compareds.append(co.compare_mutations(phos, gene, gene2))
compareds.append(co.compare_mutations(phos, gene2, gene))

compareds.append(co.compare_mutations_full(phos, gene))
compareds.append(co.compare_mutations_full(phos, gene2))

compareds.append(co.compare_mutations_full(phos, gene, gene2))
compareds.append(co.compare_mutations_full(phos, gene2, gene))

# Test other utilities
compareds.append(co.compare_clinical(prot, 'Subsite'))
compareds.append(co.compare_phosphosites(gene))
compareds.append(co.compare_gene(prot, phos, gene))

for df in compareds:
    print(df)

print(len(compareds))

