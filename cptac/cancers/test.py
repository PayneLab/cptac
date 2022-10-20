import cptac

br = cptac.Brca()

#print(br.get_genotype_all_vars("TP53", 'washu', 'washu'))

df1 = br.get_dataframe('proteomics', 'umich')
df2 = br.get_dataframe('transcriptomic', 'washu')

print(br.quick_join('umich proteomics', 'washu transcriptomics'))

