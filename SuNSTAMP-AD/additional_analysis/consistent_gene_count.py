from itertools import combinations


import pandas as pd

# Load both gene lists
rep1 = pd.read_csv('./replicate_13_1/result/final_gene_list_13_1.csv')  # columns: 'Gene', 'Value'
rep2 = pd.read_csv('./replicate_13_2/result/final_gene_list_13_2.csv')

# Take top 20 genes from each
top20_1 = rep1.sort_values('Value', ascending=False).head(20)['Gene']
top20_2 = rep2.sort_values('Value', ascending=False).head(20)['Gene']


set1 = set(top20_1)
set2 = set(top20_2)


# Make gene → rank mappings
rank1 = {gene: rank for rank, gene in enumerate(top20_1)}
rank2 = {gene: rank for rank, gene in enumerate(top20_2)}

# Get common genes
common_genes = list(set1 & set2)

# Count consistent pairs
total = 0
consistent = 0
for gene1, gene2 in combinations(common_genes, 2):
    # Ordering in rep1
    order1 = rank1[gene1] - rank1[gene2]
    order2 = rank2[gene1] - rank2[gene2]
    if order1 * order2 > 0:  # same sign = consistent
        consistent += 1
    total += 1

print(f"Consistent pairs: {consistent}/{total}")
