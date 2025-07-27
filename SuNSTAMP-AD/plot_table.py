import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations

# Load both gene lists
rep1 = pd.read_csv('./replicate_13_1/result/final_gene_list_13_1.csv')  # columns: 'Gene', 'Value'
rep2 = pd.read_csv('./replicate_13_2/result/final_gene_list_13_2.csv')

# Take top 20 genes from each
top20_1 = rep1.sort_values('Value', ascending=False).head(20)['Gene']
top20_2 = rep2.sort_values('Value', ascending=False).head(20)['Gene']

# Common genes and rankings
set1 = set(top20_1)
set2 = set(top20_2)
common_genes = sorted(list(set1 & set2))  # fixed order

rank1 = {gene: rank for rank, gene in enumerate(top20_1)}
rank2 = {gene: rank for rank, gene in enumerate(top20_2)}

# Initialize consistency matrix (fill with np.nan)
n = len(common_genes)
matrix = np.full((n, n), np.nan)

# Fill lower triangle with 1 (consistent) or 0 (inconsistent)
red = 0
green = 0
for i in range(n):
    for j in range(i):
        gene_i = common_genes[i]
        gene_j = common_genes[j]
        order1 = rank1[gene_i] - rank1[gene_j]
        order2 = rank2[gene_i] - rank2[gene_j]
        if order1 * order2 > 0:
            matrix[i, j] = 1
            green += 1
        else:
            matrix[i, j] = 0
            red += 1

# Convert to DataFrame for seaborn
df_matrix = pd.DataFrame(matrix, index=common_genes, columns=common_genes)

# Plot settings
plt.figure(figsize=(10, 8))
soothing_cmap = sns.color_palette(["#b0cbe9", "#2b8cbe"])  # light blue to mid-blue

sns.heatmap(
    df_matrix,
    cmap=soothing_cmap,
    linewidths=0.5,
    linecolor='gray',
    square=True,
    cbar=False,
    annot=True,
    fmt='.0f',
    mask=np.triu(np.ones_like(df_matrix, dtype=bool)),
    annot_kws={"fontsize": 11, "weight": "bold"}
)

# Customize ticks and title
plt.title('Pairwise Ordering Consistency Between Common Genes', fontsize=16, weight='bold')
plt.xticks(rotation=45, ha='right', fontsize=12, weight='bold')
plt.yticks(rotation=0, fontsize=12, weight='bold')

plt.tight_layout()
plt.savefig("table_genes.png", dpi=300)
plt.show()

# Stats summary
print(f"red (inconsistent): {red}")
print(f"green (consistent): {green}")
total = red + green
print(f"total pairs: {total}")
consistency = green / total
print(f"percentage consistent: {consistency:.2%}")
