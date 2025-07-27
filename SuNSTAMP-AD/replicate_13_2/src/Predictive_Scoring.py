import pandas as pd
import numpy as np
import sys
from sklearn.metrics import roc_auc_score
from scipy.spatial import KDTree

top_n_genes = int(sys.argv[1])  # takes top gene count from command line
output_file = sys.argv[2]       # takes CSV filename

# === Config ===
months = 13
replicate_now = 2
genes_from_replicate = 1
which_rank_from = 10
which_lambda_from = 3e-05

# === File paths ===
gene_expr_file = f"../data/gene_exp_{months}_{replicate_now}.csv"
spatial_file = f"../data/supervision_{months}_{replicate_now}.csv"
gene_score_file = f"../../replicate_{months}_{genes_from_replicate}/result/H_matrices_{months}_{genes_from_replicate}/output_H_matrix_first_row_rank_{which_rank_from}_{which_lambda_from}_{months}_{genes_from_replicate}.csv"
plaque_file = f"../data/plaque_{months}_{replicate_now}.csv"

# === Load data ===
gene_exp_df = pd.read_csv(gene_expr_file)
spatial_df = pd.read_csv(spatial_file)
gene_score_df = pd.read_csv(gene_score_file)
plaque_df = pd.read_csv(plaque_file)

# === Filter for common cells ===
common_cells = list(set(gene_exp_df['Cell_Id']) & set(spatial_df['Cell_Id']))
gene_exp_df = gene_exp_df[gene_exp_df['Cell_Id'].isin(common_cells)].set_index('Cell_Id')
spatial_df = spatial_df[spatial_df['Cell_Id'].isin(common_cells)].set_index('Cell_Id')

# === Get top-N gene scores ===
gene_score_df = gene_score_df.sort_values(by='Value', ascending=False).head(top_n_genes).set_index('Gene')
selected_genes = [g for g in gene_score_df.index if g in gene_exp_df.columns]
if not selected_genes:
    raise ValueError("No selected genes found in gene expression matrix!")

# === Weighted expression vector ===
gene_weights = gene_score_df.loc[selected_genes, 'Value']
gene_weights = gene_weights / gene_weights.sum()
X = gene_exp_df[selected_genes]
predictive_score = (X * gene_weights).sum(axis=1)

# === AUC computation ===
score_df = pd.DataFrame({
    'X': spatial_df['X-scaled'],
    'Y': spatial_df['Y-scaled'],
    'Score': predictive_score
})

# Label cells near plaques
plaque_tree = KDTree(plaque_df[['m.cx', 'm.cy']].values)
score_df['near_plaque'] = [
    len(plaque_tree.query_ball_point(coord, r=30)) > 0
    for coord in score_df[['X', 'Y']].values
]
auc = roc_auc_score(score_df['near_plaque'], score_df['Score'])

# Append result to CSV
with open(output_file, 'a') as f:
    f.write(f"{top_n_genes},{auc:.5f}\n")
