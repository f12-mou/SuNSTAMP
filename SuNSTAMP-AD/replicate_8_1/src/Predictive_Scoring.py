import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score
from scipy.spatial import KDTree


months = 8
replicate_now = 1
genes_from_replicate = 2
which_rank_from = 14
which_lambda_from = 1.75e-05

# === Parameters ===
gene_expr_file = f"../data/gene_exp_{months}_{replicate_now}.csv"
spatial_file = f"../data/supervision_{months}_{replicate_now}.csv"
gene_score_file = f"../../replicate_{months}_{genes_from_replicate}/result/H_Matrices_{months}_{genes_from_replicate}/output_H_matrix_first_row_rank_{which_rank_from}_{which_lambda_from}_{months}_{genes_from_replicate}.csv"
plaque_file = f"../data/plaque_{months}_{replicate_now}.csv"
top_n_genes = 20 # Change as needed

# === Step 1: Load data ===
gene_exp_df = pd.read_csv(gene_expr_file)
spatial_df = pd.read_csv(spatial_file)
gene_score_df = pd.read_csv(gene_score_file)
plaque_df = pd.read_csv(plaque_file)

# === Step 2: Find common cells ===
cells_expr = set(gene_exp_df['Cell_Id'])
cells_spatial = set(spatial_df['Cell_Id'])
common_cells = list(cells_expr & cells_spatial)

# === Step 3: Subset both dataframes to common cells ===
gene_exp_df = gene_exp_df[gene_exp_df['Cell_Id'].isin(common_cells)].set_index('Cell_Id')
spatial_df = spatial_df[spatial_df['Cell_Id'].isin(common_cells)].set_index('Cell_Id')

# === Step 4: Select top-N genes from gene score file ===
gene_score_df = gene_score_df.sort_values(by='Value', ascending=False).head(top_n_genes)
gene_score_df = gene_score_df.set_index('Gene')

# Keep only genes present in the gene expression data
selected_genes_present = [g for g in gene_score_df.index if g in gene_exp_df.columns]
if not selected_genes_present:
    raise ValueError("None of the selected genes are found in gene expression data!")

# Subset weights and normalize
gene_weights = gene_score_df.loc[selected_genes_present, 'Value']
gene_weights_normalized = gene_weights / gene_weights.sum()

# Subset expression matrix
X = gene_exp_df[selected_genes_present]

# === Step 5: Compute weighted predictive score ===
predictive_score = (X * gene_weights_normalized).sum(axis=1)

# === Step 6: Merge with Total_Attraction and coordinates ===
score_df = pd.DataFrame({
    'Predictive_Score': predictive_score,
    'Total_Attraction': spatial_df['Total_Attraction'],
    'X': spatial_df['X-scaled'],
    'Y': spatial_df['Y-scaled']
})

# === Step 7: Correlation ===
pearson_r, pearson_p = pearsonr(score_df['Predictive_Score'], score_df['Total_Attraction'])
spearman_r, spearman_p = spearmanr(score_df['Predictive_Score'], score_df['Total_Attraction'])

print(f"\n✅ Pearson correlation:  r = {pearson_r:.4f},  p = {pearson_p:.4e}")
print(f"✅ Spearman correlation: r = {spearman_r:.4f},  p = {spearman_p:.4e}")


# === Step 8: Predictive score heatmap (draw high scores last) ===
score_df = score_df.sort_values(by='Predictive_Score', ascending=True)

# === Step 8: Predictive score heatmap ===
plt.figure(figsize=(8, 6))
sc = plt.scatter(
    score_df['X'], score_df['Y'], 
    c=score_df['Predictive_Score'], 
    cmap='viridis', s=10
)
plt.colorbar(sc, label='Predictive Score')
plt.xlabel('X-scaled')
plt.ylabel('Y-scaled')
plt.title(f"Spatial Heatmap of Predictive Score (Top {top_n_genes} genes)")
plt.tight_layout()
plt.savefig("../figure/predictive_score.png", dpi=300)
plt.show()

# === Step 9: Plaque vs Cell Binary Heatmap ===
plt.figure(figsize=(8, 6))
# Cell points (label = 0, faint gray)
plt.scatter(score_df['X'], score_df['Y'], c='lightgray', s=8, label='Cells')

# Plaque points (label = 1, bright color)
plt.scatter(plaque_df['m.cx'], plaque_df['m.cy'], c='red', s=12, label='Plaques')

plt.xlabel('X')
plt.ylabel('Y')
plt.title("Spatial Heatmap of Plaque (Red) and Cells (Gray)")
plt.legend()
plt.tight_layout()
plt.savefig("../figure/plaque_vs_cells.png", dpi=300)
plt.show()

# === Step 7b: AUC score for spatial match with plaque ===

# Step 1: Build KDTree of plaque coordinates
plaque_tree = KDTree(plaque_df[['m.cx', 'm.cy']].values)

# Step 2: Define neighborhood radius (can tune this!)
radius = 30  # in coordinate units, adjust as needed

# Step 3: For each cell, check if it's within radius of any plaque
cell_coords = score_df[['X', 'Y']].values
is_near_plaque = [
    len(plaque_tree.query_ball_point(coord, r=radius)) > 0 for coord in cell_coords
]

# Step 4: Assign labels and compute AUC
score_df['near_plaque'] = is_near_plaque
auc = roc_auc_score(score_df['near_plaque'], score_df['Predictive_Score'])

print(f"✅ AUC (predictive score → plaque proximity): {auc:.4f}")
