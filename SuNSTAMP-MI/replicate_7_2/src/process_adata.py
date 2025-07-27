import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

months = 7
replicate = 2
nuwaisir_vai_score_file = "../data/scores_D7_2_alpha_0.0001_lambda_0.1.csv"

# Load the AnnData object
adata = sc.read_h5ad(f"../data/adata_{months}_{replicate}.h5ad")

# Convert expression matrix to dense if sparse
expr_matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X

# Create DataFrame for expression data
expr_df = pd.DataFrame(expr_matrix, index=adata.obs.index, columns=adata.var.index)

# Insert 'Cell_Id' at top-left corner
expr_df.index.name = "Cell_Id"

# Save expression data
expr_df.to_csv(f"../data/gene_exp_{months}_{replicate}.csv")
print("Gene expression data saved to gene_exp.csv")

# Extract spatial coordinates
spatial_coords = adata.obsm["spatial"]

# Plotting
plt.figure(figsize=(6, 6))
plt.scatter(spatial_coords[:, 0], spatial_coords[:, 1], s=10, alpha=0.7)
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.title("Spatial Distribution of Spots")
plt.savefig(f"../figure/scatter_plot_{months}_{replicate}.png")
plt.show()

# Spatial DataFrame with Cell_Id
spatial_df = pd.DataFrame(spatial_coords, columns=["X", "Y"], index=adata.obs.index)
spatial_df.index.name = "Cell_Id"
spatial_df.to_csv(f"../data/spatial_{months}_{replicate}.csv")
print(f"Spot coordinates saved to spatial_{months}_{replicate}.csv")

# Load external data
external_df = pd.read_csv(nuwaisir_vai_score_file)

# Treat first unnamed column as Cell_Id
if external_df.columns[0] != "Cell_Id":
    external_df.rename(columns={external_df.columns[0]: "Cell_Id"}, inplace=True)

# Set Cell_Id as index for easy merging
external_df.set_index("Cell_Id", inplace=True)
spatial_df_for_merge = spatial_df.copy()

# Merge on Cell_Id (index)
merged_df = spatial_df_for_merge.join(external_df[["remodeling_score", "neg_log_pvalue"]], how="inner")

# Rename columns
merged_df = merged_df.rename(columns={
    "X": "X-scaled",
    "Y": "Y-scaled",
    "remodeling_score": "Total_Attraction",
    "neg_log_pvalue": "Normalized_Total_Attraction"
})

# Reorder columns
merged_df = merged_df[["X-scaled", "Y-scaled", "Total_Attraction", "Normalized_Total_Attraction"]]
merged_df.reset_index(inplace=True)  # Make Cell_Id a column again
merged_df = merged_df[["Cell_Id", "X-scaled", "Y-scaled", "Total_Attraction", "Normalized_Total_Attraction"]]

# Save supervision CSV
merged_df.to_csv(f"../data/supervision_{months}_{replicate}.csv", index=False)
print(f"Supervision file saved to supervision_{months}_{replicate}.csv")
