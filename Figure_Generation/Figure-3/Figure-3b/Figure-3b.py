import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the two gene expression matrices
# Assuming the first column is cell IDs and the first row contains gene names
file_path_1 = "gene_exp_8_1.csv"  # Replace with your file path
file_path_2 = "gene_exp_13_1.csv"  # Replace with your file path

data_1 = pd.read_csv(file_path_1, index_col=0)
data_2 = pd.read_csv(file_path_2, index_col=0)

print(data_1.shape)
print(data_2.shape)

# Specify the gene to compare
gene_name = "Apoe"  # Replace with your gene of interest
if gene_name not in data_1.columns:
    print("not in 8 months")
if gene_name not in data_2.columns:
    print("not in 13 motnhs")

if gene_name not in data_1.columns or gene_name not in data_2.columns:
    print(f"Gene {gene_name} not found in one or both datasets!")
else:
    # Extract gene expression values
    gene_values_1 = data_1[gene_name]
    gene_values_2 = data_2[gene_name]
    print(len(gene_values_1), len(gene_values_2))
    

    # Create a DataFrame for Seaborn
    combined_data = pd.DataFrame({
        "Expression": pd.concat([gene_values_1, gene_values_2]),
        "Dataset": ["8 months"] * len(gene_values_1) + ["13 months"] * len(gene_values_2)
    })

    # Plot the distribution
    plt.figure(figsize=(6, 4))
    sns.kdeplot(data=combined_data, x="Expression", hue="Dataset", fill=True,alpha=0.5, common_norm=False, bw_adjust=0.5)
    plt.title(f"Comparison of {gene_name} Expression in Two Datasets")
    plt.xlabel("Expression Level")
    plt.ylabel("Density")
    plt.savefig("Apoe.jpg")
    plt.show()
