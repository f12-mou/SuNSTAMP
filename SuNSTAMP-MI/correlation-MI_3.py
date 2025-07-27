import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import roc_auc_score

def analyze_gene_signal(gene_exp_path, supervision_path, gene_names_path, num_genes=8):
    # Load data
    gene_exp_df = pd.read_csv(gene_exp_path)
    supervision_df = pd.read_csv(supervision_path)
    gene_names_df = pd.read_csv(gene_names_path)

    # Select top N genes and their weights
    selected_genes_df = gene_names_df.head(num_genes)
    selected_genes = selected_genes_df['Gene'].tolist()
    weights = selected_genes_df['Value'].tolist()

    # Ensure gene columns exist in gene_exp_df
    missing_genes = [g for g in selected_genes if g not in gene_exp_df.columns]
    if missing_genes:
        raise ValueError(f"The following genes are missing in gene_exp.csv: {missing_genes}")

    # Find common cells
    common_cells = set(gene_exp_df['Cell_Id']).intersection(set(supervision_df['Cell_Id']))
    gene_exp_df = gene_exp_df[gene_exp_df['Cell_Id'].isin(common_cells)].copy()
    supervision_df = supervision_df[supervision_df['Cell_Id'].isin(common_cells)].copy()

    # Sort to ensure alignment
    gene_exp_df = gene_exp_df.sort_values(by='Cell_Id')
    supervision_df = supervision_df.sort_values(by='Cell_Id')

    # Compute weighted gene expression score for each cell
    gene_matrix = gene_exp_df[selected_genes].values
    weights_array = pd.Series(weights).values
    predictive_scores = gene_matrix.dot(weights_array)

    # Extract true target
    true_scores = supervision_df['Total_Attraction'].values

    # Calculate metrics
    pearson_corr, _ = pearsonr(predictive_scores, true_scores)
    spearman_corr, _ = spearmanr(predictive_scores, true_scores)
    
    try:
        auc = roc_auc_score(true_scores, predictive_scores)
    except ValueError:
        auc = None  # e.g., if true_scores are not binary or only one class exists

    return {
        'Pearson': pearson_corr,
        'Spearman': spearman_corr,
        'AUC-ROC': auc
    }

# Example usage:
if __name__ == "__main__":
    scores = analyze_gene_signal(
        gene_exp_path="gene_exp_3_3.csv",
        supervision_path="supervision_3_3.csv",
        gene_names_path="H_3_1.csv",
        num_genes=8
    )
    print("Correlation Results:")
    for k, v in scores.items():
        print(f"{k}: {v}")
