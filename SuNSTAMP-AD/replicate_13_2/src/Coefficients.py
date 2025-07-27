import pandas as pd
import numpy as np
import os
import re
from scipy.stats import spearmanr, pearsonr

months = 13
replicate = 2
rank = 14
# also change below in stripping and 25 number line

# Paths
w_folder = f'../result/W_matrices_{months}_{replicate}'
supervision_file = f'../data/supervision_{months}_{replicate}.csv'
output_file = f'../result/Supervised_Factor_correlation_results_{months}_{replicate}.csv'

# Load the supervision data
supervision_df = pd.read_csv(supervision_file)

# Initialize list to collect results
results = []

# Loop over all W matrix files
for fname in os.listdir(w_folder):
    if fname.startswith(f'output_W_matrix_rank_{rank}_') and fname.endswith('_13_2.csv'):
        file_path = os.path.join(w_folder, fname)

        # Extract lambda from filename using regex
        match = re.search(r'rank_(\d+)_([\deE\.-]+)_13_2\.csv', fname)
        if not match:
            print(f"Filename doesn't match pattern: {fname}")
            continue
        rank = int(match.group(1))
        lambda_val = float(match.group(2))

        # Read W matrix
        w_df = pd.read_csv(file_path)

        # Merge on Cell_Id
        merged_df = pd.merge(w_df, supervision_df, on='Cell_Id')

        # Replace inf/-inf with 0
        merged_df.replace([np.inf, -np.inf], 0, inplace=True)

        # Check if Component_1 and Total_Attraction exist and are valid
        if 'Component_1' in merged_df.columns and 'Total_Attraction' in merged_df.columns:
            # Calculate correlations
            try:
                pcc, _ = pearsonr(merged_df['Component_1'], merged_df['Total_Attraction'])
                scc, _ = spearmanr(merged_df['Component_1'], merged_df['Total_Attraction'])
            except Exception as e:
                print(f"Correlation failed for {fname}: {e}")
                pcc, scc = np.nan, np.nan
        else:
            print(f"Required columns missing in {fname}")
            pcc, scc = np.nan, np.nan

        # Append result
        results.append({
            'rank': rank,
            'lambda': lambda_val,
            'pearson_corr': pcc,
            'spearman_corr': scc
        })

# Convert to DataFrame and save
results_df = pd.DataFrame(results)
results_df.sort_values(['rank', 'lambda'], inplace=True)
results_df.to_csv(output_file, index=False)

print(f"Correlation results saved to: {output_file}")
