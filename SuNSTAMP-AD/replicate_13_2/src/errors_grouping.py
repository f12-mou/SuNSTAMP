import pandas as pd
import os

# Initialize an empty list to store DataFrames
all_dfs = []

# Define the folder path and other constants
months = 13
replicate = 2
folder_path = f'../result/mou_errors_{months}_{replicate}'

# Iterate over i from 0 to 20 (i represents rank 5 to 25)
for i in range(21):
    rank = i + 5
    filename = f'mou_errors_{i}_{months}_{replicate}.csv'
    filepath = os.path.join(folder_path, filename)
    
    try:
        # Read file without header and assign column names
        df = pd.read_csv(filepath, header=None, names=['file_rank', 'lambda', 'loss'])

        # Replace file's rank column with correct rank
        df['rank'] = rank

        # Keep only relevant columns
        df = df[['rank', 'lambda', 'loss']]

        all_dfs.append(df)

    except Exception as e:
        print(f"Error reading {filepath}: {e}")

# Concatenate all DataFrames
combined_df = pd.concat(all_dfs, ignore_index=True)

# Group by rank and lambda, then sum the loss
grouped_df = combined_df.groupby(['rank', 'lambda'], as_index=False)['loss'].sum()

# Save the final DataFrame to CSV
output_file = f'../result/summed_loss_by_rank_lambda_{months}_{replicate}.csv'
grouped_df.to_csv(output_file, index=False)

print(f"Summed loss values have been saved to '{output_file}'")
