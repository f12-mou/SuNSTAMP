import pandas as pd
from scipy.stats import mannwhitneyu


months = 13
replicate = 2

# Step 1: Read the grouped loss CSV
df = pd.read_csv(f'../result/summed_loss_by_rank_lambda_{months}_{replicate}.csv')

# Step 2: Initialize dictionary for losses per rank
loss_by_rank = {rank: [] for rank in range(5, 26)}  # ranks 5 to 25 inclusive

# Step 3: Populate the dictionary
for _, row in df.iterrows():
    rank = row['rank']
    loss = row['loss']
    if rank in loss_by_rank:
        loss_by_rank[rank].append(loss)

# Optional: quick check
for rank in range(5, 7):
    print(f"Rank {rank}: {loss_by_rank[rank]}")

# Step 4: Perform Mann-Whitney U Test between consecutive ranks
results = []
for rank in range(5, 25):  # compare rank vs rank+1, so up to 24
    losses_rank1 = loss_by_rank[rank]
    losses_rank2 = loss_by_rank[rank + 1]

    # Only compare if both lists have at least 2 values
    if len(losses_rank1) > 1 and len(losses_rank2) > 1:
        stat, p_value = mannwhitneyu(losses_rank1, losses_rank2, alternative='two-sided')
        results.append({
            'Rank1': rank,
            'Rank2': rank + 1,
            'U_Statistic': stat,
            'P_Value': p_value
        })

# Step 5: Save to CSV
results_df = pd.DataFrame(results)
results_df.to_csv(f'../result/mann_whitney_results_summed_{months}_{replicate}.csv', index=False)

print("Results have been saved to 'mann_whitney_results_summed.csv'")
