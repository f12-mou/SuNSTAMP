import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load both results and label source
df1 = pd.read_csv('./replicate_13_1/result/AUC_scores_13_1.csv')
df1['Replicate'] = '13_1'

df2 = pd.read_csv('./replicate_13_2/result/AUC_scores_13_2.csv')
df2['Replicate'] = '13_2'

# Combine both
df = pd.concat([df1, df2], ignore_index=True)

# Filter top 50 genes
df = df[df['GeneCount'] <= 50]

# Plot
plt.figure(figsize=(9, 6))
sns.lineplot(data=df, x='GeneCount', y='AUC', hue='Replicate', style='Replicate',
             markers=True, dashes=False, linewidth=2.5, palette='Set1')

# Customize
plt.title("AUC vs Number of Top Genes Used (Replicate 13_1 vs 13_2)", fontsize=16, fontweight='bold')
plt.xlabel("Number of Top Genes", fontsize=14, fontweight='bold')
plt.ylabel("AUC Score", fontsize=14, fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.4)

# Horizontal line at max AUC
max_auc = df['AUC'].max()
plt.axhline(max_auc, color='green', linestyle='--', label=f'Max AUC = {max_auc:.3f}')

# Vertical line at top N genes (if any significant point)
plt.axvline(16, color='gray', linestyle='--', alpha=0.7, label='Top 16 genes')

plt.legend(fontsize=12, title='Replicate')
plt.tight_layout()
plt.savefig("auc_curve_top50_both.png", dpi=300)
plt.show()
