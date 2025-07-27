import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load results
df = pd.read_csv('../result/AUC_scores.csv')

# Subset to top 50 genes
df = df[df['GeneCount'] <= 50]

# Plot
plt.figure(figsize=(8, 6))
sns.lineplot(data=df, x='GeneCount', y='AUC', marker='o', linewidth=2.5, color='navy')

# Customize
plt.title("AUC vs Number of Top Genes Used", fontsize=16, fontweight='bold')
plt.xlabel("Number of Top Genes", fontsize=14, fontweight='bold')
plt.ylabel("AUC Score", fontsize=14, fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.4)
plt.axvline(16, color='blue', linestyle='--', label='Top 16 genes')
plt.axhline(0.97795, color='green', linestyle='--', label='Maximum AUC score')
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("../figure/auc_curve_top50.png", dpi=300)
plt.show()
