import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

# -----------------------------
# Load and preprocess
# -----------------------------
supervision_df = pd.read_csv("./data/supervision_13_1.csv")
w_df = pd.read_csv("./data/W_13_1.csv")

supervision_df = supervision_df[['Cell_Id', 'Normalized_Total_Attraction']]
w_df = w_df[['Cell_Id', 'Component_1']]

merged_df = pd.merge(supervision_df, w_df, on='Cell_Id', how='inner')

# -----------------------------
# Normalize the values
# -----------------------------
scaler = MinMaxScaler()
merged_df[['Component_1_norm', 'Normalized_Total_Attraction_norm']] = scaler.fit_transform(
    merged_df[['Component_1', 'Normalized_Total_Attraction']]
)

# -----------------------------
# Plot normalized scatter
# -----------------------------
plt.figure(figsize=(7,5))

plt.scatter(
    merged_df['Component_1_norm'],
    merged_df['Normalized_Total_Attraction_norm'],
    s=10,
    color='blue',
    alpha=0.7,
)

# -----------------------------
# Formatting
# -----------------------------
# Customize font sizes here
label_fontsize = 16   # x/y labels
legend_fontsize = 14  # legend

plt.xlabel("Normalized Supervision Signal",
           fontsize=label_fontsize,
           fontweight='bold',
           style='italic',
           fontfamily='serif')

plt.ylabel("Normalized Supervised Factor",
           fontsize=label_fontsize,
           fontweight='bold',
           style='italic',
           fontfamily='serif')

plt.grid(True, linestyle='--', alpha=0.5)

# Legend with bigger size
plt.legend(fontsize=legend_fontsize, frameon=False, loc="best")

# Remove excessive whitespace
plt.tight_layout()

plt.savefig("Figure-2d.png", dpi=300, bbox_inches='tight')
plt.show()
