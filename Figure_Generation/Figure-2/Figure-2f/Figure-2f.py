import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from sklearn.preprocessing import MinMaxScaler
# -----------------------------
# Load data
# -----------------------------
w_df = pd.read_csv('./data/W_13_1.csv')  # contains Cell_Id, Component_1, etc.
coordinates_df = pd.read_csv('./data/spatial_13_1.csv')

# Merge spatial coordinates with Component_1
merged_df = pd.merge(w_df[['Cell_Id', 'Component_1']], coordinates_df, on='Cell_Id')
merged_df.rename(columns={'X-scaled_x': 'X-scaled', 'Y-scaled_x': 'Y-scaled'}, inplace=True)

# -----------------------------
# Normalize Component_1
# -----------------------------
scaler = MinMaxScaler()
merged_df['Component_1_norm'] = scaler.fit_transform(merged_df[['Component_1']])
# -----------------------------
# Number of bins
# -----------------------------
num_bins = 5
merged_df["Bin"] = pd.qcut(
    merged_df["Component_1_norm"], 
    q=num_bins, 
    labels=False
)

# -----------------------------
# Define custom start & end colors (RGB)
# -----------------------------
start_rgb = (0, 255, 0)   # bright green
end_rgb   = (24, 61, 27)  # dark green

# Normalize RGB values to [0,1] for matplotlib
start_rgb = tuple([v/255 for v in start_rgb])
end_rgb   = tuple([v/255 for v in end_rgb])

# Create custom colormap
cmap = LinearSegmentedColormap.from_list("custom_cmap", [start_rgb, end_rgb], N=num_bins)

# Generate discrete colors for bins
colors = [cmap(i) for i in range(num_bins)]

# Map bins to colors
merged_df["Color"] = merged_df["Bin"].map(lambda x: colors[int(x)])

# Sort so lower Component_1 plotted first
merged_df = merged_df.sort_values(by="Component_1_norm", ascending=True)

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(7, 6))
plt.scatter(merged_df["X-scaled"], merged_df["Y-scaled"], 
            c=merged_df["Color"], s=12, alpha=1)

# -----------------------------
# Add colorbar with styled label
# -----------------------------
sm = plt.cm.ScalarMappable(
    cmap=cmap,
    norm=plt.Normalize(
        vmin=merged_df["Component_1_norm"].min(),
        vmax=merged_df["Component_1_norm"].max()
    )
)
cbar = plt.colorbar(sm)
cbar.set_label("Normalized Supervised Factor", 
               fontsize=14,
               fontweight='bold',
               style='italic',
               fontfamily='serif')

# -----------------------------
# Clean plot (no ticks or spines)
# -----------------------------
plt.xticks([])
plt.yticks([])
for spine in plt.gca().spines.values():
    spine.set_visible(False)



# Save high-res
plt.savefig("Figure-2f.png", bbox_inches='tight', dpi=600)
plt.show()
