import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# Input Data
# -----------------------------
data = {
    "Pathway": [
        "Iron Homeostasis WP1596",
        "Regulation Of Actin Cytoskeleton WP523",
        "Translation Factors WP307",
        "Cytoplasmic Ribosomal Proteins WP163",
        "Myometrial Relaxation And Contraction Pathways WP385",
        "EGFR1 Signaling Pathway WP572",
        "Focal Adhesion WP85"
    ],
    "P_value": [0.00004753, 0.005503, 0.03541, 0.05908, 0.1061, 0.1235, 0.1302],
    "Adjusted_pval": [0.0003327, 0.01926, 0.08263, 0.1034, 0.1302, 0.1302, 0.1302],
    "Odds_Ratio": [256.06, 20.62, 30.30, 17.77, 9.57, 8.13, 7.69],
    "Combined_Score": [2548.92, 107.28, 101.23, 50.28, 21.47, 17.01, 15.67]
}

df = pd.DataFrame(data)

# -----------------------------
# Compute -log10(P-value)
# -----------------------------
df["-log10(P-value)"] = -np.log10(df["P_value"])

# -----------------------------
# Horizontal Bar Plot
# -----------------------------
fig, ax = plt.subplots(figsize=(10,6))

bars = ax.barh(df["Pathway"], df["-log10(P-value)"],
               color="#760185", edgecolor="black", linewidth=1.2)  # single solid color

# Title
# ax.set_title("-log10(P-value) of Enriched Pathways", 
#              fontsize=18, fontweight="bold", fontstyle="italic", fontfamily="serif", loc="left")


# ax.text(-0.5, 1.05, "-log10(P-value) of Enriched Pathways", 
#         transform=ax.transAxes, fontsize=20, fontweight="bold",
#         fontstyle="italic", fontfamily="serif", ha="left", va="center")

fig.text(0.10, 0.95, "-log10(P-value) of Enriched Pathways (3d post-MI)",
         fontsize=18, fontweight="bold", fontstyle="italic", fontfamily="serif",
         ha="left", va="center")
# Remove axis labels
ax.set_xlabel("")
ax.set_ylabel("")

# Bigger tick labels, bold + italic
for label in ax.get_yticklabels():
    label.set_fontweight("bold")
    label.set_fontstyle("italic")
    label.set_fontsize(18)
for label in ax.get_xticklabels():
    label.set_fontweight("bold")
    label.set_fontstyle("italic")
    label.set_fontsize(16)

# Gridlines on x-axis
ax.grid(axis="x", linestyle="--", alpha=0.5)

plt.tight_layout()
plt.savefig("Figure-4g.png", dpi=300)
plt.show()
