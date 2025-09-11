import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# Input Data
# -----------------------------
data = {
    "Pathway": [
        "Inflammatory Response Pathway WP458",
        "Dysregulated miRNA Targeting In Insulin PI3K AKT Signaling WP3855",
        "Focal Adhesion WP85",
        "Focal Adhesion PI3K Akt mTOR Signaling Pathway WP2841",
        "Comprehensive IL 17A Signaling WP5242",
        "Regulation Of Actin Cytoskeleton WP523",
        "ApoE And miR 146 In Inflammation And Atherosclerosis WP3592",
        "Osteoblast Signaling WP238",
        "Burn Wound Healing WP5056",
        "Statin Pathway WP1"
    ],
    "P_value": [5.548e-9, 8.782e-7, 8.932e-6, 7.824e-5, 0.002640,
                0.005503, 0.005985, 0.007476, 0.01119, 0.01416],
    "Adjusted_pval": [7.767e-8, 6.147e-6, 4.168e-5, 0.0002738, 0.007392,
                      0.01197, 0.01197, 0.01308, 0.01741, 0.01982],
    "Odds_Ratio": [279.15, 216.98, 39.79, 22.49, 30.29, 20.62, 203.86, 158.54, 101.89, 79.23],
    "Combined_Score": [5306.53, 3025.86, 462.56, 212.65, 179.82, 107.28, 1043.43, 776.21, 457.73, 337.32]
}

df = pd.DataFrame(data)


# -----------------------------
df = df.head(7)

# -----------------------------
# Compute -log10(P-value)
# -----------------------------
df["-log10(P-value)"] = -np.log10(df["P_value"])

# -----------------------------
# Horizontal Bar Plot
# -----------------------------
fig, ax = plt.subplots(figsize=(12,7))

bars = ax.barh(df["Pathway"], df["-log10(P-value)"],
               color="#760185", edgecolor="black", linewidth=1.2)  # single solid color

# Manual left-aligned title
fig.text(0.10, 0.96, "-log10(P-value) of Enriched Pathways (7d post-MI)",
         fontsize=17, fontweight="bold", fontstyle="italic", fontfamily="serif",
         ha="left", va="center")

# Remove axis labels
ax.set_xlabel("")
ax.set_ylabel("")

# Bigger tick labels, bold + italic
for label in ax.get_yticklabels():
    label.set_fontweight("bold")
    label.set_fontstyle("italic")
    label.set_fontsize(17)
for label in ax.get_xticklabels():
    label.set_fontweight("bold")
    label.set_fontstyle("italic")
    label.set_fontsize(15)

# Gridlines on x-axis
ax.grid(axis="x", linestyle="--", alpha=0.5)

plt.tight_layout()
plt.savefig("Figure-4h.png", dpi=300)
plt.show()
