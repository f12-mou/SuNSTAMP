import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# -----------------------------
# Input Data
# -----------------------------
data = {
    "Pathway": [
        "Microglia Pathogen Phagocytosis Pathway",
        "ApoE And miR 146 In Inflammation And Atherosclerosis",
        "Iron Homeostasis",
        "Complement Activation Classical Pathway",
        "Statin Pathway",
        "White Fat Cell Differentiation",
        "P38 Mapk Signaling Pathway",
        "Oxidative Damage Response",
        "Dravet Syndrome Scn1a Mutation Model",
        "Alzheimer 39S Disease"
    ],
    "Pval": [0.0007610, 0.007973, 0.01391, 0.01687, 0.01884,
             0.03153, 0.03250, 0.04023, 0.06962, 0.07149],
    "Adjusted_pval": [0.01218, 0.06028, 0.06028, 0.06028, 0.06028,
                      0.07429, 0.07429, 0.08046, 0.1144, 0.1144],
    "Combined_Score": [407.96, 725.58, 345.57, 268.08, 231.83,
                       117.08, 112.42, 84.30, 39.33, 37.87],
    "Overlap": [10, 8, 7, 6, 6, 5, 5, 4, 3, 3]
}

df = pd.DataFrame(data)
df["neglogP"] = -np.log10(df["Pval"])   # <-- compute -log10(p-value)

# -----------------------------
# Horizontal Barplot (-log10 Pval)
# -----------------------------
plt.figure(figsize=(10,6))
sns.barplot(
    x="neglogP", y="Pathway", data=df,
    color="purple", edgecolor="black"   # <-- single color
)

plt.xlabel("")
plt.ylabel("")

# Style text
ax = plt.gca()



# Custom title placement
ax.set_title(
    "-log10(P-value) of Enriched Pathways",
    fontsize=16, fontweight="bold",
    loc="left",     # "left", "center", or "right"
    x=0.0,          # 0.0 = far left, 0.5 = center, 1.0 = far right
    y=1.05          # raise/lower title position
)

for label in ax.get_yticklabels():
    label.set_fontweight("bold")
    label.set_fontstyle("italic")
    label.set_fontsize(18)
for label in ax.get_xticklabels():
    label.set_fontweight("bold")
    label.set_fontstyle("italic")
    label.set_fontsize(14)

plt.tight_layout()
plt.savefig("Figure-3g.png", dpi=600)
plt.show()
