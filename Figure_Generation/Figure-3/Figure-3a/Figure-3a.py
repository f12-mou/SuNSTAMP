import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv("H_genes.csv")   # csv: Gene, Value
df_top10 = df.head(10)

# -----------------------------
# Plot
# -----------------------------
fig, ax = plt.subplots(figsize=(10,6))

bars = ax.bar(df_top10['Gene'], df_top10['Value'],
              color="#38A0EA", edgecolor="black", linewidth=1.2)  # solid color

# Add value labels on top of bars (bigger fontsize now)
for bar in bars:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, height + 0.02*max(df_top10['Value']),
            f"{height:.2f}", ha="center", va="bottom",
            fontsize=16, fontfamily="serif", fontstyle="italic", fontweight="bold")

# -----------------------------
# Styling
# -----------------------------
ax.set_title("Values of the Top 10 Genes from the First Row of H Matrix",
             fontsize=20, fontfamily="serif",
             fontstyle="italic", fontweight="bold", pad=20)

# Remove x and y axis labels (Gene / Value)
ax.set_xlabel("")
ax.set_ylabel("")

# Bigger tick values
ax.tick_params(axis="x", labelrotation=30, labelsize=20)
ax.tick_params(axis="y", labelsize=20)

# Apply font styling to ticks
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontfamily("serif")
    tick.set_fontstyle("italic")
    tick.set_fontweight("bold")

# Set y-axis limit (upper bound = 0.30)
ax.set_ylim(0, 0.30)

# Light grid
ax.grid(axis="y", linestyle="--", alpha=0.5)

plt.tight_layout()
plt.savefig("Figure-3a.png", dpi=300)
plt.show()
