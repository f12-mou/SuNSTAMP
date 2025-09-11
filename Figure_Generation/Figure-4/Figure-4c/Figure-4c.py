import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv("H_3_1.csv")   # csv: Gene, Value
df_top10 = df.head(10)

# -----------------------------
# Plot
# -----------------------------
fig, ax = plt.subplots(figsize=(11,7))

bars = ax.bar(df_top10['Gene'], df_top10['Value'],
              color="#38A0EA", edgecolor="black", linewidth=1.2)  # solid color

# -----------------------------
# Styling
# -----------------------------
ax.set_title("Top 10 Genes from H Matrix (3 days post-MI)",
             fontsize=24, fontfamily="serif",
             fontstyle="italic", fontweight="bold", pad=10)

# Remove x and y axis labels (Gene / Value)
ax.set_xlabel("")
ax.set_ylabel("")

# Bigger tick values
ax.tick_params(axis="x", labelrotation=30, labelsize=24)
ax.tick_params(axis="y", labelsize=24)

# Apply font styling to ticks
for tick in ax.get_xticklabels() + ax.get_yticklabels():
    tick.set_fontfamily("serif")
    tick.set_fontstyle("italic")
    tick.set_fontweight("bold")

# Set y-axis limit (upper bound = 0.30)
ax.set_ylim(0, 0.0098)

# Light grid
ax.grid(axis="y", linestyle="--", alpha=0.5)

plt.tight_layout()
plt.savefig("Figure-4c.png", dpi=300)
plt.show()
