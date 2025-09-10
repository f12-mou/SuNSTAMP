import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Pearson values for supervised and unsupervised approaches
# Each row corresponds to a rank (1 to 10), columns represent supervised and unsupervised
pearson_values = [
    [0.919532, -0.189289],
    [-0.129293, -0.225143],
    [-0.227380, 0.293001],
    [-0.027361, 0.090439],
    [-0.182549, -0.052712],
    [-0.137277, 0.022208],
    [-0.092424, -0.156556],
    [0.020621, 0.088067],
    [-0.134442, -0.139305],
    [-0.251070, -0.080750]
]

# Convert to a numpy array
pearson_values = np.array(pearson_values)

# Set font style
plt.rcParams.update({'font.family': 'serif', 'font.size': 12, 'font.weight': 'bold'})

# Create a heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(
    pearson_values,
    annot=True,  # Show the Pearson values in the cells
    fmt=".2f",  # Format values to 2 decimal places
    cmap="YlOrRd",  # Yellow to Red colormap
    xticklabels=["Supervised", "Unsupervised"],
    yticklabels=[f"Factor {i}" for i in range(1, 11)],
    cbar_kws={"label": "Pearson Correlation"}
)

# Add labels and title
#plt.xlabel("Approach")
#plt.ylabel("Rank")
#plt.title("Heatmap of Pearson Correlation for Supervised and Unsupervised Approaches")

# Show the plot
plt.tight_layout()
plt.savefig("comparison_page_1_5.png")
plt.show()
