import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the main CSV file
file_path = './data/spatial_13_1.csv'  
data = pd.read_csv(file_path)

# Load the additional CSV file for plaques
plaque_file_path = './data/plaque_13_1.csv'  
plaque_data = pd.read_csv(plaque_file_path)

# Set font styles
plt.rcParams['font.family'] = 'serif'  # Options: 'sans-serif', 'serif', 'monospace', etc.
plt.rcParams['font.size'] = 14         # Base font size
plt.rcParams['font.weight'] = 'bold'   # Font weight
plt.rcParams['font.style'] = 'italic'  # Font style: 'normal', 'italic', 'oblique'


cell_types = data['top_level_cell_type'].unique()
palette = sns.color_palette('hsv', len(cell_types))
color_map = {cell_type: palette[i] for i, cell_type in enumerate(cell_types)}

# Plot
plt.figure(figsize=(10, 8))

# Plot cell types
for cell_type, color in color_map.items():
    subset = data[data['top_level_cell_type'] == cell_type]
    plt.scatter(subset['X-scaled'], subset['Y-scaled'], label=cell_type, s=10, color=color)

# Plot plaques as black dots
for _, row in plaque_data.iterrows():
    plt.scatter(
        row['m.cx'],  # X coordinate
        row['m.cy'],  # Y coordinate
        s=row['s.radius.mean'] * 1.8,  # Scale the radius for better visibility
        color='black',
        alpha=0.9,
        label='Plaques' if 'Plaques' not in plt.gca().get_legend_handles_labels()[1] else None  # Avoid duplicate legend entries
    )


handles, labels = plt.gca().get_legend_handles_labels()
labels.append('Plaques')

# Add legend
plt.legend(
    handles=handles,
    labels=labels,
    title='Annotations',
    bbox_to_anchor=(1.03, 0.8),
    loc='upper left',
    scatterpoints=1,
    markerscale=2.5,  # Uniform marker scale
    fontsize=12,
    title_fontsize=14
)

# Remove axes and grid
plt.gca().set_xticks([])
plt.gca().set_yticks([])
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)

# Add title
plt.title('Scatter Plot of Cell Types and Plaques', fontsize=16, weight='bold', style='italic')

# Adjust layout
plt.tight_layout()

# Save and show the plot
plt.savefig("Figure-2b.png", dpi=600)
plt.show()


# Count the number of cells for each top-level cell type
cell_type_counts = data['top_level_cell_type'].value_counts()

# Print the counts for all cell types
print("Cell counts for each top-level cell type:\n", cell_type_counts)

# Print the number of cells in the top (most frequent) cell type
top_cell_type = cell_type_counts.idxmax()
top_cell_count = cell_type_counts.max()
print(f"\nTop cell type: '{top_cell_type}' with {top_cell_count} cells.")
