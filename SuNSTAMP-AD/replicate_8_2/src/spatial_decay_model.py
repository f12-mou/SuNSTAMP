import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Constants
PI = np.pi
months = 8
replicate = 2
scale_factor = 300

# Load cell data
cell_df = pd.read_csv(f'../data/spatial_{months}_{replicate}.csv', skiprows=[1]) 
# Load plaque data
plaque_df = pd.read_csv(f'../data/plaque_{months}_{replicate}.csv')

# Function to calculate squared root distance between two points
def euclid_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2)


# Function to calculate ideal impact value
def impact_value(radius, euclid_dist):
    res = (4/3.0) * PI * (radius) * (radius) * (radius) 
    ret_val = res * math.exp(-euclid_dist/scale_factor)
    return ret_val

# List to store the results
results = []
start=1
# Calculate the sum of attraction values for each cell
for _, cell in cell_df.iterrows():
    total_attraction = 0
    for _, plaque in plaque_df.iterrows():
        # print(cell['X-scaled'])
        # print(cell['Y-scaled'])
        # print(plaque['m.cx'])
        # print(plaque['m.cy'])
        # print(plaque['s.radius.mean'])
        euclid_dist = euclid_distance(cell['X-scaled'], cell['Y-scaled'], plaque['m.cx'], plaque['m.cy'])
        attraction = impact_value(plaque['s.radius.mean'], euclid_dist)
        total_attraction += attraction
    print (start)
    start=start+1
    results.append({'Cell_Id': cell['NAME'], 'X-scaled': cell['X-scaled'], 'Y-scaled': cell['Y-scaled'], 'Total_Attraction': total_attraction})

# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Min-Max normalization
min_val = results_df['Total_Attraction'].min()
max_val = results_df['Total_Attraction'].max()
results_df['Normalized_Total_Attraction'] = ((results_df['Total_Attraction'] - min_val) / (max_val - min_val) if max_val > min_val else 0)


# Save the DataFrame to a CSV file
results_df.to_csv(f'../data/supervision_{months}_{replicate}.csv', index=False)


# Filtering DataFrame for plotting
#filtered_df = results_df[results_df['Normalized_Total_Attraction'] >= 0.2]
filtered_df = results_df
# Plotting the normalized attraction values
plt.figure(figsize=(12, 10))  # Set the size of the figure
plt.scatter(filtered_df['X-scaled'], filtered_df['Y-scaled'], c=filtered_df['Total_Attraction'], cmap='YlOrRd', s=5)

# Plotting the normalized attraction values
# plt.figure(figsize=(12, 10))  # Set the size of the figure
# plt.scatter(results_df['X-scaled'], results_df['Y-scaled'], c=results_df['Normalized_Total_Attraction'], cmap='YlOrRd', s=5)
plt.colorbar(label='Normalized Attraction Value')
plt.title('Cells Colormap for Attraction Values')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.savefig(f'../figure/normalized_attraction_values_{months}_{replicate}.png')  # Save the figure
plt.show()  # Display the figure



# === Plotting Cells and Plaques with Radius ===
plt.figure(figsize=(12, 10))

# Plot cells
plt.scatter(results_df['X-scaled'], results_df['Y-scaled'],
            c='blue', s=5, label='Cells', alpha=0.6)

# Add plaques as red circles with radius scaled
ax = plt.gca()  # Get current axis
for _, plaque in plaque_df.iterrows():
    circle = patches.Circle(
        (plaque['m.cx'], plaque['m.cy']), 
        radius=plaque['s.radius.mean'],
        color='red',
        fill=False,
        linewidth=1.2,
        alpha=0.8
    )
    ax.add_patch(circle)

plt.title('Cells and Plaques with Radius Overlay')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.legend(['Cells', 'Plaques (radius shown)'])
plt.grid(True)
plt.axis('equal')  # Equal scaling for x and y
plt.savefig(f'../figure/cells_and_plaques_with_radii_{months}_{replicate}.png')  # Save figure
plt.show()

