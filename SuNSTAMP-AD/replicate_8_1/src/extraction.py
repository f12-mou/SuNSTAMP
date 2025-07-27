import pandas as pd


# The input file name is constant here
file_name_1 = "/kaggle/input/ad-dataset/expression_matrix_normalized.csv"
start_index_spot = 59045
end_index_spot = 68206
# Output_File
output_file_1 = "gene_exp_8_1.csv"





# Loading the CSV file
df = pd.read_csv(file_name_1, index_col=0, dtype="str")  # The first column contains gene names and should be set as index
#print(df)
# Filter columns based on the range of IDs interested in
# Convert column names to integers for comparison, if necessary
filtered_columns = [col for col in df.columns if col.isdigit() and start_index_spot <= int(col) <= end_index_spot]
filtered_df = df[filtered_columns]

# Transpose the DataFrame
transposed_df = filtered_df.T

# Save the transposed DataFrame to a new CSV file
transposed_df.to_csv(output_file_1)

print("Extraction and transposition completed.")
