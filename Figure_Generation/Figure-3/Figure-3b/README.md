# Gene Expression Comparison (8 vs 13 Months old AD)

## Overview
This script compares the expression distribution of a selected gene between two datasets (8 months vs 13 months). It generates a Kernel density plot using Seaborn.

## Required Files
- `gene_exp_8_1.csv` – Gene expression matrix for the 8-month dataset.  
- `gene_exp_13_1.csv` – Gene expression matrix for the 13-month dataset.  


## Output
- `Apoe.jpg` – Example output showing the distribution of **Apoe** expression across the two datasets.  

## Usage
Update the `gene_name` variable in the script to any gene of interest, then run:

```bash
python compare_expression.py
