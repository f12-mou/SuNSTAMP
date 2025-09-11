# Pathway Enrichment Analysis (Top 20 Genes)

This folder contains the files and outputs related to the pathway enrichment analysis of the top 20 genes obtained from SuNSTAMP for 13 month old AD mouse model (replicate-1).

## Overview

- The input file `Top_Genes_13_1.csv` contains the ranked gene list.  
- The **top 20 genes** were selected and analyzed using the [Enrichr](https://maayanlab.cloud/Enrichr/) tool for pathway enrichment.  
- The resulting pathways were plotted using **-log10(P-value)** to visualize their significance.  

## Contents

- `Top_Genes_13_1.csv` – Ranked list of genes used as input.  
- `Figure-3g.py` – Python script for plotting the enrichment results (-log10(P)).  
- `Figure-3g.png` – Final figure showing the top pathways enriched for the top 20 genes.  

## How to Reproduce

1. Select the top 20 genes from `Top_Genes_13_1.csv`.  
2. Upload them to [Enrichr](https://maayanlab.cloud/Enrichr/).  
3. Download the enrichment results table from Enrichr.  
4. Use the python file to generate the figure.