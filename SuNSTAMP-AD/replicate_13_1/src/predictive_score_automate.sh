#!/bin/bash

echo "GeneCount,AUC" > ../result/AUC_scores.csv  # Header
for i in $(seq 1 100); do
    echo "Running for top $i genes"
    python3 Predictive_scoring.py $i ../result/AUC_scores.csv
done
