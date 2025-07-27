#!/bin/bash

# Fixed rank value
RANK=14

# List of lambda values
lambdas=(
    3e-6 2.75e-6 2.5e-6 2e-6 1.75e-6 1.5e-6 1.25e-6 1e-6
    3e-5 2.75e-5 2.5e-5 2e-5 1.75e-5 1.5e-5 1.25e-5 1e-5
    3e-4 2.75e-4 2.5e-4 2e-4 1.75e-4 1.5e-4 1.25e-4 1e-4
    1e-3 1e-2
)

# Loop over each lambda and run the Python script
for lambda in "${lambdas[@]}"; do
    echo "Running with lambda = $lambda"
    python3 new_supervised_final_run.py --rank $RANK --lambda $lambda
done
