#!/bin/bash

# Fixed rank value
RANK=12

# List of lambda values
lambdas=(
    0.01 0.02 0.05 0.1 0.2 0.5
    1 2 5
    8 10 12
    14 15 16
    18 20 25
    30 40 50
    70 100 200
    500 1000
)


# lambdas=(
#     200
#     500 1000
# )


# Loop over each lambda and run the Python script
for lambda in "${lambdas[@]}"; do
    echo "Running with lambda = $lambda"
    python3 new_supervised_final_run.py --rank $RANK --lambda $lambda
done
