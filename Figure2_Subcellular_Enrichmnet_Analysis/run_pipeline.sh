#!/bin/bash

# Get the directory where the script is located
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo "Running All-in-One GO Analysis and Plotting Pipeline..."

# The python script now handles everything:
# - Finding clusters
# - Running analysis for each
# - Combining results
# - Generating the final plot

python "${SCRIPT_DIR}/src/simplifiedGO.py" \
    --input_file "${SCRIPT_DIR}/data/Microsome_ANS_vsDMSO_Col0_ufm1_Limma_with_log2counts.tsv" \
    --output_dir "${SCRIPT_DIR}/results"

echo "Pipeline finished. Check the '${SCRIPT_DIR}/results' directory."
