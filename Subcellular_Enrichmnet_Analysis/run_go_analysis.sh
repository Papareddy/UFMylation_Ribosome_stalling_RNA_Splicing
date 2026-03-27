#!/bin/bash

# This script creates a pipeline to run Gene Ontology (GO) enrichment analysis.
# It first prepares the input data by selecting the 'gene_id' and 'cluster' columns
# from the source file. Then, it runs the simplifiedGO.py script and moves the
# output to the 'results' directory.

set -e

# Define file paths
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
INPUT_FILE="$SCRIPT_DIR/data/Microsome_ANS_vsDMSO_Col0_ufm1_Limma_with_log2counts.tsv"
TEMP_INPUT_FILE="$SCRIPT_DIR/data/temp_for_go.tsv"
PYTHON_SCRIPT="$SCRIPT_DIR/src/simplifiedGO.py"
RESULTS_DIR="$SCRIPT_DIR/results"

# Create a temporary input file with only gene_id and cluster columns
echo "Preparing input file for GO analysis..."
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {for(i=1;i<=NF;i++) {if($i=="gene_id") gid=i; if($i=="cluster") cid=i}; print "gene_id","cluster"; next} {print $gid, $cid}' "$INPUT_FILE" > "$TEMP_INPUT_FILE"

# Run the GO analysis script
echo "Running GO enrichment analysis..."
python "$PYTHON_SCRIPT" "$TEMP_INPUT_FILE"

# The python script creates a directory based on the input file name.
# Let's find it and move the contents.
TEMP_OUTPUT_DIR="temp_for_go_results" 

if [ -d "$TEMP_OUTPUT_DIR" ]; then
    echo "Moving results to the '$RESULTS_DIR' directory..."
    # Ensure the results directory exists
    mkdir -p "$RESULTS_DIR"
    # Move the contents of the temporary output directory to the results directory
    # The -n option prevents overwriting existing files.
    mv -n "$TEMP_OUTPUT_DIR"/* "$RESULTS_DIR/"
    # Remove the now-empty temporary output directory
    rmdir "$TEMP_OUTPUT_DIR"
    echo "Results have been moved."
else
    echo "Warning: Output directory '$TEMP_OUTPUT_DIR' not found. Nothing to move."
fi

# Clean up the temporary input file
echo "Cleaning up temporary files..."
rm "$TEMP_INPUT_FILE"

echo "Pipeline finished."
