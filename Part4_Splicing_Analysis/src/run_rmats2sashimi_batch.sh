#!/bin/bash

# ==============================================================================
# BATCH RMATS2SASHIMIPLOT GENERATOR (4-GROUP MODE) - FIXED .GF FORMAT
# ==============================================================================
# Purpose: Generate sashimi plots for conserved targets with all 4 groups.
# Groups: WT_CTRL, WT_ANS, UFM1_CTRL, UFM1_ANS
# ==============================================================================

# --- 1. CONFIGURATION ---
BAM_DIR="/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/merged_bams"
RMATS_FILTERED="miscellaneous/sashimi_plots_output/Conserved_RI.MATS.JCEC.txt"
OUTPUT_DIR="miscellaneous/sashimi_plots_output/rmats_plots_4group"
TOOLS_DIR="miscellaneous/tools"

# Colors (User Requested)
# Order: WT_CTRL, WT_ANS, UFM1_CTRL, UFM1_ANS
C_WT_CTRL="#2196F3"    # Blue
C_WT_ANS="#F44336"     # Red
C_UFM_CTRL="#4CAF50"   # Green
C_UFM_ANS="#FF9800"    # Orange

mkdir -p "$OUTPUT_DIR"

echo "[INFO] Starting 4-Group rMATS2Sashimi Batch Process..."

# --- 2. ENVIRONMENT SETUP ---
TOOL_DIR="$TOOLS_DIR/rmats2sashimiplot"
EXEC_SCRIPT=$(find "$TOOL_DIR" -name "rmats2sashimiplot.py" | head -n 1)

# --- 3. BAM GROUPING AND .GF GENERATION ---

# Define the specific BAM files we want in order
# Ordering here is CRITICAL for indices 1, 2, 3, 4
BAM_WT_CTRL="$BAM_DIR/wt_ctrl_NUCL_merged.bam"
BAM_WT_ANS="$BAM_DIR/wt_ans_NUCL_merged.bam"
BAM_UFM_CTRL="$BAM_DIR/ufm1_ctrl_dox_NUCL_merged.bam"
BAM_UFM_ANS="$BAM_DIR/ufm1_ans_dox_NUCL_merged.bam"

# Check existence
for b in "$BAM_WT_CTRL" "$BAM_WT_ANS" "$BAM_UFM_CTRL" "$BAM_UFM_ANS"; do
    if [ ! -f "$b" ]; then echo "[ERROR] Missing BAM: $b"; exit 1; fi
done

# Prepare --b1 and --b2
# rmats2sashimiplot merges b1 and b2 internally. 
# Total list index will be: 1:b1[0], 2:b1[1], 3:b2[0], 4:b2[1]
B1_LIST="$BAM_WT_CTRL,$BAM_WT_ANS"
B2_LIST="$BAM_UFM_CTRL,$BAM_UFM_ANS"

# Create Group Info File (.gf)
# Format required by MISO: GroupName: Index (1-based)
GF_FILE="$OUTPUT_DIR/groups.gf"
cat <<EOF > "$GF_FILE"
WT_CTRL:1
WT_ANS:2
UFM1_CTRL:3
UFM1_ANS:4
EOF

# Color string must match the order of groups in .gf
# Since we defined WT_CTRL first, WT_ANS second, etc.
COLOR_STR="$C_WT_CTRL,$C_WT_ANS,$C_UFM_CTRL,$C_UFM_ANS"

echo "[INFO] Using .gf format:"
cat "$GF_FILE"

# --- 4. EXECUTION ---
echo "[INFO] Running rmats2sashimiplot with 4 groups..."

python3 "$EXEC_SCRIPT" \
    --b1 "$B1_LIST" \
    --b2 "$B2_LIST" \
    --event-type RI \
    -e "$RMATS_FILTERED" \
    --l1 "WT" \
    --l2 "UFM1" \
    --group-info "$GF_FILE" \
    --color "$COLOR_STR" \
    --exon_s 1 \
    --intron_s 1 \
    -o "$OUTPUT_DIR"

echo "[INFO] Done. Results in $OUTPUT_DIR"
