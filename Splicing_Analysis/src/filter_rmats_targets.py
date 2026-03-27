import pandas as pd
import os

# Inputs
bed_file = "mammalian_RI_dpsi01_fdr05/cross_species_visualization/IGV_BEDs/Conserved_Targets_Human.bed"
rmats_file = "data/human/nucleus/WT_ctrl_vs_ANS_NUCL/RI.MATS.JCEC.txt"
out_filtered = "miscellaneous/sashimi_plots_output/Conserved_RI.MATS.JCEC.txt"

# Ensure output directory
os.makedirs("miscellaneous/sashimi_plots_output", exist_ok=True)

# 1. Get Target Gene List
print(f"Reading targets from {bed_file}...")
# BED format: chrom, start, end, name, score, strand
# Name is usually Gene:Chr:Start-End. We extract Gene.
try:
    bed_df = pd.read_csv(bed_file, sep='\t', header=None)
    # Extract gene name. Assuming name format "Gene:..."
    bed_df['gene'] = bed_df[3].apply(lambda x: x.split(':')[0])
    targets = set(bed_df['gene'].unique())
    # Add uppercase just in case
    targets_upper = set([x.upper() for x in targets])
    print(f"Found {len(targets)} unique target genes.")
except Exception as e:
    print(f"Error reading BED: {e}")
    exit(1)

# 2. Filter rMATS file
print(f"Filtering {rmats_file}...")
try:
    # rMATS is tab separated
    rmats_df = pd.read_csv(rmats_file, sep='\t')
    
    # Clean gene symbol (remove formatting if any)
    # The head output showed "TAFAZZIN" with quotes maybe? 
    # Pandas read_csv handles standard CSV quotes, but checked 'head' output: "TAFAZZIN" (with quotes).
    # Pandas usually strips quotes if generic, but let's strict check.
    
    # Normalize geneSymbol column
    if 'geneSymbol' not in rmats_df.columns:
        print("Column 'geneSymbol' not found!")
        print(rmats_df.columns)
        exit(1)
        
    rmats_df['gene_upper'] = rmats_df['geneSymbol'].astype(str).str.upper().str.strip('"')
    
    # Filter
    filtered = rmats_df[rmats_df['gene_upper'].isin(targets_upper)]
    
    print(f"Filtered down to {len(filtered)} events matching conserved targets.")
    
    # Save
    filtered.drop(columns=['gene_upper']).to_csv(out_filtered, sep='\t', index=False, quoting=3) # QUOTE_NONE? rMATS uses quotes?
    # The original file had quotes around strings.
    # To match input format exactly for rmats2sashimiplot, we might need quotes?
    # rmats2sashimiplot parses standard rMATS.
    # Let's save standard TSV.
    # Actually, rmats2sashimiplot is Python based, likely pandas based reading. Standard TSV is fine.
    
    # Check if empty
    if len(filtered) == 0:
        print("WARNING: No events found matching targets!")
    else:
        print(f"Saved filtered events to {out_filtered}")
        
except Exception as e:
    print(f"Error processing rMATS file: {e}")
    exit(1)
