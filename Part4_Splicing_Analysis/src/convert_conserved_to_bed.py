import pandas as pd
import os

# Paths
base_dir = "mammalian_RI_dpsi01_fdr05"
summary_file = os.path.join(base_dir, "cross_species_visualization", "Conserved_Mammalian_Splicing_Summary.tsv")
human_file = os.path.join(base_dir, "human", "nucleus", "step01_data_prep", "UFM1_dependent.tsv")
mouse_file = os.path.join(base_dir, "mouse", "total", "step01_data_prep", "UFM1_dependent.tsv")

out_dir = os.path.join(base_dir, "cross_species_visualization", "IGV_BEDs")
os.makedirs(out_dir, exist_ok=True)

def create_bed(species_name, dependent_file, target_genes, output_path):
    print(f"Processing {species_name}...")
    try:
        df = pd.read_csv(dependent_file, sep='\t')
    except Exception as e:
        print(f"Error reading {dependent_file}: {e}")
        return

    # Filter for conserved genes
    # Check if 'geneSymbol' is in the target list
    # Targets might be mixed case, so normalise
    target_genes_upper = set([x.upper() for x in target_genes])
    
    # Filter
    # Ensure geneSymbol column exists
    if 'geneSymbol' not in df.columns:
        print(f"Column 'geneSymbol' not found in {dependent_file}")
        return
        
    df['geneSymbol_upper'] = df['geneSymbol'].astype(str).str.upper()
    subset = df[df['geneSymbol_upper'].isin(target_genes_upper)].copy()
    
    print(f"Found {len(subset)} events for {len(target_genes)} conserved genes.")
    
    # Extract Coordinates for Retained Intron
    # Intron starts after upstream exon end (upstreamEE)
    # Intron ends before downstream exon start (downstreamES)
    # BED is 0-based start, 1-based end (half-open)
    # rMATS upstreamEE is 1-based end of exon? Or 0-based?
    # rMATS is 0-based.
    # Start: upstreamEE (0-based end of exon IS 0-based start of intron)
    # End: downstreamES (0-based start of exon IS 0-based end of intron)
    
    # Columns usually: upstreamEE, downstreamES
    # Inmerged file they have suffix .WT usually?
    
    # Detect Coordinate Columns
    if 'upstreamEE.WT' in subset.columns:
        # Merged format
        subset['chrom'] = subset['chr.WT']
        subset['start'] = subset['upstreamEE.WT']
        subset['end'] = subset['downstreamES.WT']
        subset['strand'] = subset['strand.WT']
        subset['dPSI'] = subset['dPSI_num.WT'] # Optional score
        subset['gene'] = subset['geneSymbol']
    elif 'upstreamEE' in subset.columns:
        # Standard format
        subset['chrom'] = subset['chr']
        subset['start'] = subset['upstreamEE']
        subset['end'] = subset['downstreamES']
        subset['strand'] = subset['strand']
        subset['dPSI'] = 0
        subset['gene'] = subset['geneSymbol']
    else:
        print("Could not identify coordinate columns (upstreamEE/downstreamES)")
        return

    # Sanity Check
    # Valid intron: start < end
    subset = subset[subset['start'] < subset['end']]
    
    # Create BED fields
    # chrom, start, end, name, score, strand
    subset['score'] = (subset['dPSI'] * 1000).fillna(0).astype(int) # Fake score from dPSI
    
    # Add coordinates to name for uniqueness if multiple events per gene
    subset['name'] = subset['gene'] + ":" + subset['chrom'].astype(str) + ":" + \
                     subset['start'].astype(str) + "-" + subset['end'].astype(str)

    bed_df = subset[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    
    # Ensure chrom starts with chr if not? (IGV is flexible, but good practice)
    # subset['chrom'] = subset['chrom'].apply(lambda x: x if str(x).startswith('chr') else 'chr'+str(x))
    
    bed_df.to_csv(output_path, sep='\t', header=False, index=False)
    print(f"Saved {output_path}")

# 1. Load Summary to get gene list
print("Loading Summary...")
summary = pd.read_csv(summary_file, sep='\t')
conserved_genes = summary['geneSymbol_Upper'].unique().tolist()
print(f"Conserved Genes ({len(conserved_genes)}): {conserved_genes}")

# 2. Create Human BED
create_bed("Human", human_file, conserved_genes, os.path.join(out_dir, "Conserved_Targets_Human.bed"))

# 3. Create Mouse BED
create_bed("Mouse", mouse_file, conserved_genes, os.path.join(out_dir, "Conserved_Targets_Mouse.bed"))
