#!/usr/bin/env python3
"""analyze_anchored_density.py
Calculates the density of splicing events (overlap coverage) anchored at Start and Stop codons.
Input: Per-event TSV (with transcript IDs and coordinates) and GTF (for CDS coordinates).
Output: TSV of density values suitable for plotting.
"""

import pandas as pd
import numpy as np
import argparse
import os
import gzip
from collections import defaultdict

def parse_gtf_cds(gtf_path):
    """
    Parses GTF to find Coding Start and Stop genomic coordinates for each transcript.
    Returns: dict {transcript_id: {'strand': +/-1, 'cds_min': int, 'cds_max': int}}
    """
    print(f"[INFO] Parsing GTF: {gtf_path}")
    tx_data = defaultdict(lambda: {'strand': None, 'cds_coords': []})
    
    open_func = gzip.open if gtf_path.endswith('.gz') else open
    
    with open_func(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            feature = parts[2]
            if feature != 'CDS': continue
            
            # Parse attributes
            attr_block = parts[8]
            tid = None
            for x in attr_block.split(';'):
                x = x.strip()
                if x.startswith('transcript_id'):
                    tid = x.split(' ')[1].replace('"', '')
                    break
            
            if not tid: continue
            
            strand = parts[6]
            start = int(parts[3])
            end = int(parts[4])
            
            # Strip version from tid
            tid_base = tid.split('.')[0]
            
            tx_data[tid_base]['strand'] = strand
            tx_data[tid_base]['cds_coords'].append(start)
            tx_data[tid_base]['cds_coords'].append(end)
            
    # Compile Start/Stop
    final = {}
    for tid, data in tx_data.items():
        if not data['cds_coords']: continue
        cds_min = min(data['cds_coords'])
        cds_max = max(data['cds_coords'])
        strand = data['strand']
        
        # Determine Start/Stop based on strand
        final[tid] = {
            'strand': strand,
            'cds_start_genomic': cds_min if strand == '+' else cds_max,
            'cds_stop_genomic': cds_max if strand == '+' else cds_min
        }
        
    print(f"[INFO] Found CDS info for {len(final)} transcripts.")
    return final

def parse_transcript_ids(val):
    if pd.isna(val): return []
    val = str(val)
    # Cleanup strict standard
    ids = []
    for x in val.replace(';', ',').split(','):
        x = x.strip()
        if not x: continue
        # Remove version if needed? No, usually keep strict matching or try both.
        ids.append(x)
    return ids

def calculate_density(events_df, tx_map, anchor_type, window=500):
    """
    anchor_type: 'Start_Codon' or 'Stop_Codon'
    window: +/- bp to cover
    events_df must have: RI coordinates (parsed from ID or cols)
    """
    # Grid: -window to +window
    # We'll use a dictionary to sum coverage: pos -> count
    density = defaultdict(int)
    n_events_used = 0
    
    for _, row in events_df.iterrows():
        # Parse Event Coords
        # Assuming rMATS format in ID: chr:start:end... or explicit columns
        # per_event_compact has 'event_id'
        # RMATS RI ID: chr:start:end:upstreamEE:downstreamES...
        # Wait, per_event_compact_for_plotting might NOT have coords columns explicitly?
        # But 'event_id' usually contains them.
        
        eid = str(row['event_id'])
        parts = eid.split('|') # RMATS ID format?
        # Standard rMATS ID (from rMATS output): chr:strand:start:end... 
        # But here checking previous output "chrX|+|154412085|154413567..."
        # It seems delimiter is '|'.
        # Format: chr|strand|ri_start|ri_end|upstreamEE|downstreamES
        # Let's verify from `prepare_rmats_data.R` or checking file.
        # User output showed: "chrX|+|154412085|154413567|..."
        # Index 0: chrX|+|... NO.
        # It's "chrX|+|154412085|154413567..." ?
        # Actually usually: chr|strand|ri_start|ri_end|...
        # Let's try splitting by '|'.
        
        # User showed: chrX|+|154412085|154413567|154412085|154412214|154413206|154413567
        # 0: chrX
        # 1: +
        # 2: 154412085 (starts)
        # 3: 154413567 (ends)
        # ...
        
        try:
            p = eid.split('|')
            chrom = p[0]
            strand = p[1]
            e_start = int(p[2])
            e_end = int(p[3])
        except:
            continue

        # Find Transcript
        tids = parse_transcript_ids(row.get('advanced_transcript_ids'))
        if not tids: continue
        
        # Pick best transcript (first one that's in our map)
        ref_tx = None
        for t in tids:
            # Try removing version
            base = t.split('.')[0]
            if t in tx_map: ref_tx = tx_map[t]; break
            if base in tx_map: ref_tx = tx_map[base]; break
            
        if not ref_tx: continue
        if ref_tx['strand'] != strand: continue # specific check
        
        # Get Anchor Genomic Coord
        anchor_g = 0
        if anchor_type == 'Start_Codon':
            anchor_g = ref_tx['cds_start_genomic']
        else: # Stop
            anchor_g = ref_tx['cds_stop_genomic']
            
        # Calculate Relative Coords
        # Strand +: Rel = Genomic - Anchor
        # Strand -: Rel = Anchor - Genomic
        
        if strand == '+':
            rel_start = e_start - anchor_g
            rel_end = e_end - anchor_g
        else:
            rel_start = anchor_g - e_end
            rel_end = anchor_g - e_start
            
        # Add to density coverage
        # Range [rel_start, rel_end) (0-based half open? rMATS coords?)
        # rMATS 0-based start, 1-based end?
        # Let's assume [start, end)
        
        # Optimize: Only iterate overlap with [-window, window]
        # Intersection
        s = max(rel_start, -window)
        e = min(rel_end, window)
        
        if s < e:
            n_events_used += 1
            for i in range(s, e):
                density[i] += 1
                
    return density, n_events_used

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_table", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--outfile", default=None, help="Custom output filename (default: anchored_density_plot_data.tsv inside outdir)")
    parser.add_argument("--window", type=int, default=500, help="Window size (+/- bp) around anchor point")
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # Load inputs
    df = pd.read_csv(args.input_table, sep='\t')
    tx_map = parse_gtf_cds(args.gtf)
    
    # We want density for UFM1_dependent vs Independent (Label/dataset column)
    # df has 'dataset' column? From `protein_primary_sequence_impact.py` output.
    if 'dataset' not in df.columns:
        print("[ERROR] Input table missing 'dataset' column.")
        return

    datasets = df['dataset'].unique()
    anchors = ['Start_Codon', 'Stop_Codon']
    
    all_data = []
    
    for ds in datasets:
        sub = df[df['dataset'] == ds]
        print(f"[INFO] Processing {ds} (n={len(sub)})...")
        
        for ans in anchors:
            den_map, n_used = calculate_density(sub, tx_map, ans, args.window)
            print(f"  > {ans}: Used associated transcripts for {n_used} events.")
            
            # Normalize? "Data used to generate such plot"
            # Usually users want 'Fraction of events' or just raw density + N.
            # Let's output density (count) and N (total used). User can normalize in R.
            # Or provide Normalized = Count / N_Used
            
            for pos in range(-args.window, args.window):
                cnt = den_map.get(pos, 0)
                norm = cnt / n_used if n_used > 0 else 0
                all_data.append({
                    "Dataset": ds,
                    "Anchor": ans,
                    "Position": pos,
                    "Count": cnt,
                    "Total_Events": n_used,
                    "Normalized_Density": norm
                })
                
    # Save
    out_df = pd.DataFrame(all_data)
    if args.outfile:
        out_file = args.outfile
    else:
        out_file = os.path.join(args.outdir, "anchored_density_plot_data.tsv")
        
    out_df.to_csv(out_file, sep='\t', index=False)
    print(f"[DONE] Saved density data to {out_file}")
    
    # Optional: Quick Python Plot to verify?
    # Keeping it simple as requested: "output the data... so i can manually plot later"
    # But a quick PNG helps verification.
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        # Smooth for Plotting (Sliding Window)
        # Python rolling mean
        out_df['Smoothed'] = out_df.groupby(['Dataset', 'Anchor'])['Normalized_Density'].transform(lambda x: x.rolling(window=50, center=True).mean())
        
        g = sns.FacetGrid(out_df, col="Anchor", hue="Dataset", height=5, aspect=1.5, sharey=False)
        g.map(sns.lineplot, "Position", "Smoothed")
        g.add_legend()
        plt.savefig(os.path.join(args.outdir, "anchored_density_RI.pdf"), format='pdf', bbox_inches='tight')
        print("[INFO] Saved preview plot.")
    except Exception as e:
        print(f"[WARN] Failed to plot preview: {e}")

if __name__ == "__main__":
    main()
