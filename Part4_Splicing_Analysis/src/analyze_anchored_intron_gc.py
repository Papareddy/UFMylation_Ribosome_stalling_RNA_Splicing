#!/usr/bin/env python3
"""
analyze_anchored_intron_gc.py

Calculates the GC content of UFM1-dependent, UFM1-independent, and Constitutive introns
anchored at Start and Stop codons.

Logic Update (Gene-based Anchoring):
1. Load GTF: Identify the "canonical" (longest CDS) transcript for each Gene. Record its Start/Stop codons.
2. Load Events: Get Intron coordinates and associated Gene ID from TSV.
3. Anchoring: For each intron, find its Gene's canonical transcript. Use that transcript's Start/Stop for anchoring.
"""

import argparse
import sys
import os
import pandas as pd
import numpy as np
import pysam
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import random

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze Anchored Intron GC Content")
    parser.add_argument("--dependent", required=True, help="TSV of UFM1-dependent events")
    parser.add_argument("--independent", required=True, help="TSV of UFM1-independent events")
    parser.add_argument("--gtf", required=True, help="GTF file")
    parser.add_argument("--genome", required=True, help="Genome FASTA")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--window", type=int, default=1000, help="Window size (bp)")
    return parser.parse_args()

def load_gtf_data(gtf_path):
    """
    Parses GTF to identify the longest CDS transcript for each gene.
    Returns:
      gene_map: {gene_id: {'chrom': c, 'strand': s, 'cds_start': start, 'cds_stop': stop}}
      all_introns: List of constitutive intron candidates.
    """
    print(f"[INFO] Parsing GTF: {gtf_path}")
    
    # Store transcript CDS info temporarily
    tx_data = defaultdict(lambda: {'start_codons': [], 'stop_codons': [], 'exons': [], 'gene_id': None, 'strand': None, 'chrom': None})
    
    if gtf_path.endswith('.gz'):
        import gzip
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
        
    with open_func(gtf_path, mode) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            chrom = parts[0]
            feature = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attr = parts[8]
            
            if feature not in ['exon', 'start_codon', 'stop_codon', 'CDS']: continue
            
            # Extract Transcript ID and Gene ID
            m_tid = re.search(r'transcript_id "([^"]+)"', attr)
            if not m_tid: continue
            tid = m_tid.group(1).split('.')[0]
            
            m_gid = re.search(r'gene_id "([^"]+)"', attr)
            gid = m_gid.group(1).split('.')[0] if m_gid else None
            
            t = tx_data[tid]
            t['gene_id'] = gid
            t['chrom'] = chrom
            t['strand'] = strand
            
            if feature == 'start_codon': t['start_codons'].append((start, end))
            elif feature == 'stop_codon': t['stop_codons'].append((start, end))
            elif feature == 'exon': t['exons'].append((start, end))
            elif feature == 'CDS': 
                t.setdefault('cds_len', 0)
                t['cds_len'] += (end - start + 1)

    print(f"[INFO] Loaded {len(tx_data)} transcripts.")
    
    # Select Best Transcript per Gene
    gene_map = {}
    all_introns = []
    
    for tid, data in tx_data.items():
        gid = data['gene_id']
        if not gid: continue
        
        # Calculate CDS length/info
        if 'cds_len' not in data: data['cds_len'] = 0
        
        # Determine Start/Stop
        if not data['start_codons'] or not data['stop_codons']:
            continue
            
        if data['strand'] == '+':
            cds_start = min([x[0] for x in data['start_codons']])
            cds_stop = max([x[1] for x in data['stop_codons']])
        else:
            cds_start = max([x[1] for x in data['start_codons']])
            cds_stop = min([x[0] for x in data['stop_codons']])
            
        # Update Gene Map (Keep longest CDS)
        if gid not in gene_map or data['cds_len'] > gene_map[gid]['cds_len']:
            gene_map[gid] = {
                'chrom': data['chrom'],
                'strand': data['strand'],
                'cds_start': cds_start,
                'cds_stop': cds_stop,
                'cds_len': data['cds_len'],
                'tid': tid
            }
            
        # Collect Introns for Constitutive Pool
        exons = sorted(data['exons'])
        if len(exons) > 1:
            for i in range(len(exons) - 1):
                all_introns.append({
                    'chrom': data['chrom'],
                    'start': exons[i][1] + 1,
                    'end': exons[i+1][0] - 1,
                    'strand': data['strand'],
                    'gene_id': gid
                })

    print(f"[INFO] Identified best transcripts for {len(gene_map)} genes.")
    return gene_map, all_introns

def get_event_coords(tsv_path):
    """
    Parses TSV to get events (chrom, start, end, strand, gene_id).
    """
    print(f"[INFO] Loading Events from {tsv_path}")
    df = pd.read_csv(tsv_path, sep='\t')
    events = []
    
    # Helper to find col
    def find_col(aliases):
        for a in aliases:
            if a in df.columns: return a
        return None

    # Mapping
    c_col = find_col(['chr', 'seqid', 'chrom', 'chr.WT', 'chr.UFM'])
    s_col = find_col(['strand', 'STRAND', 'strand.WT', 'strand.UFM'])
    uee_col = find_col(['upstreamEE', 'upstreamEE.WT', 'upstreamEE.UFM'])  
    des_col = find_col(['downstreamES', 'downstreamES.WT', 'downstreamES.UFM'])
    
    # Try multiple GeneID columns
    gene_col = find_col(['GeneID', 'gene_id', 'ENSG', 'GeneID.WT']) 
    if not gene_col:
        gene_col = find_col(['geneSymbol', 'geneSymbol.WT']) # Fallback (will try to map to gene_id later? No, map must match GTF)

    # If parsing failed or critical columns missing
    if not (c_col and s_col and uee_col and des_col and gene_col):
        print(f"[WARN] Missing columns in {tsv_path}. Cols: {list(df.columns)}")
        return events

    for _, row in df.iterrows():
        try:
            chrom = str(row[c_col])
            strand = row[s_col]
            # rMATS coordinates -> GTF 1-based Intron
            istart = int(row[uee_col]) + 1
            iend = int(row[des_col])   
            
            # Extract GeneID
            gid_raw = str(row[gene_col])
            gid = gid_raw.split('.')[0] # Remove version if present
            
            # Normalize chrom prefix
            if not chrom.startswith('chr'): chrom = 'chr' + chrom
            
            events.append({
                'chrom': chrom, 'start': istart, 'end': iend, 'strand': strand, 'gene_id': gid
            })
        except Exception as e:
            continue
            
    print(f"[INFO] Loaded {len(events)} events from {tsv_path}.")
    return events

def calculate_gc_percentages(intron_list, gene_map, genome, window):
    """
    Calculates GC % profile anchored at Start/Stop of the associated gene's best transcript.
    """
    
    # Allocate arrays
    start_gc_sum = np.zeros(2*window + 1, dtype=np.float64)
    start_depth = np.zeros(2*window + 1, dtype=np.float64)
    stop_gc_sum = np.zeros(2*window + 1, dtype=np.float64)
    stop_depth = np.zeros(2*window + 1, dtype=np.float64)
    
    processed_introns = 0
    fasta = pysam.FastaFile(genome)
    
    for intron in intron_list:
        gid = intron.get('gene_id')
        if not gid or gid not in gene_map:
             # Try without version stripping? (already stripped)
             continue
             
        tx_info = gene_map[gid]
        cds_start = tx_info['cds_start']
        cds_stop = tx_info['cds_stop']
        strand = tx_info['strand']
        chrom = tx_info['chrom'] # Use chrom from Gene Map (GTF source)
        
        # We don't need intron start/end for Gene-Anchored Analysis (centered on CDS Start/Stop)
        # We just need to fetch sequence around the CDS Start/Stop.
        
        # Fetch logic handles chrom prefix mismatch internally if needed

        
        # Fetch Window around CDS Start
        # Genomic coords:
        try:
            # Handle Map to Start
            if strand == '+':
                center = cds_start
                # Region: [center - window, center + window]
                # Indices: 0-based fetch.
                reg_start = center - window - 1
                reg_end = center + window
                # Fetch
                reg_seq = fasta.fetch(chrom, reg_start, reg_end).upper() if chrom in fasta.references else fasta.fetch(chrom.replace('chr',''), reg_start, reg_end).upper()
            else:
                center = cds_start # Max coord for negative strand start
                # Region: [center - window, center + window]?
                # For Neg strand: Start is the HIGHER coord.
                # Upstream is HIGHER. Downstream is LOWER.
                # Window: [center - window, center + window]
                # Then Reverse Complement.
                reg_start = center - window - 1
                reg_end = center + window
                reg_seq = fasta.fetch(chrom, reg_start, reg_end).upper()
                reg_seq = reg_seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                
            if len(reg_seq) == 2*window + 1:
                gc = np.fromiter((1 if b in 'GC' else 0 for b in reg_seq), dtype=np.float64)
                start_gc_sum += gc
                start_depth += 1
                
            # Handle Map to Stop
            if strand == '+':
                center = cds_stop
                reg_start = center - window - 1
                reg_end = center + window
                reg_seq = fasta.fetch(chrom, reg_start, reg_end).upper()
            else:
                center = cds_stop # Min coord for negative strand stop
                reg_start = center - window - 1
                reg_end = center + window
                reg_seq = fasta.fetch(chrom, reg_start, reg_end).upper()
                reg_seq = reg_seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

            if len(reg_seq) == 2*window + 1:
                gc = np.fromiter((1 if b in 'GC' else 0 for b in reg_seq), dtype=np.float64)
                stop_gc_sum += gc
                stop_depth += 1
                
            processed_introns += 1

        except Exception as e:
            continue

    print(f"Processed {processed_introns} events/genes for GC profile.")
    with np.errstate(divide='ignore', invalid='ignore'):
        start_profile = start_gc_sum / start_depth
        stop_profile = stop_gc_sum / stop_depth
    return start_profile, stop_profile, start_depth, stop_depth

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    gene_map, all_introns = load_gtf_data(args.gtf)
    print(f"[INFO] Gene Map Size: {len(gene_map)}")
    
    dep_events = get_event_coords(args.dependent)
    indep_events = get_event_coords(args.independent)
    
    # Constitutive: Random sample of introns -> map to genes
    # Or just random sample of GENES that are not dependent/independent?
    # Better: Use 'all_introns' to get GeneIDs.
    
    # Identify Dependent and Independent GENES
    dep_genes = set(e['gene_id'] for e in dep_events if e['gene_id'])
    indep_genes = set(e['gene_id'] for e in indep_events if e['gene_id'])
    
    # Filter dep/indep events to unique genes to avoid over-counting if a gene has multiple events?
    # Usually "per event" is fine, but if we answer "GC content around Start Codon", we should do "per Gene".
    # Let's do PER GENE to avoid bias from multi-intron genes?
    # Actually, RIs are events. If a gene has 2 RIs, it might be relevant twice?
    # But for Start Codon analysis, it's the SAME start codon.
    # So we should use UNIQUE GENES.
    
    # Create lists of Gene IDs for each group
    group_dep = list(dep_genes)
    group_indep = list(indep_genes)
    
    # Const Genes = All Genes in gene_map - Dep - Indep
    all_genes = set(gene_map.keys())
    const_genes = list(all_genes - dep_genes - indep_genes)
    
    # Subsample Const
    if len(const_genes) > 20000:
        print(f"[INFO] Subsampling constitutive genes from {len(const_genes)} to 20000")
        random.seed(42)
        const_genes = random.sample(const_genes, 20000)
        
    print(f"[INFO] Gene Counts: Dependent={len(group_dep)}, Independent={len(group_indep)}, Constitutive={len(const_genes)}")
    
    # To reuse calculate_gc_percentages, pass dummy events with gene_id
    # We just need [{'gene_id': g, 'chrom': ...}]
    # But we need chrom for the fetch call in loop (though gene_map has it).
    # Refactor loop to use gene_map's chrom.
    
    def make_list(genes):
        l = []
        for g in genes:
            if g in gene_map:
                l.append({'gene_id': g, 'chrom': gene_map[g]['chrom']})
        return l
        
    l_dep = make_list(group_dep)
    l_indep = make_list(group_indep)
    l_const = make_list(const_genes)
    
    print("[INFO] Calculating GC Profile for Dependent...")
    d_start, d_stop, _, _ = calculate_gc_percentages(l_dep, gene_map, args.genome, args.window)
    
    print("[INFO] Calculating GC Profile for Independent...")
    i_start, i_stop, _, _ = calculate_gc_percentages(l_indep, gene_map, args.genome, args.window)
    
    print("[INFO] Calculating GC Profile for Constitutive...")
    c_start, c_stop, _, _ = calculate_gc_percentages(l_const, gene_map, args.genome, args.window)
    
    # Save Data
    x_axis = np.arange(-args.window, args.window + 1)
    
    data = []
    for i, x in enumerate(x_axis):
        data.append({'Position': x, 'Group': 'Lost', 'Feature': 'Start_Codon', 'GC_Percent': d_start[i]})
        data.append({'Position': x, 'Group': 'Lost', 'Feature': 'Stop_Codon', 'GC_Percent': d_stop[i]})
        data.append({'Position': x, 'Group': 'Preserved', 'Feature': 'Start_Codon', 'GC_Percent': i_start[i]})
        data.append({'Position': x, 'Group': 'Preserved', 'Feature': 'Stop_Codon', 'GC_Percent': i_stop[i]})
        data.append({'Position': x, 'Group': 'Constitutive', 'Feature': 'Start_Codon', 'GC_Percent': c_start[i]})
        data.append({'Position': x, 'Group': 'Constitutive', 'Feature': 'Stop_Codon', 'GC_Percent': c_stop[i]})
        
    df = pd.DataFrame(data)
    outfile = os.path.join(args.outdir, "anchored_intron_gc_data.tsv")
    df.to_csv(outfile, sep='\t', index=False)
    print(f"[DONE] Saved data to {outfile}")
    
    # Plot
    print("[INFO] Plotting...")
    sns.set_theme(style="whitegrid")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    sns.lineplot(data=df[df['Feature']=='Start_Codon'], x='Position', y='GC_Percent', hue='Group', ax=axes[0])
    axes[0].set_title(f"Start Codon Anchored GC Content")
    axes[0].set_ylabel("GC Fraction")
    
    sns.lineplot(data=df[df['Feature']=='Stop_Codon'], x='Position', y='GC_Percent', hue='Group', ax=axes[1])
    axes[1].set_title(f"Stop Codon Anchored GC Content")
    axes[1].set_ylabel("GC Fraction")
    
    plt.tight_layout()
    plotfile = os.path.join(args.outdir, "anchored_intron_gc_plot.pdf")
    plt.savefig(plotfile, format='pdf', bbox_inches='tight')
    print(f"[DONE] Saved plot to {plotfile}")

if __name__ == "__main__":
    main()
