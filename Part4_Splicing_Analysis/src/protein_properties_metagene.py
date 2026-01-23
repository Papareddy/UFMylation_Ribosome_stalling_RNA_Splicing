#!/usr/bin/env python3
"""
Protein Property Metagene Analysis
Calculates Bio.SeqUtils.ProtParam properties along protein length for:
- Control proteins
- UFM1-dependent RI proteins  
- UFM1-independent RI proteins
"""

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

def parse_gtf(gtf_file):
    """Extract CDS regions per transcript"""
    import gzip
    from collections import defaultdict
    cds_dict = defaultdict(list)
    
    # Handle gzipped files
    open_func = gzip.open if gtf_file.endswith('.gz') else open
    mode = 'rt' if gtf_file.endswith('.gz') else 'r'
    
    with open_func(gtf_file, mode) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS': continue
            
            chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            attrs = dict(item.strip().split(' ', 1) for item in parts[8].rstrip(';').split(';') if ' ' in item.strip())
            
            tx_id = attrs.get('transcript_id', '').strip('"')
            gene_id = attrs.get('gene_id', '').strip('"').split('.')[0]  # Strip version
            
            if tx_id and gene_id:
                cds_dict[tx_id].append({
                    'chr': chrom, 'start': start, 'end': end, 
                    'strand': strand, 'gene_id': gene_id
                })
    
    return cds_dict

def extract_cds_sequence(cds_list, genome_fasta):
    """Extract and translate CDS to protein"""
    if not cds_list: return None
    
    # Sort by position
    strand = cds_list[0]['strand']
    cds_sorted = sorted(cds_list, key=lambda x: x['start'], reverse=(strand == '-'))
    
    # Extract sequences
    seq_parts = []
    for cds in cds_sorted:
        try:
            # genome_fasta dict returns SeqRecord, need .seq
            seq_record = genome_fasta[cds['chr']]
            seq = str(seq_record.seq[cds['start']-1:cds['end']])
            seq_parts.append(seq)
        except:
            return None
    
    cds_seq = ''.join(seq_parts)
    if strand == '-':
        cds_seq = str(Seq(cds_seq).reverse_complement())
    
    # Translate
    if len(cds_seq) % 3 != 0 or len(cds_seq) < 30:
        return None
        
    try:
        protein = str(Seq(cds_seq).translate(to_stop=True))
        return protein if len(protein) >= 10 else None  # Min 10 AA
    except:
        return None

def calculate_properties_window(seq, window_size=30, step=10):
    """Calculate ProtParam properties in sliding windows"""
    results = []
    
    for i in range(0, len(seq) - window_size + 1, step):
        window = seq[i:i+window_size]
        
        # Skip if too many X or non-standard AA
        if window.count('X') > 3 or any(aa not in 'ACDEFGHIKLMNPQRSTVWY' for aa in window):
            continue
        
        try:
            analyzer = ProteinAnalysis(window)
            
            results.append({
                'position': i + window_size/2,  # Center of window
                'molecular_weight': analyzer.molecular_weight(),
                'aromaticity': analyzer.aromaticity(),
                'instability': analyzer.instability_index(),
                'isoelectric_point': analyzer.isoelectric_point(),
                'gravy': analyzer.gravy(),  # Hydropathy
                'helix_fraction': analyzer.secondary_structure_fraction()[0],
                'turn_fraction': analyzer.secondary_structure_fraction()[1],
                'sheet_fraction': analyzer.secondary_structure_fraction()[2],
                'charge_at_pH7': sum([
                    window.count('K') + window.count('R') -  # Positive
                    window.count('D') - window.count('E')     # Negative
                ]) / len(window),
                'aliphatic_index': analyzer.aliphatic_index(),
            })
        except:
            continue
    
    return pd.DataFrame(results)

def normalize_to_metagene(df_list, nbins=100):
    """Normalize protein positions to 0-100% for metagene"""
    metagene_data = defaultdict(list)
    
    for df in df_list:
        if df.empty: continue
        
        # Normalize positions to 0-1
        df['rel_pos'] = df['position'] / df['position'].max()
        
        # Bin into 100 bins
        df['bin'] = pd.cut(df['rel_pos'], bins=nbins, labels=False)
        
        # Average per bin
        binned = df.groupby('bin').mean()
        
        for col in binned.columns:
            if col not in ['position', 'rel_pos', 'bin']:
                metagene_data[col].append(binned[col].values)
    
    # Average across all proteins
    result = {}
    for prop, values in metagene_data.items():
        stacked = np.array([v for v in values if len(v) == nbins])
        if len(stacked) > 0:
            result[prop] = np.nanmean(stacked, axis=0)
            result[f'{prop}_sem'] = np.nanstd(stacked, axis=0) / np.sqrt(len(stacked))
    
    return result

def main():
    parser = argparse.ArgumentParser(description='Protein property metagene analysis')
    parser.add_argument('--gtf', required=True, help='GTF annotation file')
    parser.add_argument('--genome', required=True, help='Genome FASTA file')
    parser.add_argument('--dep_genes', required=True, help='UFM1-dependent gene list (one per line)')
    parser.add_argument('--indep_genes', required=True, help='UFM1-independent gene list (one per line)')
    parser.add_argument('--outdir', default='.', help='Output directory')
    parser.add_argument('--window', type=int, default=30, help='Sliding window size (AA)')
    parser.add_argument('--step', type=int, default=10, help='Step size (AA)')
    parser.add_argument('--nbins', type=int, default=100, help='Number of bins for metagene')
    args = parser.parse_args()
    
    print("=" * 70)
    print("PROTEIN PROPERTY METAGENE ANALYSIS")
    print("=" * 70)
    
    # Load gene lists
    print("\n[1/6] Loading gene lists...")
    with open(args.dep_genes) as f:
        dep_genes = set(line.strip().split('.')[0] for line in f)
    with open(args.indep_genes) as f:
        indep_genes = set(line.strip().split('.')[0] for line in f)
    
    print(f"  Dependent genes: {len(dep_genes)}")
    print(f"  Independent genes: {len(indep_genes)}")
    
    # Parse GTF
    print("\n[2/6] Parsing GTF...")
    cds_dict = parse_gtf(args.gtf)
    print(f"  Found {len(cds_dict)} transcripts with CDS")
    
    # Group transcripts by gene
    gene_to_tx = defaultdict(list)
    for tx_id, cds_list in cds_dict.items():
        if cds_list:
            gene_id = cds_list[0]['gene_id']
            gene_to_tx[gene_id].append(tx_id)
    
    # Load genome
    print("\n[3/6] Loading genome...")
    genome = SeqIO.to_dict(SeqIO.parse(args.genome, 'fasta'))
    print(f"  Loaded {len(genome)} chromosomes")
    
    # Extract proteins
    print("\n[4/6] Extracting proteins (longest CDS per gene)...")
    
    def get_longest_protein(gene_id, group_name):
        """Get longest protein isoform for a gene"""
        if gene_id not in gene_to_tx:
            return None
        
        best_protein = None
        max_len = 0
        
        for tx_id in gene_to_tx[gene_id]:
            protein = extract_cds_sequence(cds_dict[tx_id], genome)
            if protein and len(protein) > max_len:
                best_protein = protein
                max_len = len(protein)
        
        return best_protein
    
    proteins_dep = {}
    proteins_indep = {}
    proteins_control = {}
    
    # Dependent
    for i, gene in enumerate(dep_genes, 1):
        if i % 50 == 0:
            print(f"  Processing dependent: {i}/{len(dep_genes)}")
        prot = get_longest_protein(gene, 'Dependent')
        if prot:
            proteins_dep[gene] = prot
    
    # Independent
    for i, gene in enumerate(indep_genes, 1):
        if i % 50 == 0:
            print(f"  Processing independent: {i}/{len(indep_genes)}")
        prot = get_longest_protein(gene, 'Independent')
        if prot:
            proteins_indep[gene] = prot
    
    # Control (all other genes)
    all_genes = set(gene_to_tx.keys())
    control_genes = all_genes - dep_genes - indep_genes
    
    # Sample control to match dependent count (avoid bias)
    import random
    random.seed(42)
    control_sample = random.sample(list(control_genes), min(len(proteins_dep) * 3, len(control_genes)))
    
    for i, gene in enumerate(control_sample, 1):
        if i % 100 == 0:
            print(f"  Processing control: {i}/{len(control_sample)}")
        prot = get_longest_protein(gene, 'Control')
        if prot:
            proteins_control[gene] = prot
    
    print(f"\n  Extracted proteins:")
    print(f"    Dependent: {len(proteins_dep)}")
    print(f"    Independent: {len(proteins_indep)}")
    print(f"    Control: {len(proteins_control)}")
    
    # Calculate properties
    print(f"\n[5/6] Calculating properties (window={args.window}, step={args.step})...")
    
    prop_dep = [calculate_properties_window(p, args.window, args.step) for p in proteins_dep.values()]
    prop_indep = [calculate_properties_window(p, args.window, args.step) for p in proteins_indep.values()]
    prop_control = [calculate_properties_window(p, args.window, args.step) for p in proteins_control.values()]
    
    # Create metagenes
    print("\n[6/6] Creating metagene profiles...")
    meta_dep = normalize_to_metagene(prop_dep, args.nbins)
    meta_indep = normalize_to_metagene(prop_indep, args.nbins)
    meta_control = normalize_to_metagene(prop_control, args.nbins)
    
    # Plot
    print("\nGenerating plots...")
    properties = ['gravy', 'aromaticity', 'instability', 'isoelectric_point', 
                  'charge_at_pH7', 'helix_fraction', 'sheet_fraction', 'aliphatic_index']
    
    fig, axes = plt.subplots(4, 2, figsize=(14, 16))
    axes = axes.flatten()
    
    x = np.linspace(0, 100, args.nbins)
    colors = {'Control': '#9E9E9E', 'Dependent': '#E57373', 'Independent': '#64B5F6'}
    
    for idx, prop in enumerate(properties):
        ax = axes[idx]
        
        # Plot each group
        for name, meta, color in [('Control', meta_control, colors['Control']),
                                    ('Dependent', meta_dep, colors['Dependent']),
                                    ('Independent', meta_indep, colors['Independent'])]:
            if prop in meta:
                y = meta[prop]
                sem = meta.get(f'{prop}_sem', np.zeros_like(y))
                
                ax.plot(x, y, label=name, color=color, linewidth=2)
                ax.fill_between(x, y-sem, y+sem, color=color, alpha=0.2)
        
        ax.set_xlabel('Protein Position (%)')
        ax.set_ylabel(prop.replace('_', ' ').title())
        ax.legend()
        ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{args.outdir}/protein_properties_metagene.pdf')
    plt.savefig(f'{args.outdir}/protein_properties_metagene.png', dpi=300)
    print(f"\n✓ Saved plots to {args.outdir}/")
    
    # Save data
    for name, meta in [('Control', meta_control), ('Dependent', meta_dep), ('Independent', meta_indep)]:
        df = pd.DataFrame({'Position_pct': x})
        for prop in properties:
            if prop in meta:
                df[prop] = meta[prop]
                df[f'{prop}_SEM'] = meta.get(f'{prop}_sem', np.nan)
        
        df.to_csv(f'{args.outdir}/metagene_{name}.tsv', sep='\t', index=False)
    
    print(f"✓ Saved data tables to {args.outdir}/")
    print("\nDone!")

if __name__ == '__main__':
    main()
