#!/usr/bin/env python3
"""
ese_motif_scanning.py
Step B: ESE, RBP, and CAG Motif Scanning.

Performs:
1. CAG trinucleotide counting (PRIMARY METRIC)
2. FIMO-based RBP motif scanning using Human CIS-BP database
"""

import pandas as pd
import argparse
import os
import subprocess
import sys
from collections import defaultdict
from scipy import stats
import re


def count_cag_in_sequences(fasta_path, exon_window=(-30, -1)):
    """
    Count CAG motifs in the exonic portion of ROI_1 sequences.
    
    For 5'SS ±30bp sequences (61bp total):
    - Exonic portion: positions 0-29 (first 30 bases, representing -30 to -1)
    - Intronic portion: positions 30-60 (last 31 bases, representing 0 to +30)
    
    Returns:
        dict with 'total', 'with_cag', 'with_terminal_cag', 'cag_counts'
    """
    results = {
        'total': 0,
        'with_exonic_cag': 0,
        'with_terminal_cag': 0,  # CAG at positions -3, -2, -1
        'cag_counts': []
    }
    
    if not os.path.exists(fasta_path):
        return results
    
    with open(fasta_path) as f:
        seq_id = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                seq_id = line[1:]
                continue
            
            if not line:
                continue
            
            seq = line.upper()
            results['total'] += 1
            
            # For ±30bp sequences, exonic portion is first 30 bases
            # The sequence is [exon(-30 to -1)] | [intron(0 to +30)]
            exon_portion = seq[:30]  # First 30 bases are exonic
            
            # Count CAG in exonic portion
            cag_count = exon_portion.count('CAG')
            results['cag_counts'].append(cag_count)
            
            if cag_count > 0:
                results['with_exonic_cag'] += 1
            
            # Check for terminal CAG (at positions -3, -2, -1 = last 3 bases of exon)
            # In our 30-base exon portion, terminal = positions 27-29
            if exon_portion[-3:] == 'CAG':
                results['with_terminal_cag'] += 1
    
    return results


def calculate_enrichment_pvalue(test_count, test_total, control_count, control_total):
    """Calculate Fisher's exact test p-value for enrichment."""
    if test_total == 0 or control_total == 0:
        return 1.0
    
    # Contingency table
    table = [
        [test_count, test_total - test_count],
        [control_count, control_total - control_count]
    ]
    
    _, pvalue = stats.fisher_exact(table, alternative='two-sided')
    return pvalue


def generate_cag_stats_table(groups_data, output_path):
    """
    Generate CAG statistics table (PRIMARY OUTPUT).
    
    Columns: Group, Total_Sequences, N_Exonic_CAG, Pct_Exonic_CAG, 
             N_Terminal_CAG, Pct_Terminal_CAG, Pvalue_vs_Control
    """
    rows = []
    control_data = groups_data.get('Constitutive', {})
    
    for group_name, data in groups_data.items():
        total = data.get('total', 0)
        with_cag = data.get('with_exonic_cag', 0)
        with_terminal = data.get('with_terminal_cag', 0)
        
        pct_cag = (with_cag / total * 100) if total > 0 else 0
        pct_terminal = (with_terminal / total * 100) if total > 0 else 0
        
        # P-value vs control
        if group_name == 'Constitutive':
            pvalue = float('nan')
        else:
            control_total = control_data.get('total', 0)
            control_cag = control_data.get('with_exonic_cag', 0)
            pvalue = calculate_enrichment_pvalue(with_cag, total, control_cag, control_total)
        
        rows.append({
            'Group': group_name,
            'Total_Sequences': total,
            'N_Exonic_CAG': with_cag,
            'Pct_Exonic_CAG': round(pct_cag, 2),
            'N_Terminal_CAG': with_terminal,
            'Pct_Terminal_CAG': round(pct_terminal, 2),
            'Pvalue_vs_Control': round(pvalue, 6) if not pd.isna(pvalue) else 'NA'
        })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"[INFO] CAG Stats Table written to: {output_path}")
    
    # Print table to console
    print("\n" + "="*80)
    print("CAG TRINUCLEOTIDE ANALYSIS - PRIMARY METRICS")
    print("="*80)
    print(df.to_string(index=False))
    print("="*80 + "\n")
    
    return df


def run_fimo(fasta_path, motif_db, output_dir, threshold=1e-4):
    """
    Run FIMO to scan sequences against RBP motif database.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Use meme_env for FIMO
    cmd = f"mamba run -n meme_env fimo --oc {output_dir} --thresh {threshold} --verbosity 1 {motif_db} {fasta_path}"
    
    print(f"[EXEC] Running FIMO on {os.path.basename(fasta_path)}...")
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"[WARN] FIMO failed for {fasta_path}: {e}")
        return None
    
    # Parse FIMO output
    fimo_tsv = os.path.join(output_dir, "fimo.tsv")
    if os.path.exists(fimo_tsv):
        return pd.read_csv(fimo_tsv, sep='\t', comment='#')
    return None



def calculate_motif_enrichment(fimo_results_dict, groups=['UFM1_dependent', 'UFM1_independent', 'Constitutive']):
    """
    Calculate motif enrichment ratios vs control.
    
    Returns DataFrame with:
    - Motif name
    - Hits in each group
    - Log2 enrichment (Dependent/Control, Independent/Control)
    """
    if 'Constitutive' not in fimo_results_dict or fimo_results_dict['Constitutive'] is None:
        print("[WARN] No control FIMO results, skipping enrichment calculation.")
        return None
    
    # Count hits per motif per group
    motif_counts = defaultdict(lambda: defaultdict(int))
    seq_counts = {}
    
    for group, df in fimo_results_dict.items():
        if df is None or df.empty:
            seq_counts[group] = 0
            continue
        
        seq_counts[group] = df['sequence_name'].nunique()
        for motif in df['motif_id'].unique():
            motif_counts[motif][group] = len(df[df['motif_id'] == motif])
    
    # Calculate enrichment
    records = []
    control_seqs = seq_counts.get('Constitutive', 1)
    
    for motif, counts in motif_counts.items():
        control_hits = counts.get('Constitutive', 0)
        dep_hits = counts.get('UFM1_dependent', 0)
        indep_hits = counts.get('UFM1_independent', 0)
        
        dep_seqs = seq_counts.get('UFM1_dependent', 1)
        indep_seqs = seq_counts.get('UFM1_independent', 1)
        
        # Normalize by sequence count
        control_rate = (control_hits / control_seqs) if control_seqs > 0 else 0.001
        dep_rate = (dep_hits / dep_seqs) if dep_seqs > 0 else 0.001
        indep_rate = (indep_hits / indep_seqs) if indep_seqs > 0 else 0.001
        
        # Add pseudo-count to avoid log(0)
        import math
        log2_dep = math.log2(max(dep_rate, 0.001) / max(control_rate, 0.001))
        log2_indep = math.log2(max(indep_rate, 0.001) / max(control_rate, 0.001))
        
        records.append({
            'Motif': motif,
            'Hits_Dependent': dep_hits,
            'Hits_Independent': indep_hits,
            'Hits_Control': control_hits,
            'Log2_Enrichment_Dep_vs_Control': round(log2_dep, 3),
            'Log2_Enrichment_Indep_vs_Control': round(log2_indep, 3)
        })
    
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description="ESE Motif Scanning (Step B)")
    parser.add_argument("--extraction_dir", required=True, help="Directory with extracted FASTA files from Step A")
    parser.add_argument("--motif_db", required=True, help="Path to CIS-BP MEME motif database")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    groups = ['UFM1_dependent', 'UFM1_independent', 'Constitutive']
    
    # === STEP B.1: CAG ANALYSIS (PRIMARY METRIC) ===
    print("\n" + "="*80)
    print("STEP B.1: CAG TRINUCLEOTIDE ANALYSIS")
    print("="*80)
    
    cag_results = {}
    for group in groups:
        roi1_fasta = os.path.join(args.extraction_dir, f"{group}.roi1_5ss.fa")
        cag_results[group] = count_cag_in_sequences(roi1_fasta)
        print(f"[INFO] {group}: {cag_results[group]['with_exonic_cag']}/{cag_results[group]['total']} sequences with exonic CAG")
    
    # Generate CAG stats table
    cag_table_path = os.path.join(args.outdir, "CAG_Stats_Table.csv")
    generate_cag_stats_table(cag_results, cag_table_path)
    
    # === STEP B.2: FIMO RBP SCANNING ===
    print("\n" + "="*80)
    print("STEP B.2: FIMO RBP MOTIF SCANNING")
    print("="*80)
    
    fimo_results = {}
    for group in groups:
        for roi in ['roi1_5ss', 'roi2_3ss']:
            fasta_path = os.path.join(args.extraction_dir, f"{group}.{roi}.fa")
            if not os.path.exists(fasta_path):
                continue
            
            fimo_out = os.path.join(args.outdir, f"fimo_{group}_{roi}")
            result = run_fimo(fasta_path, args.motif_db, fimo_out)
            fimo_results[f"{group}_{roi}"] = result
    
    # Calculate enrichment for ROI1 (5'SS exonic region)
    roi1_results = {g: fimo_results.get(f"{g}_roi1_5ss") for g in groups}
    enrichment_df = calculate_motif_enrichment(roi1_results)
    
    if enrichment_df is not None and not enrichment_df.empty:
        enrichment_path = os.path.join(args.outdir, "RBP_Motif_Enrichment.csv")
        enrichment_df.to_csv(enrichment_path, index=False)
        print(f"[INFO] RBP Enrichment results written to: {enrichment_path}")
    
    print("\n[INFO] Step B complete.")
    print(f"[INFO] Output directory: {args.outdir}")


if __name__ == "__main__":
    main()
