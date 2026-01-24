#!/usr/bin/env python3
"""
ese_consensus_scanning.py
Enhanced ESE & Branch Point Scanning using Consensus Sequences

Includes:
- ESEfinder motifs (SF2/ASF, SC35, SRp40, SRp55)
- RESCUE-ESE hexamers (238 validated ESEs)
- Branch Point consensus (YNYURAY)
- CAG trinucleotide analysis
"""

import os
import re
import pandas as pd
from collections import defaultdict
from scipy import stats
import argparse

# =============================================================================
# ESE CONSENSUS SEQUENCES
# =============================================================================

# ESEfinder motifs (Cartegni et al., 2003) - simplified consensus patterns
ESEFINDER_MOTIFS = {
    'SF2/ASF': ['SRSASGA', 'RGAAGAAC', 'CARCCSAG', 'GGACG'],  # Also known as SRSF1
    'SC35': ['GRYYCSYR', 'UGCUGUU', 'AGSAGAG'],  # Also known as SRSF2
    'SRp40': ['ACGAGAGAY', 'AGGAGAA', 'TGCGTU'],  # Also known as SRSF5
    'SRp55': ['GGYRCRG', 'CUCKUCY', 'UGGAG'],  # Also known as SRSF6
}

# RESCUE-ESE hexamers (Fairbrother et al., 2002) - most validated subset
# Full list has 238; here are the most significant ones
RESCUE_ESE_HEXAMERS = [
    'GAAGAA', 'AAGAAG', 'AGAAGA', 'GAAGAG', 'AAGAAC',  # GA-rich
    'GAGGAG', 'AGGAGG', 'GGAGGA', 'GAGGAC', 'AGGAGA',  # Purine-rich
    'CAGAAG', 'AGAAGG', 'GAAGGC', 'AAGGAA', 'AGGAAG',  # Mixed
    'AAGAAA', 'AGAAAA', 'GAAAAG', 'AAAAGA', 'AAAGAA',  # A-rich ESEs
    'ACAAGA', 'CAAGAA', 'AAGAAC', 'AGAACA', 'GAACAA',  # CA-containing
    'GAAGAT', 'AAGATG', 'AGATGA', 'GATGAA', 'ATGAAG',  # GAT variants
    'CTGAAG', 'TGAAGC', 'GAAGCT', 'AAGCTG', 'AGCTGA',  # CTG variants
    'AAGACC', 'AGACCA', 'GACCAA', 'ACCAAG', 'CCAAGA',  # ACC variants
]

# Branch Point Consensus (mammalian)
# YNYURAY where Y=pyrimidine (C/T), R=purine (A/G), N=any
BRANCH_POINT_PATTERNS = [
    r'[CT][ACGT][CT][T][AG]A[CT]',  # YNYURAY (U=T in DNA)
    r'[CT]TACT[AG]',                 # More stringent consensus
    r'[CT][CT]T[AG]A[CT]',           # Relaxed pyrimidine tract
]

# =============================================================================
# SCANNING FUNCTIONS
# =============================================================================

def expand_iupac(pattern):
    """Expand IUPAC codes to regex pattern."""
    iupac = {
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
        'U': 'T',  # DNA equivalent
    }
    result = pattern.upper()
    for code, expansion in iupac.items():
        result = result.replace(code, expansion)
    return result


def count_motif_in_sequences(fasta_path, motif_patterns, region_filter=None):
    """
    Count motif occurrences in FASTA sequences.
    
    Args:
        fasta_path: Path to FASTA file
        motif_patterns: List of motif patterns (IUPAC or regex)
        region_filter: Tuple (start, end) for filtering region (0-indexed)
    
    Returns:
        dict with counts and details
    """
    results = {
        'total_sequences': 0,
        'sequences_with_motif': 0,
        'total_hits': 0,
        'hit_positions': []
    }
    
    if not os.path.exists(fasta_path):
        return results
    
    # Compile regex patterns
    compiled_patterns = []
    for pattern in motif_patterns:
        regex = expand_iupac(pattern)
        compiled_patterns.append(re.compile(regex, re.IGNORECASE))
    
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
            results['total_sequences'] += 1
            
            # Apply region filter if specified
            if region_filter:
                start, end = region_filter
                seq = seq[max(0, start):min(len(seq), end)]
            
            # Search for motifs
            found_any = False
            for pattern in compiled_patterns:
                for match in pattern.finditer(seq):
                    results['total_hits'] += 1
                    results['hit_positions'].append(match.start())
                    found_any = True
            
            if found_any:
                results['sequences_with_motif'] += 1
    
    return results


def analyze_esefinder(fasta_path, exon_region=(0, 50)):
    """Analyze ESEfinder motifs in exonic portion."""
    all_results = {}
    
    for sr_protein, motifs in ESEFINDER_MOTIFS.items():
        results = count_motif_in_sequences(fasta_path, motifs, exon_region)
        all_results[sr_protein] = results
    
    return all_results


def analyze_rescue_ese(fasta_path, exon_region=(0, 50)):
    """Analyze RESCUE-ESE hexamers in exonic portion."""
    return count_motif_in_sequences(fasta_path, RESCUE_ESE_HEXAMERS, exon_region)


def analyze_branch_point(fasta_path):
    """Analyze branch point consensus in ROI_3 (full window)."""
    return count_motif_in_sequences(fasta_path, BRANCH_POINT_PATTERNS)


def analyze_cag(fasta_path, exon_region=(0, 50)):
    """Analyze CAG trinucleotide in exonic portion."""
    results = {
        'total': 0,
        'with_exonic_cag': 0,
        'with_terminal_cag': 0,
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
            
            # Exonic portion (first 50 bases for ±50bp window)
            exon_portion = seq[exon_region[0]:exon_region[1]]
            
            cag_count = exon_portion.count('CAG')
            results['cag_counts'].append(cag_count)
            
            if cag_count > 0:
                results['with_exonic_cag'] += 1
            
            # Terminal CAG (last 3 bases of exon portion)
            if len(exon_portion) >= 3 and exon_portion[-3:] == 'CAG':
                results['with_terminal_cag'] += 1
    
    return results


def fisher_test_vs_control(test_count, test_total, control_count, control_total):
    """Calculate Fisher's exact test p-value."""
    if test_total == 0 or control_total == 0:
        return 1.0
    
    table = [
        [test_count, test_total - test_count],
        [control_count, control_total - control_count]
    ]
    
    _, pvalue = stats.fisher_exact(table, alternative='two-sided')
    return pvalue


def main():
    parser = argparse.ArgumentParser(description="Enhanced ESE & Branch Point Scanning")
    parser.add_argument("--extraction_dir", required=True, help="Directory with FASTA files")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    groups = ['UFM1_dependent', 'UFM1_independent', 'Constitutive']
    
    # Results storage
    cag_results = {}
    esefinder_results = {}
    rescue_ese_results = {}
    branch_point_results = {}
    
    print("\n" + "="*80)
    print("ENHANCED ESE & BRANCH POINT ANALYSIS")
    print("="*80)
    
    for group in groups:
        roi1_fasta = os.path.join(args.extraction_dir, f"{group}.roi1_5ss.fa")
        roi3_fasta = os.path.join(args.extraction_dir, f"{group}.roi3_bp.fa")
        
        # CAG Analysis (exonic portion of 5'SS, first 50bp)
        cag_results[group] = analyze_cag(roi1_fasta, exon_region=(0, 50))
        
        # ESEfinder Analysis
        esefinder_results[group] = analyze_esefinder(roi1_fasta, exon_region=(0, 50))
        
        # RESCUE-ESE Analysis
        rescue_ese_results[group] = analyze_rescue_ese(roi1_fasta, exon_region=(0, 50))
        
        # Branch Point Analysis
        branch_point_results[group] = analyze_branch_point(roi3_fasta)
    
    # === GENERATE REPORTS ===
    
    # 1. CAG Stats Table
    print("\n" + "-"*60)
    print("CAG TRINUCLEOTIDE ANALYSIS")
    print("-"*60)
    
    cag_rows = []
    control = cag_results.get('Constitutive', {})
    for group in groups:
        data = cag_results[group]
        total = data.get('total', 0)
        with_cag = data.get('with_exonic_cag', 0)
        with_terminal = data.get('with_terminal_cag', 0)
        
        pct_cag = (with_cag / total * 100) if total > 0 else 0
        pct_terminal = (with_terminal / total * 100) if total > 0 else 0
        
        pvalue = 'NA' if group == 'Constitutive' else round(
            fisher_test_vs_control(with_cag, total, control.get('with_exonic_cag', 0), control.get('total', 0)), 6)
        
        cag_rows.append({
            'Group': group, 'Total': total, 'N_Exonic_CAG': with_cag,
            'Pct_Exonic_CAG': round(pct_cag, 2), 'N_Terminal_CAG': with_terminal,
            'Pct_Terminal_CAG': round(pct_terminal, 2), 'Pvalue': pvalue
        })
        print(f"{group}: {with_cag}/{total} ({pct_cag:.1f}%) with exonic CAG")
    
    cag_df = pd.DataFrame(cag_rows)
    cag_df.to_csv(os.path.join(args.outdir, "CAG_Enhanced_Stats.csv"), index=False)
    
    # 2. ESEfinder Stats Table
    print("\n" + "-"*60)
    print("ESEFINDER MOTIFS")
    print("-"*60)
    
    esefinder_rows = []
    for sr_protein in ESEFINDER_MOTIFS.keys():
        for group in groups:
            data = esefinder_results[group].get(sr_protein, {})
            total = data.get('total_sequences', 0)
            with_motif = data.get('sequences_with_motif', 0)
            hits = data.get('total_hits', 0)
            
            pct = (with_motif / total * 100) if total > 0 else 0
            
            control_data = esefinder_results.get('Constitutive', {}).get(sr_protein, {})
            pvalue = 'NA' if group == 'Constitutive' else round(
                fisher_test_vs_control(with_motif, total, 
                                      control_data.get('sequences_with_motif', 0),
                                      control_data.get('total_sequences', 0)), 6)
            
            esefinder_rows.append({
                'SR_Protein': sr_protein, 'Group': group, 'Total': total,
                'N_With_Motif': with_motif, 'Pct_With_Motif': round(pct, 2),
                'Total_Hits': hits, 'Pvalue': pvalue
            })
        
        # Print summary for each SR protein
        dep = esefinder_results['UFM1_dependent'].get(sr_protein, {})
        ind = esefinder_results['UFM1_independent'].get(sr_protein, {})
        con = esefinder_results['Constitutive'].get(sr_protein, {})
        print(f"{sr_protein}: Dep={dep.get('sequences_with_motif',0)}/{dep.get('total_sequences',0)}, "
              f"Ind={ind.get('sequences_with_motif',0)}/{ind.get('total_sequences',0)}, "
              f"Ctrl={con.get('sequences_with_motif',0)}/{con.get('total_sequences',0)}")
    
    esefinder_df = pd.DataFrame(esefinder_rows)
    esefinder_df.to_csv(os.path.join(args.outdir, "ESEfinder_Stats.csv"), index=False)
    
    # 3. RESCUE-ESE Stats
    print("\n" + "-"*60)
    print("RESCUE-ESE HEXAMERS")
    print("-"*60)
    
    rescue_rows = []
    for group in groups:
        data = rescue_ese_results[group]
        total = data.get('total_sequences', 0)
        with_ese = data.get('sequences_with_motif', 0)
        hits = data.get('total_hits', 0)
        
        pct = (with_ese / total * 100) if total > 0 else 0
        
        pvalue = 'NA' if group == 'Constitutive' else round(
            fisher_test_vs_control(with_ese, total,
                                  rescue_ese_results['Constitutive'].get('sequences_with_motif', 0),
                                  rescue_ese_results['Constitutive'].get('total_sequences', 0)), 6)
        
        rescue_rows.append({
            'Group': group, 'Total': total, 'N_With_ESE': with_ese,
            'Pct_With_ESE': round(pct, 2), 'Total_Hits': hits, 'Pvalue': pvalue
        })
        print(f"{group}: {with_ese}/{total} ({pct:.1f}%) with RESCUE-ESE hexamers")
    
    rescue_df = pd.DataFrame(rescue_rows)
    rescue_df.to_csv(os.path.join(args.outdir, "RESCUE_ESE_Stats.csv"), index=False)
    
    # 4. Branch Point Stats
    print("\n" + "-"*60)
    print("BRANCH POINT CONSENSUS (YNYURAY)")
    print("-"*60)
    
    bp_rows = []
    for group in groups:
        data = branch_point_results[group]
        total = data.get('total_sequences', 0)
        with_bp = data.get('sequences_with_motif', 0)
        hits = data.get('total_hits', 0)
        
        pct = (with_bp / total * 100) if total > 0 else 0
        
        pvalue = 'NA' if group == 'Constitutive' else round(
            fisher_test_vs_control(with_bp, total,
                                  branch_point_results['Constitutive'].get('sequences_with_motif', 0),
                                  branch_point_results['Constitutive'].get('total_sequences', 0)), 6)
        
        bp_rows.append({
            'Group': group, 'Total': total, 'N_With_BP': with_bp,
            'Pct_With_BP': round(pct, 2), 'Total_Hits': hits, 'Pvalue': pvalue
        })
        print(f"{group}: {with_bp}/{total} ({pct:.1f}%) with branch point consensus")
    
    bp_df = pd.DataFrame(bp_rows)
    bp_df.to_csv(os.path.join(args.outdir, "BranchPoint_Stats.csv"), index=False)
    
    print("\n" + "="*80)
    print(f"OUTPUT FILES written to: {args.outdir}")
    print("="*80)


if __name__ == "__main__":
    main()
