#!/usr/bin/env python3
"""
Analyze CAG enrichment in sequences with weak 5' splice sites.
Tests the hypothesis that UFM1-dependent introns have both CAG and weak 5'SS.
"""

import os
import re
from collections import defaultdict
from scipy import stats
import pandas as pd

# MaxEntScan 5'SS consensus scoring (simplified)
# Canonical 5'SS: MAG|GTRAGT (M=A/C)
# Strong: xxAG|GT
# Weak: anything else

def score_5ss_strength(seq):
    """
    Score 5' splice site strength based on the 9-mer around the junction.
    Sequence should be: EXON[-3:-1] | INTRON[0:5] = 8bp total
    For our 101bp sequences (±50bp), junction is at position 50
    Returns: 'strong' or 'weak'
    """
    if len(seq) < 101:
        return 'unknown'
    
    # Positions 47-54 (0-indexed): last 3 exon + GT + 4 intron
    junction = seq[47:56].upper()  # 9 bases around junction
    
    # Check for canonical GT
    if len(junction) >= 5 and junction[3:5] == 'GT':
        # Strong: xAG|GTxxGT or xAG|GT with consensus
        if junction[2] == 'G' and junction[1] == 'A':  # AG|GT
            return 'strong'
        elif junction[5] == 'A' and junction[7] == 'G':  # GT with AAGT consensus
            return 'strong'
        else:
            return 'weak'
    else:
        return 'non_canonical'


def has_terminal_cag(seq):
    """Check if sequence has CAG at terminal exon position (-3 to -1, positions 47-49)."""
    if len(seq) < 101:
        return False
    exon_end = seq[47:50].upper()  # Last 3 bases of exon
    return exon_end == 'CAG'


def analyze_cag_5ss_correlation(fasta_path, group_name):
    """Analyze CAG and 5'SS strength correlation."""
    results = defaultdict(int)
    
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
            ss_strength = score_5ss_strength(seq)
            has_cag = has_terminal_cag(seq)
            
            key = f"{ss_strength}_{'CAG' if has_cag else 'noCAG'}"
            results[key] += 1
            results['total'] += 1
    
    return results


def main():
    groups = {
        'UFM1_dependent': 'results/human/nucleus/step10_ese_motifs/UFM1_dependent.roi1_5ss.fa',
        'UFM1_independent': 'results/human/nucleus/step10_ese_motifs/UFM1_independent.roi1_5ss.fa',
        'Constitutive': 'results/human/nucleus/step10_ese_motifs/Constitutive.roi1_5ss.fa'
    }
    
    print('='*70)
    print('CAG + WEAK 5\' SPLICE SITE CORRELATION ANALYSIS')
    print('='*70)
    
    all_results = {}
    for group_name, fasta_path in groups.items():
        results = analyze_cag_5ss_correlation(fasta_path, group_name)
        all_results[group_name] = results
        
        total = results.get('total', 0)
        weak_cag = results.get('weak_CAG', 0)
        weak_no = results.get('weak_noCAG', 0)
        strong_cag = results.get('strong_CAG', 0)
        strong_no = results.get('strong_noCAG', 0)
        
        print(f"\n{group_name} (n={total}):")
        print(f"  Weak 5'SS + CAG:    {weak_cag:4d} ({100*weak_cag/total if total > 0 else 0:.1f}%)")
        print(f"  Weak 5'SS + no CAG: {weak_no:4d} ({100*weak_no/total if total > 0 else 0:.1f}%)")
        print(f"  Strong 5'SS + CAG:  {strong_cag:4d} ({100*strong_cag/total if total > 0 else 0:.1f}%)")
        print(f"  Strong 5'SS + no CAG: {strong_no:4d} ({100*strong_no/total if total > 0 else 0:.1f}%)")
    
    # Statistical test: Is CAG + Weak 5'SS enriched in Dependent vs Control?
    print('\n' + '='*70)
    print('STATISTICAL TEST: CAG + Weak 5\'SS Enrichment')
    print('='*70)
    
    dep = all_results['UFM1_dependent']
    ctrl = all_results['Constitutive']
    
    # Contingency table for CAG+Weak vs all others
    dep_weak_cag = dep.get('weak_CAG', 0)
    dep_other = dep.get('total', 0) - dep_weak_cag
    ctrl_weak_cag = ctrl.get('weak_CAG', 0)
    ctrl_other = ctrl.get('total', 0) - ctrl_weak_cag
    
    table = [[dep_weak_cag, dep_other], [ctrl_weak_cag, ctrl_other]]
    odds_ratio, pvalue = stats.fisher_exact(table)
    
    print(f"\nDependent: {dep_weak_cag}/{dep.get('total',0)} with Weak 5'SS + CAG")
    print(f"Control:   {ctrl_weak_cag}/{ctrl.get('total',0)} with Weak 5'SS + CAG")
    print(f"Odds Ratio: {odds_ratio:.2f}")
    print(f"P-value:    {pvalue:.2e}")
    
    if pvalue < 0.05:
        print(">>> STATISTICALLY SIGNIFICANT <<<")


if __name__ == "__main__":
    main()
