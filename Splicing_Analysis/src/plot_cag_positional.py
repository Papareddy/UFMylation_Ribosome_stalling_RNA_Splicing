#!/usr/bin/env python3
"""
Positional CAG density analysis across the exonic region (-50 to 0 bp).
Creates a plot comparing CAG frequency at each position for all groups.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

def count_cag_at_positions(fasta_path, window_size=50):
    """
    Count CAG occurrences at each position in the exonic region.
    For 101bp sequences (±50bp around junction), exon is positions 0-49.
    """
    position_counts = defaultdict(int)
    total_sequences = 0
    
    if not os.path.exists(fasta_path):
        return position_counts, 0
    
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>') or not line:
                continue
            
            seq = line.upper()
            total_sequences += 1
            
            # Scan exonic region (first 50bp) for CAG
            exon = seq[:window_size]
            for i in range(len(exon) - 2):
                if exon[i:i+3] == 'CAG':
                    position_counts[i] += 1  # Position relative to -50
    
    return position_counts, total_sequences


def main():
    groups = {
        'UFM1_dependent': 'results/human/nucleus/step10_ese_motifs/UFM1_dependent.roi1_5ss.fa',
        'UFM1_independent': 'results/human/nucleus/step10_ese_motifs/UFM1_independent.roi1_5ss.fa',
        'Constitutive': 'results/human/nucleus/step10_ese_motifs/Constitutive.roi1_5ss.fa'
    }
    
    colors = {'UFM1_dependent': '#E41A1C', 'UFM1_independent': '#377EB8', 'Constitutive': 'gray'}
    
    # Create figure
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    
    # === Panel A: CAG frequency at each position ===
    ax1 = axes[0]
    
    all_data = {}
    for group_name, fasta_path in groups.items():
        counts, total = count_cag_at_positions(fasta_path)
        
        # Convert to frequency (normalize by total sequences)
        positions = list(range(48))  # positions 0-47 (CAG can start at 0-47 to fit in 50bp)
        freq = [counts.get(i, 0) / total * 100 if total > 0 else 0 for i in positions]
        
        all_data[group_name] = {'counts': counts, 'total': total, 'freq': freq}
        
        # Convert position to relative (-50 to -3)
        x_positions = [i - 50 for i in positions]
        
        ax1.plot(x_positions, freq, label=group_name, color=colors[group_name], linewidth=2, alpha=0.8)
    
    ax1.axvline(x=-3, color='black', linestyle='--', linewidth=1, label='Terminal CAG position')
    ax1.set_xlabel('Position relative to 5\' splice site (bp)', fontsize=12)
    ax1.set_ylabel('% Sequences with CAG', fontsize=12)
    ax1.set_title('CAG Trinucleotide Position in Exon (5\'SS ±50bp)', fontsize=14, fontweight='bold')
    ax1.legend(loc='upper left')
    ax1.set_xlim(-50, 0)
    ax1.grid(True, alpha=0.3)
    
    # === Panel B: Cumulative CAG count (bar plot for total in exon) ===
    ax2 = axes[1]
    
    group_names = list(groups.keys())
    totals = [all_data[g]['total'] for g in group_names]
    cag_totals = [sum(all_data[g]['counts'].values()) for g in group_names]
    avg_cag_per_seq = [cag_totals[i] / totals[i] if totals[i] > 0 else 0 for i in range(len(group_names))]
    
    bars = ax2.bar(group_names, avg_cag_per_seq, color=[colors[g] for g in group_names], edgecolor='black')
    
    for bar, val in zip(bars, avg_cag_per_seq):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, 
                f'{val:.2f}', ha='center', va='bottom', fontsize=11)
    
    ax2.set_ylabel('Average CAG count per sequence', fontsize=12)
    ax2.set_title('Total CAG Content in Exonic Region (-50 to 0 bp)', fontsize=14, fontweight='bold')
    ax2.set_ylim(0, max(avg_cag_per_seq) * 1.2)
    
    plt.tight_layout()
    plt.savefig('results/human/nucleus/step10_ese_motifs/plots/CAG_Positional_Density.pdf', dpi=150)
    plt.savefig('results/human/nucleus/step10_ese_motifs/plots/CAG_Positional_Density.png', dpi=150)
    print('[INFO] CAG positional density plot saved.')
    
    # Print summary
    print('\n' + '='*60)
    print('CAG POSITIONAL ANALYSIS SUMMARY')
    print('='*60)
    for g in group_names:
        print(f"{g}: {avg_cag_per_seq[group_names.index(g)]:.2f} CAG per sequence (n={all_data[g]['total']})")


if __name__ == "__main__":
    main()
