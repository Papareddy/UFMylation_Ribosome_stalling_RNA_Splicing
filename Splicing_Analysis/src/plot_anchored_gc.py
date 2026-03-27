#!/usr/bin/env python3
"""
plot_anchored_gc.py

Standalone script to plot anchored GC content data (TSV) with customizable window size.
Usage:
    python3 src/plot_anchored_gc.py path/to/data.tsv --window 500 --output my_plot.png
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Plot Anchored Intron GC Content from TSV")
    parser.add_argument("input_tsv", help="Path to anchored_intron_gc_data.tsv")
    parser.add_argument("--window", type=int, default=1000, help="Plot window size (e.g., 500 for +/- 500bp). Max depends on generated data (usually 1000).")
    parser.add_argument("--output", default="custom_anchored_gc_plot.png", help="Output plot filename")
    parser.add_argument("--ymax", type=float, default=None, help="Optional Y-axis max limit")
    parser.add_argument("--ymin", type=float, default=None, help="Optional Y-axis min limit")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_tsv):
        print(f"Error: File not found: {args.input_tsv}")
        return

    print(f"Loading {args.input_tsv}...")
    df = pd.read_csv(args.input_tsv, sep='\t')
    
    # Filter by window
    if args.window:
        print(f"Filtering data to +/- {args.window}bp...")
        df = df[(df['Position'] >= -args.window) & (df['Position'] <= args.window)]
        
    print("Plotting...")
    sns.set_theme(style="whitegrid")
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    # Define custom palette if standard 'Lost'/'Preserved' group names are present
    palette = None
    groups = df['Group'].unique()
    if 'Lost' in groups and 'Preserved' in groups and 'Constitutive' in groups:
        palette = {'Lost': 'C1', 'Preserved': 'C0', 'Constitutive': 'C2'} # Adjust colors as needed
    
    # Start Codon
    sns.lineplot(data=df[df['Feature']=='Start_Codon'], x='Position', y='GC_Percent', hue='Group', palette=palette, ax=axes[0])
    axes[0].set_title(f"Start Codon Anchored GC Content (+/- {args.window}bp)")
    axes[0].set_ylabel("GC Fraction")
    if args.ymin is not None and args.ymax is not None:
        axes[0].set_ylim(args.ymin, args.ymax)
    
    # Stop Codon
    sns.lineplot(data=df[df['Feature']=='Stop_Codon'], x='Position', y='GC_Percent', hue='Group', palette=palette, ax=axes[1])
    axes[1].set_title(f"Stop Codon Anchored GC Content (+/- {args.window}bp)")
    axes[1].set_ylabel("GC Fraction")
    if args.ymin is not None and args.ymax is not None:
        axes[1].set_ylim(args.ymin, args.ymax)
    
    plt.tight_layout()
    plt.savefig(args.output)
    print(f"Saved plot to {args.output}")

if __name__ == "__main__":
    main()
