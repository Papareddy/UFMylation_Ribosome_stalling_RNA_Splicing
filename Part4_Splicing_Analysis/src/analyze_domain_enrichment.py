import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
import seaborn as sns
from adjustText import adjust_text

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--enrichment", required=True)
    parser.add_argument("--protein-attributes", action="store_true", help="Plot protein attributes enrichment (expects different input format)")
    parser.add_argument("--fdr", type=float, default=0.01, help="FDR threshold for significance (default: 0.01)")
    parser.add_argument("--outfile", required=True)
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.enrichment, sep='\t')
    except Exception as e:
        print(f"[ERROR] Failed to load enrichment file: {e}")
        return

    if df.empty:
        print("[WARN] Enrichment data is empty.")
        plt.figure()
        plt.text(0.5, 0.5, "No Data", ha='center')
        plt.savefig(args.outfile)
        return

    sns.set_style("whitegrid")

    if args.protein_attributes:
        # Protein Attributes Plotting Logic
        # Expects: feature, odds_ratio, pval, log2_odds, comparison, fdr
        
        # Filter for relevant comparisons if needed
        # We want to show: Lost vs Genome, Preserved vs Genome, maybe Lost vs Preserved
        
        # Structure: Grouped Bar Chart
        # X: Feature
        # Y: Log2 Odds
        # Hue: Comparison
        
        plt.figure(figsize=(10, 6))
        
        # Rename comparisons for display
        # Rename comparisons for display
        # General formatting: Replace underscores with spaces
        df['Comparison'] = df['comparison'].str.replace('_', ' ')
        
        # Bar Plot
        # Handle Inf
        df['log2_odds'] = df['log2_odds'].replace([np.inf, -np.inf], np.nan)
        
        bar = sns.barplot(data=df, x='feature', y='log2_odds', hue='Comparison', palette='muted')
        
        plt.title('Protein Attributes Enrichment (Odds Ratio)')
        plt.ylabel('Log2 Odds Ratio')
        plt.xlabel('Feature')
        plt.axhline(0, color='black', lw=0.5)
        
        # Add significance stars
        # Iterate over bars (patches)? Seaborn makes this tricky.
        # Alternative: Iterate over data and use text coordinates.
        
        # Simple annotation loop
        # Get current axis limits
        y_min, y_max = plt.ylim()
        range_y = y_max - y_min
        
        # We need to know x-position of each bar.
        # This is complex in seaborn grouped bar.
        # Simplified: Just labels or print table? User asked for "find odds ratio", implied validation.
        # The plot shows the direction.
        
        # Let's try to annotate FDR significance
        # For now, just plot is fine. User asked to "find... and plot".
        
        plt.legend(title='Comparison', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(args.outfile, dpi=300)
        print(f"Generated biophysical plot: {args.outfile}")
        return

    # ... Original Volcano Logic ...
    # Data Prep
    df['log2_odds'] = df['log2_odds'].replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=['log2_odds', 'fdr'])
    
    min_fdr = df[df['fdr'] > 0]['fdr'].min()
    df['fdr'] = df['fdr'].replace(0, min_fdr / 10)
    df['neg_log_fdr'] = -np.log10(df['fdr'])
    
    fdr_thresh = args.fdr
    log2_thresh = 1.0
    df['significant'] = (df['fdr'] < fdr_thresh) & (df['log2_odds'].abs() > log2_thresh)

    # Plotting
    unique_sets = sorted(df['set'].unique())
    num_sets = len(unique_sets)
    
    nrows = 1
    ncols = 2
    figsize = (16, 8)
    
    if num_sets > 2:
        nrows = 2
        ncols = 2
        figsize = (16, 16)
    
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharey=True, squeeze=False)
    axes_flat = axes.flatten()
    
    colors = sns.color_palette("husl", num_sets)
    
    for i, set_name in enumerate(unique_sets):
        ax = axes_flat[i]
        col = colors[i]
        title = f"{set_name} vs Genome"
        
        subset = df[df['set'] == set_name].copy()
        
        if subset.empty:
            ax.text(0.5, 0.5, f"No Data for {set_name}", ha='center')
            continue

        # Non-sig
        ax.scatter(subset[~subset['significant']]['log2_odds'], 
                   subset[~subset['significant']]['neg_log_fdr'], 
                   color='lightgray', alpha=0.5, s=20, label='Not Significant')
        
        # Sig
        sig_subset = subset[subset['significant']]
        ax.scatter(sig_subset['log2_odds'], 
                   sig_subset['neg_log_fdr'], 
                   color=col, alpha=0.8, s=40, label=f'Significant ({set_name})')

        # Labels
        texts = []
        top_hits = sig_subset.sort_values('fdr') # Label all significant
        for _, row in top_hits.iterrows():
            texts.append(ax.text(row['log2_odds'], row['neg_log_fdr'], row['domain'], fontsize=9))
            
        if texts:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
            
        ax.set_title(title, fontsize=14)
        ax.set_xlabel('Log2 Enrichment (Odds Ratio)', fontsize=12)
        ax.axvline(x=0, linestyle='--', color='black', lw=0.5)
        ax.axhline(y=-np.log10(fdr_thresh), linestyle='--', color='red', lw=0.5)
        ax.legend()
    
    # Hide unused axes
    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis('off')

    axes_flat[0].set_ylabel('-Log10 FDR', fontsize=12)
    if nrows > 1:
        axes_flat[2].set_ylabel('-Log10 FDR', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(args.outfile, dpi=300)
    print(f"Generated plot: {args.outfile}")

if __name__ == "__main__":
    main()
