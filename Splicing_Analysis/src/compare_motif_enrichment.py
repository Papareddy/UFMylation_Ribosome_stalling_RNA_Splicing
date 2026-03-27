
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import numpy as np

def load_ame(filepath):
    """Loads AME TSV, returns dict mapping Motif_ID -> E-value"""
    try:
        df = pd.read_csv(filepath, sep="\t", comment="#")
        # Ensure regex or string matching doesn't fail if column names vary slightly
        # Standard AME columns: rank, motif_DB, motif_ID, motif_alt_ID, consensus, p-value, adj_p-value, E-value, ...
        return df
    except Exception as e:
        print(f"[ERROR] Failed to read {filepath}: {e}")
        return pd.DataFrame()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dep_tsv", required=True, help="AME TSV for UFM1 Dependent")
    parser.add_argument("--indep_tsv", required=True, help="AME TSV for UFM1 Independent")
    parser.add_argument("--out_prefix", required=True, help="Output prefix for plots and tables")
    args = parser.parse_args()

    print(f"Loading Dependent: {args.dep_tsv}")
    df_dep = load_ame(args.dep_tsv)
    print(f"Loading Independent: {args.indep_tsv}")
    df_indep = load_ame(args.indep_tsv)

    if df_dep.empty or df_indep.empty:
        print("One or both input files are empty or invalid.")
        return

    # Select relevant columns and rename
    # We use 'motif_ID' as key. 'motif_alt_ID' is often more readable (Gene Name).
    
    def prep_df(df, suffix):
        # Create a display name: AltID if present, else ID
        df['Display'] = df.apply(lambda x: x['motif_alt_ID'] if pd.notna(x['motif_alt_ID']) and x['motif_alt_ID'] != "." else x['motif_ID'], axis=1)
        
        # Extract E-value, handle weird scientific notation if any (pandas usually handles it)
        # Handle cases where E-value might be 0 -> clip to min non-zero?
        # AME E-values can be very small.
        df['logE'] = -np.log10(df['E-value'].replace(0, 1e-300))
        
        return df[['motif_ID', 'Display', 'E-value', 'logE']].copy().add_suffix(f"_{suffix}")

    d_dep = prep_df(df_dep, "Dep")
    d_ind = prep_df(df_indep, "Indep")
    
    # Merge on motif_ID (removed suffix for merging key)
    # Actually I added suffix to everything. Let's fix.
    d_dep = d_dep.rename(columns={'motif_ID_Dep': 'motif_ID'})
    d_ind = d_ind.rename(columns={'motif_ID_Indep': 'motif_ID'})
    
    # Outer merge to keep all motifs found in either
    merged = pd.merge(d_dep, d_ind, on='motif_ID', how='outer')
    
    # Fill NA E-values with 1 (logE = 0) signifying not significant/not found
    merged['E-value_Dep'] = merged['E-value_Dep'].fillna(1000) # > threshold
    merged['logE_Dep'] = merged['logE_Dep'].fillna(0)
    
    merged['E-value_Indep'] = merged['E-value_Indep'].fillna(1000)
    merged['logE_Indep'] = merged['logE_Indep'].fillna(0)
    
    # Fill Display Name
    merged['Display'] = merged['Display_Dep'].combine_first(merged['Display_Indep'])
    
    # --- SCATTER PLOT ---
    plt.figure(figsize=(10, 10))
    sns.set_style("whitegrid")
    
    # Scatter points
    # Color by difference?
    merged['diff'] = merged['logE_Dep'] - merged['logE_Indep']
    
    scatter = sns.scatterplot(
        data=merged, 
        x='logE_Indep', 
        y='logE_Dep', 
        hue='diff', 
        palette='vlag', # Blue (Indep) -> Red (Dep)
        alpha=0.7,
        edgecolor='k',
        legend=False
    )
    
    # Diagonal line
    max_val = max(merged['logE_Dep'].max(), merged['logE_Indep'].max()) * 1.05
    plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    plt.xlim(-0.5, max_val)
    plt.ylim(-0.5, max_val)
    
    plt.xlabel("UFM1 Independent (-log10 E-value)")
    plt.ylabel("UFM1 Dependent (-log10 E-value)")
    plt.title("Motif Enrichment Comparison (vs Constitutive)")
    
    # Label top divergent motifs
    # Top Dep specific (High Diff)
    top_dep = merged.nlargest(10, 'diff')
    # Top Indep specific (Low Diff)
    top_ind = merged.nsmallest(10, 'diff')
    
    texts = []
    
    def add_lbl(row):
        # Only label if significant in respective dimension
        if row['diff'] > 0 and row['logE_Dep'] > 2:
             texts.append(plt.text(row['logE_Indep'], row['logE_Dep'], row['Display'], fontsize=8))
        elif row['diff'] < 0 and row['logE_Indep'] > 2:
             texts.append(plt.text(row['logE_Indep'], row['logE_Dep'], row['Display'], fontsize=8))

    for _, row in top_dep.iterrows(): add_lbl(row)
    for _, row in top_ind.iterrows(): add_lbl(row)
    
    # Use adjust_text if available, else standard
    try:
        from adjustText import adjust_text
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    except ImportError:
        print("adjustText not found, labels might overlap")
    
    out_pdf = args.out_prefix + "_scatter.pdf"
    plt.savefig(out_pdf)
    print(f"Scatter plot saved to {out_pdf}")
    
    # --- Save Data ---
    out_tsv = args.out_prefix + "_data.tsv"
    merged.to_csv(out_tsv, sep="\t", index=False)
    print(f"Merged data saved to {out_tsv}")

if __name__ == "__main__":
    main()
