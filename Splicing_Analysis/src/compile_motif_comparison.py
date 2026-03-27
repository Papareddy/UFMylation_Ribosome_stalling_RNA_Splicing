#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def parse_ame(fpath):
    if not os.path.exists(fpath): return None
    try:
        # AME TSV usually has headers. Skip comments if any.
        df = pd.read_csv(fpath, sep="\t", comment="#")
        # Columns: motif_ID, motif_alt_ID, consensus, p-value, adj_p-value, E-value, tests, FA_pos, TP_pos, FA_neg, TP_neg
        # We need motif_ID and E-value
        if "motif_ID" not in df.columns: return None
        return df
    except:
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lost_dir", required=True, help="Directory containing ame.tsv for Lost vs Background")
    parser.add_argument("--pres_dir", required=True, help="Directory containing ame.tsv for Preserved vs Background")
    parser.add_argument("--outfile", required=True)
    args = parser.parse_args()

    df_l = parse_ame(os.path.join(args.lost_dir, "ame.tsv"))
    df_p = parse_ame(os.path.join(args.pres_dir, "ame.tsv"))
    
    if df_l is None and df_p is None:
        print("No AME results found.")
        return

    # Normalize columns
    cols = ["motif_ID", "motif_alt_ID", "consensus", "E-value"]
    
    data_l = df_l[cols].copy() if df_l is not None else pd.DataFrame(columns=cols)
    data_p = df_p[cols].copy() if df_p is not None else pd.DataFrame(columns=cols)
    
    data_l = data_l.rename(columns={"E-value": "E_value_Lost"})
    data_p = data_p.rename(columns={"E-value": "E_value_Preserved"})
    
    # Merge outer
    merged = pd.merge(data_l, data_p, on=["motif_ID", "motif_alt_ID", "consensus"], how="outer")
    
    # Fill NA E-values with 1.0 (or max E-value observed?) -> Log10(1) = 0
    # Actually, E-values can be > 1? AME reports them.
    # Missing means not found/not significant? 
    # Or passed threshold? We set threshold high, so it should be there if tested.
    # If not present, maybe E-value is huge.
    # Let's fill with 1e5 for log plot purposes if missing?
    # Or just 10.
    
    merged["E_value_Lost"] = merged["E_value_Lost"].fillna(100)
    merged["E_value_Preserved"] = merged["E_value_Preserved"].fillna(100)
    
    # Log10
    # Avoid log(0) if E=0 (unlikely but possible if underflow). Clip min.
    merged["Log10_E_Lost"] = -np.log10(merged["E_value_Lost"] + 1e-300)
    merged["Log10_E_Preserved"] = -np.log10(merged["E_value_Preserved"] + 1e-300)
    
    merged.to_csv(args.outfile, sep="\t", index=False)
    
    # Plot
    plt.figure(figsize=(6,6))
    sns.scatterplot(data=merged, x="Log10_E_Lost", y="Log10_E_Preserved", alpha=0.7)
    
    # Add diagonal line
    lims = [
        min(min(merged["Log10_E_Lost"]), min(merged["Log10_E_Preserved"])),
        max(max(merged["Log10_E_Lost"]), max(merged["Log10_E_Preserved"]))
    ]
    plt.plot([lims[0], lims[1]], [lims[0], lims[1]], 'k--', alpha=0.5)
    
    plt.title("Motif Enrichment Comparison (-log10 E-value)")
    plt.xlabel("-log10 E-value (Lost)")
    plt.ylabel("-log10 E-value (Preserved)")
    
    plot_path = args.outfile.replace(".tsv", ".png")
    plt.savefig(plot_path)
    print(f"Comparison saved to {args.outfile} and {plot_path}")

if __name__ == "__main__":
    main()
