#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
from scipy.stats import fisher_exact

def main():
    parser = argparse.ArgumentParser(description="Plot SE Position Distribution (Odds Ratio)")
    parser.add_argument("--per_event_file", required=True, help="Path to per_event_compact_for_plotting.tsv")
    parser.add_argument("--output_file", required=True, help="Output plot path (PNG/PDF)")
    parser.add_argument("--background_counts", required=False, help="Path to genome background counts TSV (RI_Position | count)")
    parser.add_argument("--stats_table", required=False, help="Path to save stats table (TSV)")
    args = parser.parse_args()

    if not os.path.exists(args.per_event_file):
        print(f"[WARN] File {args.per_event_file} not found. Skipping.")
        return

    df = pd.read_csv(args.per_event_file, sep="\t")
    
    # Filter for SE events
    if "EventType" in df.columns:
        df = df[df["EventType"] == "SE"]
    
    if df.empty:
        print("[INFO] No SE events found to plot.")
        return

    # Helper: Rename RI_Position to Position if needed, or use RI_Position
    # Splicing script outputs RI_Position for compatibility.
    if "RI_Position" not in df.columns and "Event_Position" in df.columns:
         df["RI_Position"] = df["Event_Position"]

    if "RI_Position" not in df.columns:
        print("[WARN] RI_Position/Event_Position column missing. Skipping.")
        return

    valid_cats = ["Start_Codon", "Stop_Codon", "CDS", "5UTR", "3UTR", "Intronic", "NonCoding"] # Intronic rare for SE logic but possible
    
    # Calculate counts from input data
    counts = df.groupby(["dataset", "RI_Position"]).size().reset_index(name="count")
    total = df.groupby("dataset").size().reset_index(name="total")
    merged = pd.merge(counts, total, on="dataset")
    
    # --- Handle Genome Background ---
    bg_df = None
    ref_group = None
    
    if args.background_counts:
        try:
             bg_df = pd.read_csv(args.background_counts, sep="\t")
             # Normalize valid cats
             bg_df = bg_df[bg_df["RI_Position"].isin(valid_cats)]
             # Ensure numeric
             bg_df["count"] = pd.to_numeric(bg_df["count"], errors="coerce").fillna(0).astype(int)
             bg_total = bg_df["count"].sum()
             
             if bg_total > 0:
                 bg_df["dataset"] = "Genome"
                 bg_df["total"] = bg_total
                 # Append to merged logic (manually)
                 merged = pd.concat([merged, bg_df], ignore_index=True)
                 ref_group = "Genome"
        except Exception as e:
            print(f"[WARN] Failed to load background: {e}")
    
    # Fallback reference if no Genome
    if not ref_group:
        if "preserved" in merged["dataset"].unique():
            ref_group = "preserved"
            
    if not ref_group:
        print("[WARN] No reference group found (Genome or preserved). Cannot calculate Odds Ratios.")
        return
        
    print(f"[INFO] Reference Group for Odds Ratio: {ref_group}")

    # --- Comparison Logic ---
    stats_data = []
    
    # Order categories
    order = [c for c in valid_cats if c in merged["RI_Position"].unique()]
    remaining = [c for c in merged["RI_Position"].unique() if c not in valid_cats]
    order += sorted(remaining)
    
    datasets = merged["dataset"].unique()
    test_groups = [d for d in datasets if d != ref_group and not d.endswith("_psE")] # Skip psE directionals if standard available

    # Robust count retrieval
    # Re-map totals globally
    total_map = {}
    # From original DF
    t_source = df.groupby("dataset").size().to_dict()
    if ref_group == "Genome" and bg_df is not None:
        t_source["Genome"] = bg_df.iloc[0]["total"]
    
    # Counts dict: (dataset, pos) -> count
    count_map = df.groupby(["dataset", "RI_Position"]).size().to_dict()
    if ref_group == "Genome" and bg_df is not None:
        for _, r in bg_df.iterrows():
            count_map[("Genome", r["RI_Position"])] = r["count"]

    for pos in order:
        for grp in test_groups:
            c1 = count_map.get((grp, pos), 0)
            t1 = t_source.get(grp, 0)
            
            c2 = count_map.get((ref_group, pos), 0)
            t2 = t_source.get(ref_group, 0)
            
            if t1 == 0 or t2 == 0: continue
            
            # Contingency Table: [[In, Not_In], [In_Ref, Not_In_Ref]]
            table = [[c1, t1 - c1], [c2, t2 - c2]]
            odd_ratio, p_value = fisher_exact(table)
            
            sig = ""
            if p_value < 0.001: sig = "***"
            elif p_value < 0.01: sig = "**"
            elif p_value < 0.05: sig = "*"
            
            stats_data.append({
                "RI_Position": pos,
                "Group": grp,
                "Reference": ref_group,
                "Odds_Ratio": odd_ratio,
                "P_Value": p_value,
                "Significance": sig,
                "Label": f"{grp} vs {ref_group}"
            })
            
            print(f"[STATS] {pos}: {grp} vs {ref_group} -> OR={odd_ratio:.2f}, p={p_value:.2e} {sig}")

    stats_df = pd.DataFrame(stats_data)
    
    # Save Stats Table
    out_table_path = args.stats_table
    if not out_table_path:
        # Default name
        base, _ = os.path.splitext(args.output_file)
        out_table_path = base + "_Stats.tsv"
        
    stats_df.to_csv(out_table_path, sep="\t", index=False)
    print(f"[DONE] Saved stats table to {out_table_path}")

    # --- Plot Log2 Odds Ratio ---
    if stats_df.empty:
        print("[WARN] No stats generated to plot.")
        return

    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 6))
    
    # Calculate Log2 Odds Ratio
    import numpy as np
    
    # Convert to numeric
    vals = pd.to_numeric(stats_df["Odds_Ratio"], errors="coerce")
    
    plot_vals = []
    
    for v in vals:
        if v == 0:
            plot_vals.append(-6.0) # Cap at 1/64 depletion
        elif np.isinf(v):
            plot_vals.append(6.0)  # Cap at 64x enrichment
        else:
            if np.isnan(v): plot_vals.append(0.0)
            else:
                 l = np.log2(v)
                 if l > 6.0: l = 6.0
                 if l < -6.0: l = -6.0
                 plot_vals.append(l)
                 
    stats_df["Log2_OR"] = plot_vals

    ax = sns.barplot(data=stats_df, x="RI_Position", y="Log2_OR", hue="Group", order=order, palette="viridis")
    
    # Add y=0 line
    plt.axhline(0.0, color='red', linestyle='--', linewidth=1)
    
    # Annotate stars
    pos_map = {p: i for i, p in enumerate(order)}
    hue_order = sorted(stats_df["Group"].unique())
    w = 0.8 / len(hue_order)
    
    for i, row in stats_df.iterrows():
        pos = row["RI_Position"]
        if pos not in pos_map: continue
        x_idx = pos_map[pos]
        
        grp_idx = hue_order.index(row["Group"])
        x_pos = x_idx - 0.4 + w * (grp_idx + 0.5)
        # Position star
        y_val = row["Log2_OR"]
        y_pos = y_val + (0.5 if y_val >= 0 else -0.8)
        
        if y_val >= 6.0: y_pos = 6.2
        if y_val <= -6.0: y_pos = -6.5
        
        sig = row["Significance"]
        
        if sig:
            plt.text(x_pos, y_pos, sig, ha='center', va='bottom', color='black', fontsize=12, fontweight='bold')

    plt.title(f"Skipped Exon (SE) Enrichment vs {ref_group} (Log2 Odds Ratio)")
    plt.xlabel("Genomic Feature")
    plt.ylabel("Log2 Odds Ratio (Capped at +/- 6)")
    
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(args.output_file, dpi=300)
    print(f"[DONE] Saved Log2 Odds Ratio plot to {args.output_file}")

if __name__ == "__main__":
    main()
