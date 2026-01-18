#!/usr/bin/env python3
"""analyze_se_ri_adjacency.py
Analyzes the spatial relationship and dPSI correlation between Skipped Exon (SE) and Retained Intron (RI) events.
Specifically tests if Lost RI (dPSI < 0) is adjacent to Gained SE (dPSI > 0).
"""

import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, fisher_exact

def load_data(fpaths, event_type):
    dfs = []
    if isinstance(fpaths, str): fpaths = [fpaths]
    for f in fpaths:
        if os.path.exists(f):
            d = pd.read_csv(f, sep='\t')
            if "EventType" in d.columns: d = d[d["EventType"] == event_type]
            dfs.append(d)
    if not dfs: return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)

def get_col(row, base):
    if base in row: return row[base]
    if f"{base}.WT" in row: return row[f"{base}.WT"]
    if f"{base}.UFM" in row: return row[f"{base}.UFM"]
    # If merged from prepare_rmats_data.R, coordinates might be suffixed
    raise KeyError(f"Column {base} not found (checked .WT, .UFM)")

def parse_coords(row, etype):
    chrom = get_col(row, "chr")
    strand = get_col(row, "strand")
    
    if etype == "RI":
        uEE = int(get_col(row, "upstreamEE"))
        dES = int(get_col(row, "downstreamES"))
        dpsi = row.get("dPSI_num.WT") if "dPSI_num.WT" in row else row.get("dPSI_num")
        return chrom, strand, uEE, dES, dpsi
        
    elif etype == "SE":
        c_start = int(get_col(row, "exonStart_0base"))
        c_end = int(get_col(row, "exonEnd"))
        u_ee = int(get_col(row, "upstreamEE"))
        d_es = int(get_col(row, "downstreamES"))
        dpsi = row.get("dPSI_num.WT") if "dPSI_num.WT" in row else row.get("dPSI_num")
        return chrom, strand, (u_ee, c_start), (c_end, d_es), dpsi

def check_overlap(iv1, iv2):
    # iv = (start, end)
    # Overlap if max(start1, start2) < min(end1, end2)
    s = max(iv1[0], iv2[0])
    e = min(iv1[1], iv2[1])
    return (e - s) > 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--se_files", nargs='+', required=True, help="List of SE files")
    parser.add_argument("--ri_files", nargs='+', required=True, help="List of RI files")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    se_df = load_data(args.se_files, "SE")
    ri_df = load_data(args.ri_files, "RI")
    
    print(f"[INFO] Loaded {len(se_df)} SE events and {len(ri_df)} RI events.")
    
    if not se_df.empty:
        print("[DEBUG] Sample SE Row:")
        print(se_df.iloc[0])
        print("[DEBUG] Parsed SE Coords:", parse_coords(se_df.iloc[0], "SE"))
        
    if not ri_df.empty:
        print("[DEBUG] Sample RI Row:")
        print(ri_df.iloc[0])
        print("[DEBUG] Parsed RI Coords:", parse_coords(ri_df.iloc[0], "RI"))
    
    # Simple Interval Tree or just loop (datasets are small < 2000)
    # Loop acceptable.
    
    # Store RI by chrom
    ri_by_chr = {}
    for idx, row in ri_df.iterrows():
        try:
            chrom, strand, uEE, dES, dpsi = parse_coords(row, "RI")
            if chrom not in ri_by_chr: ri_by_chr[chrom] = []
            ri_by_chr[chrom].append({
                "iv": (uEE, dES),
                "dpsi": dpsi,
                "id": row.get("event_id"),
                "strand": strand
            })
        except: pass
        
    adjacency = []
    
    for idx, row in se_df.iterrows():
        try:
            chrom, strand, i1, i2, se_dpsi = parse_coords(row, "SE")
            if chrom not in ri_by_chr: continue
            
            # Check against all RIs on this chrom
            # Could optimize with bisect if needed, but N is small
            for ri_ev in ri_by_chr[chrom]:
                # Check Upstream Intron Overlap
                is_adj = False
                rel = ""
                
                # Loose overlap check
                if check_overlap(i1, ri_ev["iv"]):
                    rel = "Upstream_overlap"
                    is_adj = True
                elif check_overlap(i2, ri_ev["iv"]):
                    rel = "Downstream_overlap"
                    is_adj = True
                    
                if is_adj:
                    adjacency.append({
                         "SE_ID": row.get("event_id"),
                         "RI_ID": ri_ev["id"],
                         "SE_dPSI": se_dpsi,
                         "RI_dPSI": ri_ev["dpsi"],
                         "Relation": rel,
                         "Strand": strand
                    })
        except: pass
        
    res_df = pd.DataFrame(adjacency)
    out_tsv = os.path.join(args.outdir, "SE_RI_Adjacency.tsv")
    res_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[INFO] Found {len(res_df)} overlapping SE-RI pairs. Saved to {out_tsv}")
    
    if res_df.empty:
        print("[WARN] No adjacency found. Exiting.")
        return

    # --- Analysis ---
    # Classify Directions
    # SE Pos (>0.05), SE Neg (<-0.05)
    # RI Pos (>0.05), RI Neg (<-0.05)
    
    def get_cat(x):
        if x > 0.05: return "Pos"
        if x < -0.05: return "Neg"
        return "Neu"
        
    res_df["SE_Dir"] = res_df["SE_dPSI"].apply(get_cat)
    res_df["RI_Dir"] = res_df["RI_dPSI"].apply(get_cat)
    
    # 1. Contingency Table
    ct = pd.crosstab(res_df["SE_Dir"], res_df["RI_Dir"])
    ct_path = os.path.join(args.outdir, "SE_RI_Direction_Counts.tsv")
    ct.to_csv(ct_path, sep="\t")
    print("[INFO] Direction Counts:")
    print(ct)
    
    # 2. Scatter Plot
    plt.figure(figsize=(8, 8))
    sns.scatterplot(data=res_df, x="SE_dPSI", y="RI_dPSI", hue="Relation", alpha=0.6)
    plt.axhline(0, color='grey', linestyle='--')
    plt.axvline(0, color='grey', linestyle='--')
    plt.title("Correlation of adjacent SE and RI events")
    
    # Calc Correlation
    clean = res_df.dropna(subset=["SE_dPSI", "RI_dPSI"])
    if len(clean) > 2:
        r, p = pearsonr(clean["SE_dPSI"], clean["RI_dPSI"])
        plt.annotate(f"r={r:.2f}, p={p:.2e}", xy=(0.05, 0.95), xycoords='axes fraction')
        
    plot_path = os.path.join(args.outdir, "SE_RI_Correlation.png")
    plt.savefig(plot_path)
    print(f"[DONE] Saved plot to {plot_path}")
    
    # 3. Hypothesis Check: RI Neg + SE Pos
    # Filter for SE_Pos and RI_Neg
    target = res_df[(res_df["SE_Dir"] == "Pos") & (res_df["RI_Dir"] == "Neg")]
    print(f"[INFO] Found {len(target)} pairs satisfying 'RI Neg (Lost) + SE Pos (Gained)'.")
    if not target.empty:
        target_path = os.path.join(args.outdir, "SE_Pos_RI_Neg_Events.tsv")
        target.to_csv(target_path, sep="\t", index=False)
        print(f"Saved target events to {target_path}")

if __name__ == "__main__":
    main()
