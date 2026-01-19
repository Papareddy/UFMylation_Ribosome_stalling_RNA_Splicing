#!/usr/bin/env python3
"""analyze_stalling_motifs.py
Analyzes translation-coupled stalling features across multiple groups.
"""

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal, mannwhitneyu

def translate_frames(dna_seq):
    prots = []
    for i in range(3):
        if len(dna_seq) < i + 3: continue
        seq_to_trans = dna_seq[i : i + 3 * ((len(dna_seq)-i)//3)]
        if len(seq_to_trans) < 3: continue
        prots.append(str(Seq(seq_to_trans).translate()))
    return prots

def analyze_sequences_robust(fa_path):
    RARE_CODONS = {'UUA', 'UCG', 'CCG', 'GCG', 'ACG', 'GUA', 'CGU', 'CGA'}
    KD_SCALE = {'A': 1.8, 'R':-4.5, 'N':-3.5, 'D':-3.5, 'C': 2.5, 'Q':-3.5, 'E':-3.5, 'G':-0.4, 'H':-3.2, 'I': 4.5, 'L': 3.8, 'K':-3.9, 'M': 1.9, 'F': 2.8, 'P':-1.6, 'S':-0.8, 'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V': 4.2}

    metrics = []
    candidates = [] # (header, seq, reason)
    
    if not os.path.exists(fa_path): return pd.DataFrame(), []
    
    for record in SeqIO.parse(fa_path, "fasta"):
        dna = str(record.seq).upper()
        frames = translate_frames(dna)
        if not frames: continue
        
        best_stats = {}
        max_basic = -1
        best_prot = ""
        reason = ""
        
        for i, prot in enumerate(frames):
            l = len(prot)
            if l == 0: continue
            
            dna_frame = dna[i : i + 3 * l]
            codons = [dna_frame[j:j+3] for j in range(0, len(dna_frame), 3)]
            rare_count = sum(1 for c in codons if c in RARE_CODONS)
            rare_dens = rare_count / len(codons) if codons else 0
            
            hydro_score = sum(KD_SCALE.get(aa, 0) for aa in prot) / l
            
            k = prot.count('K')
            r = prot.count('R')
            p = prot.count('P')
            d = prot.count('D')
            e = prot.count('E')
            
            basic_dens = (k+r)/l
            
            # Check features
            is_poly_basic = 1 if ('KKK' in prot or 'RRR' in prot) else 0
            is_poly_pro = 1 if 'PPP' in prot else 0
            
            if basic_dens > max_basic:
                max_basic = basic_dens
                best_prot = prot
                
                # Tag reason
                r_list = []
                if is_poly_basic: r_list.append("PolyBasic")
                if is_poly_pro: r_list.append("PolyPro")
                reason = ",".join(r_list)
                
                best_stats = {
                    "Basic_Density": basic_dens,
                    "Proline_Density": p/l,
                    "Net_Charge_Density": ((k+r)-(d+e))/l,
                    "Rare_Codon_Density": rare_dens,
                    "Hydrophobicity": hydro_score,
                    "Poly_Basic": is_poly_basic,
                    "Poly_Pro": is_poly_pro,
                    "Poly_Gly": 1 if 'GGG' in prot else 0,
                    "DiPeptide_Slow": 1 if ('PP' in prot or 'PG' in prot or 'GP' in prot) else 0
                }
        
        if best_stats:
            metrics.append(best_stats)
            if best_stats["Poly_Basic"] == 1:
                candidates.append((record.id, best_prot, "PolyBasic"))
            elif best_stats["Poly_Pro"] == 1:
                candidates.append((record.id, best_prot, "PolyPro"))
            
    return pd.DataFrame(metrics), candidates

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--groups", nargs='+', required=True, help="List of 'Name:File' strings")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    all_data = []
    all_candidates = []
    
    print(f"{'Group':<20} | {'PolyBasic(%)':<12} | {'PolyPro(%)':<12}")
    print("-" * 50)
    
    for g_arg in args.groups:
        if ':' not in g_arg: continue
        name, fpath = g_arg.split(':', 1)
        
        df, cands = analyze_sequences_robust(fpath)
        if df.empty:
            print(f"[WARN] No sequences for {name}")
            continue
            
        df["Group"] = name
        all_data.append(df)
        
        # Store candidates with group name
        for cid, cseq, creason in cands:
            all_candidates.append(f">{cid}_{name} type={creason}\n{cseq}")
        
        pb = df["Poly_Basic"].mean() * 100
        pp = df["Poly_Pro"].mean() * 100
        print(f"{name:<20} | {pb:<12.2f} | {pp:<12.2f}")
        
    full_df = pd.concat(all_data, ignore_index=True)
    full_df.to_csv(os.path.join(args.outdir, "stalling_analysis_full.tsv"), sep="\t", index=False)
    
    # Save Candidates
    cand_out = os.path.join(args.outdir, "stalling_candidates.fa")
    with open(cand_out, "w") as f:
        f.writelines(c + "\n" for c in all_candidates)
    
    # Calculate Stats (Fisher's Exact Test for Poly_Basic)
    from scipy.stats import fisher_exact
    
    stats_rows = []
    groups = full_df["Group"].unique()
    
    # Updated Stats Logic: Compare against 'Genome' if present, otherwise just pairwise
    # Or specifically "UFM1_*" against "Genome"
    
    background_group = "Genome"
    if background_group in groups:
         experimental_groups = [g for g in groups if g != background_group]
         comparisons = [(g, background_group) for g in experimental_groups]
    else:
         # Fallback to simple lost vs preserved if Genome not found
         lost_g = [g for g in groups if "UFM1_dependent" in g]
         pres_g = [g for g in groups if "UFM1_independent" in g]
         comparisons = []
         for l in lost_g:
              for p in pres_g:
                   comparisons.append((l, p))
    
    for g1_name, g2_name in comparisons:
            # Get counts
            g1_data = full_df[full_df["Group"]==g1_name]["Poly_Basic"]
            g2_data = full_df[full_df["Group"]==g2_name]["Poly_Basic"]
            
            g1_poly = g1_data.sum()
            g1_tot = len(g1_data)
            g2_poly = g2_data.sum()
            g2_tot = len(g2_data)
            
            # Fisher matrix: [[Poly, NonPoly], [Poly, NonPoly]]? 
            # Or [[Group1_Poly, Group2_Poly], [Group1_Non, Group2_Non]]?
            
            # Table for fisher_exact: [[G1_Yes, G1_No], [G2_Yes, G2_No]]
            table = [[g1_poly, g1_tot - g1_poly], [g2_poly, g2_tot - g2_poly]]
            stat, pval = fisher_exact(table)
            
            stats_rows.append({
                "Group1": g1_name,
                "Group2": g2_name,
                "P_Value": pval,
                "Odds_Ratio": stat,
                "G1_Pct": (g1_poly/g1_tot)*100 if g1_tot > 0 else 0,
                "G2_Pct": (g2_poly/g2_tot)*100 if g2_tot > 0 else 0
            })
            print(f"[STATS] {g1_name} vs {g2_name}: P={pval:.2e} (OR={stat:.2f})")

    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(os.path.join(args.outdir, "stalling_stats.tsv"), sep="\t", index=False)

    # Plotting
    continuous_metrics = ["Basic_Density", "Proline_Density", "Net_Charge_Density", "Rare_Codon_Density", "Hydrophobicity"]
    binary_metrics = ["Poly_Basic", "Poly_Pro", "Poly_Gly", "DiPeptide_Slow"]
    
    # 1. Continuous Metrics (Boxplots)
    for met in continuous_metrics:
        plt.figure(figsize=(6, 5))
        sns.boxplot(data=full_df, x="Group", y=met, showfliers=False)
        sns.stripplot(data=full_df, x="Group", y=met, color="black", alpha=0.3, size=2, jitter=True)
        
        # Stats (approximate if multiple groups, but pairwise for title if 2)
        if len(stats_rows) >= 1:
             # Add stats for each comparison available in stats_rows (which now has specific pairs)
             # Just plot p-values for each pair? Too messy on title.
             # Just print them to console/log, users can see tsv.
             # Or if we have Genome, check diff against Genome for the primary group (UFM1_dependent)
             
             targets = [s for s in stats_rows if s["Group1"] == "UFM1_dependent"]
             if targets:
                 # Re-run Mann-Whitney for UFM1_dependence vs Genome
                 g1 = targets[0]["Group1"]
                 g2 = targets[0]["Group2"] # Likely Genome
                 
                 v1 = full_df[full_df["Group"]==g1][met]
                 v2 = full_df[full_df["Group"]==g2][met]
                 
                 try:
                     u_stat, u_p = mannwhitneyu(v1, v2)
                     plt.title(f"{met}\n{g1} vs {g2}\nP={u_p:.2e}")
                 except: pass

        plt.ylabel(met)
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, f"Boxplot_{met}.png"))
        plt.close()

    # 2. Binary Metrics (Barplots)
    for met in binary_metrics:
        # Calculate prevalence per group
        summary = full_df.groupby("Group")[met].mean().reset_index()
        summary["Percent"] = summary[met] * 100
        
        plt.figure(figsize=(6, 5))
        sns.barplot(data=summary, x="Group", y="Percent")
        plt.ylabel(f"% Sequence with {met}")
        
        # Calc Stats for this metric specifically - Re-run for display
        # We will use the defined 'comparisons' list from earlier
        # Assuming comparisons list is accessible here (it is local above) - better to re-derive
        
        # Re-derive comparison pairs same as above
        if "Genome" in groups:
             comp_pairs = [(g, "Genome") for g in groups if g != "Genome"]
        else:
             # fallback
             comp_pairs = [] # Skip title stats if simple logic fails or just make title generic
             pass

        if comp_pairs:
             # Pick the first one involving UFM1_dependent for the title space?
             # Or just title "Prevalence"
             relevant = [p for p in comp_pairs if p[0] == "UFM1_dependent"]
             if relevant:
                  g1, g2 = relevant[0]
                  v1 = full_df[full_df["Group"]==g1][met]
                  v2 = full_df[full_df["Group"]==g2][met]
                  
                  c1_yes = v1.sum()
                  c1_no = len(v1) - c1_yes
                  c2_yes = v2.sum()
                  c2_no = len(v2) - c2_yes
                  
                  table = [[c1_yes, c1_no], [c2_yes, c2_no]]
                  stat, pval = fisher_exact(table)
                  plt.title(f"{met}\n{g1} vs {g2}: P={pval:.2e}")

        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, f"Barplot_{met}.png"))
        plt.close()

    print(f"[DONE] Saved comparison to {args.outdir}")

if __name__ == "__main__":
    main()
