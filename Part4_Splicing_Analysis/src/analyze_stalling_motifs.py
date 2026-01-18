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
    
    # Simple pairwise: Lost vs Preserved
    # Assume groups contain "Lost" and "Preserved" substrings
    lost_g = [g for g in groups if "UFM1_dependent" in g]
    pres_g = [g for g in groups if "UFM1_independent" in g]
    
    for l_name in lost_g:
        for p_name in pres_g:
            # Get counts
            l_data = full_df[full_df["Group"]==l_name]["Poly_Basic"]
            p_data = full_df[full_df["Group"]==p_name]["Poly_Basic"]
            
            l_poly = l_data.sum()
            l_tot = len(l_data)
            p_poly = p_data.sum()
            p_tot = len(p_data)
            
            # Fisher matrix: [[Poly, NonPoly], [Poly, NonPoly]]? 
            # Or [[Group1_Poly, Group2_Poly], [Group1_Non, Group2_Non]]?
            # Scipy: [[a, b], [c, d]]
            # a: Group1 Poly, b: Group1 Non
            # c: Group2 Poly, d: Group2 Non
            
            table = [[l_poly, l_tot - l_poly], [p_poly, p_tot - p_poly]]
            stat, pval = fisher_exact(table)
            
            stats_rows.append({
                "Group1": l_name,
                "Group2": p_name,
                "P_Value": pval,
                "Odds_Ratio": stat,
                "G1_Pct": (l_poly/l_tot)*100,
                "G2_Pct": (p_poly/p_tot)*100
            })
            print(f"[STATS] {l_name} vs {p_name}: P={pval:.2e} (OR={stat:.2f})")

    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(os.path.join(args.outdir, "stalling_stats.tsv"), sep="\t", index=False)

    # Plotting
    features = ["Basic_Density", "Proline_Density", "Poly_Basic"]
    
    # Bar plot for Poly_Basic %
    plt.figure(figsize=(10, 6))
    summary = full_df.groupby("Group")["Poly_Basic"].mean().reset_index()
    
    ax = sns.barplot(data=summary, x="Group", y="Poly_Basic")
    plt.title("Poly-Basic Motif Prevalence by Group")
    plt.ylabel("Fraction of Exons with Poly-Basic Stretch")
    
    # Annotate stats on plot if simple Lost vs Preserved
    if len(stats_rows) == 1:
        row = stats_rows[0]
        plt.title(f"Poly-Basic Motif Prevalence\n{row['Group1']} vs {row['Group2']}: p={row['P_Value']:.2e}")
        
    plt.savefig(os.path.join(args.outdir, "Poly_Basic_Prevalence.png"))
    
    print(f"[DONE] Saved comparison to {args.outdir}")

if __name__ == "__main__":
    main()
