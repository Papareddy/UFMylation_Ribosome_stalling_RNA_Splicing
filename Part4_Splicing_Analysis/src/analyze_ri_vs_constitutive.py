#!/usr/bin/env python3
"""analyze_ri_vs_constitutive.py
Comparisons of RI Introns (Lost/Preserved) against a Constitutive Intron Background.
Uses:
1. generate_constitutive_introns.py (to get bg BED)
2. bedtools getfasta (to get bg FASTA)
"""

import os
import sys
import argparse
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ranksums, mannwhitneyu
import pandas as pd
from Bio import SeqIO

def run_cmd(cmd, description):
    print(f"[EXEC] {description}...")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed: {description}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--genome_fasta", required=True)
    parser.add_argument("--lost_tsv", required=True)
    parser.add_argument("--preserved_tsv", required=True)
    parser.add_argument("--lost_fa", required=True, help="FASTA of lost introns (from extract_ri_motifs.py)")
    parser.add_argument("--preserved_fa", required=True, help="FASTA of preserved introns")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--script_dir", required=True, help="Directory containing generate_constitutive_introns.py")
    parser.add_argument("--motif_db", help="Path to MEME Motif DB for AME analysis", required=False)
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # Files
    const_bed = os.path.join(args.outdir, "constitutive_introns.bed")
    const_fa = os.path.join(args.outdir, "constitutive_introns.fa")
    
    # 1. Generate Constitutive Background
    if not os.path.exists(const_bed):
        gen_script = os.path.join(args.script_dir, "generate_constitutive_introns.py")
        cmd_gen = f"python3 {gen_script} --gtf {args.gtf} --out_bed {const_bed} --exclude_files {args.lost_tsv} {args.preserved_tsv} --n_sample 5000"
        run_cmd(cmd_gen, "Generating Constitutive Intron Background")
    
    # 2. Extract Background Sequences (Intron, 5'SS, 3'SS, Wide Variants)
    kinds = ['intron', '5ss', '3ss', '5ss_wide', '3ss_wide']
    
    # BED paths inferred from generate_constitutive_introns logic
    # It produced base_out.[kind].bed
    base_out = const_bed.replace(".bed", "")
    if base_out.endswith(".intron"): base_out = base_out.replace(".intron", "")
    
    beds = {k: f"{base_out}.{k}.bed" for k in kinds}
    fastas = {k: f"{base_out}.{k}.fa" for k in kinds}
    
    for k in kinds:
        if not os.path.exists(fastas[k]):
            if os.path.exists(beds[k]):
                cmd_fa = f"bedtools getfasta -s -fi {args.genome_fasta} -bed {beds[k]} -fo {fastas[k]} -name"
                run_cmd(cmd_fa, f"Extracting Constitutive {k} sequences")
    
    # 3. Score Background SS (MaxEntScan)
    # 5ss
    score_ss(fastas['5ss'], '5ss')
    score_ss(fastas['3ss'], '3ss')

    # --- Deep Dive Analyses ---
    
    # ... (Previous Extract/Score Steps)
    
    # --- Deep Dive Analyses ---
    stats_rows = []
    
    # 5A: MaxEnt Comparison (Wilcoxon)
    # Need .scores from Lost/Preserved (assumed in same dir as input fastas)
    # Note: extract_ri_motifs.py puts them alongside input fasta.
    lost_5ss = args.lost_fa + ".scores" # extract_ri_motifs convention: fasta.scores
    lost_3ss = args.lost_fa.replace("intron.fa", "3ss.fa.scores") 
    pres_5ss = args.preserved_fa.replace("intron.fa", "5ss.fa.scores")
    pres_3ss = args.preserved_fa.replace("intron.fa", "3ss.fa.scores")
    
    # Fix pathing logic if needed. extract_ri_motifs usually makes:
    # prefix.5ss.fa, prefix.5ss.fa.scores
    # args.lost_fa is "prefix.intron.fa".
    # So "prefix.5ss.fa.scores" is correct.
    lost_5ss = args.lost_fa.replace(".intron.fa", ".5ss.fa.scores")
    lost_3ss = args.lost_fa.replace(".intron.fa", ".3ss.fa.scores")
    pres_5ss = args.preserved_fa.replace(".intron.fa", ".5ss.fa.scores")
    pres_3ss = args.preserved_fa.replace(".intron.fa", ".3ss.fa.scores")

    const_5ss = fastas['5ss'] + ".scores"
    const_3ss = fastas['3ss'] + ".scores"
    
    print("[INFO] Analysis 5A: MaxEntScan Statistical Comparison...")
    stats_rows.extend(compare_scores(lost_5ss, const_5ss, "UFM1_dependent", "Constitutive", "5'SS MaxEnt Score", args.outdir))
    stats_rows.extend(compare_scores(lost_3ss, const_3ss, "UFM1_dependent", "Constitutive", "3'SS MaxEnt Score", args.outdir))
    stats_rows.extend(compare_scores(pres_5ss, const_5ss, "UFM1_independent", "Constitutive", "5'SS MaxEnt Score", args.outdir))
    stats_rows.extend(compare_scores(pres_3ss, const_3ss, "UFM1_independent", "Constitutive", "3'SS MaxEnt Score", args.outdir))

    # 5C: GC Content & Length
    print("[INFO] Analysis 5C: GC Content & Length Analysis...")
    stats_rows.extend(compare_seq_features(args.lost_fa, fastas['intron'], "UFM1_dependent", "Constitutive", args.outdir))
    stats_rows.extend(compare_seq_features(args.preserved_fa, fastas['intron'], "UFM1_independent", "Constitutive", args.outdir))

    # Save Stats
    stats_df = pd.DataFrame(stats_rows)
    stats_out = os.path.join(args.outdir, "summary_stats.tsv")
    stats_df.to_csv(stats_out, sep="\t", index=False)
    print(f"[INFO] Summary Stats saved to {stats_out}")

    # --- Consolidated Feature Analysis ---
    print("[INFO] Generating Master Feature File...")
    
    # helper to load features
    def load_feat(fa_path, kind):
        if not os.path.exists(fa_path): return {}
        # load length, gc
        lens, gcs = get_gc_len(fa_path)
        # load scores if they exist
        scores5 = load_scores(fa_path.replace(f".{kind}.fa", ".5ss.fa.scores")) 
        scores3 = load_scores(fa_path.replace(f".{kind}.fa", ".3ss.fa.scores"))
        
        # We need to map them by index if possible, or just lists?
        # The user wants distributions. Lists are fine if we just concat.
        # But for "Master File" usually implies a table with ID. 
        # get_gc_len reads fasta sequentially. 
        # score files are sequential.
        # So we can assume order is preserved.
        
        # MaxEntScan output might have fewer lines if errors? Assuming 1-to-1.
        data = []
        for i in range(len(lens)):
            s5 = scores5[i] if i < len(scores5) else "NaN"
            s3 = scores3[i] if i < len(scores3) else "NaN"
            data.append({
                "Length": lens[i],
                "GC": gcs[i],
                "MaxEnt5": s5, 
                "MaxEnt3": s3
            })
        return data

    master_data = []
    
    # 1. Constitutive
    c_data = load_feat(fastas['intron'], "intron")
    for d in c_data: d["Group"] = "Constitutive"; master_data.append(d)
    
    # 2. Lost -> UFM1_dependent
    l_data = load_feat(args.lost_fa, "intron")
    for d in l_data: d["Group"] = "UFM1_dependent"; master_data.append(d)
    
    # 3. Preserved -> UFM1_independent
    p_data = load_feat(args.preserved_fa, "intron")
    for d in p_data: d["Group"] = "UFM1_independent"; master_data.append(d)
    
    master_df = pd.DataFrame(master_data)
    master_tsv = os.path.join(args.outdir, "master_features.tsv")
    master_df.to_csv(master_tsv, sep="\t", index=False)
    
    # Call R Plotter
    print("[INFO] Plotting Consolidated Features...")
    plot_r = os.path.join(args.script_dir, "plot_features.R")
    plot_out = os.path.join(args.outdir, "combined_features.pdf")
    cmd_plot = f"mamba run -n splicing-functional Rscript {plot_r} --input={master_tsv} --out={plot_out}"
    try:
        run_cmd(cmd_plot, "Plotting Features")
    except:
        print("[WARN] Feature plotting failed.")

    # 5B: AME (Known Motifs) - Loop over regions
    regions = ['intron', '5ss', '3ss']
    
    if args.motif_db and os.path.exists(args.motif_db):
        for r in regions:
            print(f"[INFO] Analysis 5B: Running AME for {r} with {args.motif_db}...")
            
            # Infer paths
            # args.lost_fa defines the prefix mostly.
            # lost_fa input is "path/lost.intron.fa"
            # so we replace ".intron.fa" with ".{r}.fa"
            
            # Determine Wide variant for SS
            if r in ['5ss', '3ss']:
                r_use = f"{r}_wide"
            else:
                r_use = r
                
            l_fa = args.lost_fa.replace(".intron.fa", f".{r_use}.fa")
            p_fa = args.preserved_fa.replace(".intron.fa", f".{r_use}.fa")
            bg_fa = fastas[r_use] # constitutive
            
            if not os.path.exists(l_fa) or not os.path.exists(p_fa) or not os.path.exists(bg_fa):
                print(f"[WARN] Missing input files for AME region {r} (using {r_use}). Skipping.")
                continue
                
            # Run AME - Dependent
            out_ame_l = os.path.join(args.outdir, f"ame_{r}_UFM1_dependent_vs_constitutive")
            cmd_ame_l = f"mamba run -n meme_env ame --verbose 1 --evalue-report-threshold 1000 --control {bg_fa} --oc {out_ame_l} {l_fa} {args.motif_db}"
            try: run_cmd(cmd_ame_l, f"AME {r}: UFM1_dependent vs Constitutive")
            except: print(f"[WARN] AME {r} UFM1_dependent failed.")
            
            # Run AME - Independent
            out_ame_p = os.path.join(args.outdir, f"ame_{r}_UFM1_independent_vs_constitutive")
            cmd_ame_p = f"mamba run -n meme_env ame --verbose 1 --evalue-report-threshold 1000 --control {bg_fa} --oc {out_ame_p} {p_fa} {args.motif_db}"
            try: run_cmd(cmd_ame_p, f"AME {r}: UFM1_independent vs Constitutive")
            except: print(f"[WARN] AME {r} UFM1_independent failed.")
            
            # Plot Comparison
            tsv_l = os.path.join(out_ame_l, "ame.tsv")
            tsv_p = os.path.join(out_ame_p, "ame.tsv")
            
            if os.path.exists(tsv_l) and os.path.exists(tsv_p):
                print(f"[INFO] Generating Motif Comparison Scatter Plot for {r}...")
                plot_r_script = os.path.join(args.script_dir, "plot_ame_comparison.R")
                plot_pdf = os.path.join(args.outdir, f"motif_enrichment_comparison_{r}.pdf")
                
                if os.path.exists(plot_r_script):
                    cmd_plot_r = f"Rscript {plot_r_script} --dep={tsv_l} --indep={tsv_p} --out={plot_pdf}"
                    try:
                        run_cmd(cmd_plot_r, f"Plotting Motif Comparison ({r})")
                    except:
                        print(f"[WARN] Plotting failed for {r}.")

    elif args.motif_db:
        print(f"[WARN] Motif DB not found at {args.motif_db}. Skipping AME.")


# ... (Imports) ...

def compare_scores(f1, f2, l1, l2, feature_name, outdir):
    s1 = load_scores(f1)
    s2 = load_scores(f2)
    if not s1 or not s2: return []
    
    stat, p = ranksums(s1, s2)
    mean1, mean2 = pd.Series(s1).mean(), pd.Series(s2).mean()
    print(f"[STATS] {feature_name}: {l1} ({mean1:.2f}) vs {l2} ({mean2:.2f}) -> p={p:.2e}")
    
    # Return row for table
    row = {
        "Feature": feature_name,
        "Group1": l1, "Mean1": mean1,
        "Group2": l2, "Mean2": mean2,
        "P_Value": p
    }
    
    # Plotting code (simplified or preserved)
    # ... (Keep plotting) ...
    res = {l1: s1, l2: s2}
    data = []
    for k, v in res.items():
        for x in v: data.append({"Group": k, "Score": x})
    df = pd.DataFrame(data)
    plt.figure(figsize=(6, 4))
    sns.boxplot(data=df, x="Group", y="Score")
    plt.title(f"{feature_name}\np={p:.2e}")
    plt.savefig(os.path.join(outdir, f"Comparison_{feature_name.replace(' ','_')}_{l1}_vs_{l2}.png"))
    plt.close()
    
    return [row]

def compare_seq_features(f1, f2, l1, l2, outdir):
    len1, gc1 = get_gc_len(f1)
    len2, gc2 = get_gc_len(f2)
    if not len1 or not len2: return []
    
    rows = []
    
    # Length Stats
    s, p_len = mannwhitneyu(len1, len2)
    m_l1, m_l2 = pd.Series(len1).mean(), pd.Series(len2).mean()
    rows.append({
        "Feature": "Length (bp)",
        "Group1": l1, "Mean1": m_l1,
        "Group2": l2, "Mean2": m_l2,
        "P_Value": p_len
    })
    
    # GC Stats
    s, p_gc = mannwhitneyu(gc1, gc2)
    m_g1, m_g2 = pd.Series(gc1).mean(), pd.Series(gc2).mean()
    rows.append({
        "Feature": "GC Content",
        "Group1": l1, "Mean1": m_g1,
        "Group2": l2, "Mean2": m_g2,
        "P_Value": p_gc
    })
    
    # Plots (Length log scale)
    df_len = pd.DataFrame([{"Group": l1, "Length": x} for x in len1] + [{"Group": l2, "Length": x} for x in len2])
    plt.figure(figsize=(6, 4))
    sns.boxplot(data=df_len, x="Group", y="Length")
    plt.yscale("log")
    plt.title(f"Length (bp)\np={p_len:.2e}")
    plt.savefig(os.path.join(outdir, f"Comparison_Length_{l1}_vs_{l2}.png"))
    plt.close()

    # Plots (GC)
    df_gc = pd.DataFrame([{"Group": l1, "GC": x} for x in gc1] + [{"Group": l2, "GC": x} for x in gc2])
    plt.figure(figsize=(6, 4))
    sns.boxplot(data=df_gc, x="Group", y="GC")
    plt.title(f"GC Content\np={p_gc:.2e}")
    plt.savefig(os.path.join(outdir, f"Comparison_GC_{l1}_vs_{l2}.png"))
    plt.close()

    return rows

# New score_ss using native python scorer
import maxent_scorer

def score_ss(fa_file, kind):
    if not os.path.exists(fa_file): return
    out_file = fa_file + ".scores"
    if os.path.exists(out_file): return
    
    print(f"[INFO] Scoring {kind} using Native MaxEnt Scorer...")
    try:
        with open(fa_file) as f, open(out_file, "w") as out:
            for line in f:
                if line.startswith(">"): continue
                seq = line.strip()
                if not seq: continue
                
                # Clean seq just in case
                score = 0
                if kind == '5ss':
                    score = maxent_scorer.score5(seq)
                else:
                    score = maxent_scorer.score3(seq)
                
                out.write(f"{seq}\t{score:.2f}\n")
    except Exception as e:
        print(f"[ERROR] Scoring failed for {fa_file}: {e}")

def load_scores(fpath):
    if not os.path.exists(fpath): return []
    try:
        # MaxEnt output: sequence \t score
        df = pd.read_csv(fpath, sep="\t", header=None, names=["seq", "score"])
        return df["score"].dropna().tolist()
    except: return []

def get_gc_len(fpath):
    lens = []
    gcs = []
    if not os.path.exists(fpath): return [], []
    for rec in SeqIO.parse(fpath, "fasta"):
        seq = str(rec.seq).upper()
        l = len(seq)
        if l == 0: continue
        gc = (seq.count("G") + seq.count("C")) / l
        lens.append(l)
        gcs.append(gc)
    return lens, gcs
if __name__ == "__main__":
    main()
