#!/usr/bin/env python3
"""analyze_se_vs_constitutive.py
Comparisons of SE Exons (Lost/Preserved) against a Constitutive Exon Background.
Analyze regions: Exon, 3'SS, 5'SS, Upstream Intron, Downstream Intron.
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
        # sys.exit(1) # Don't hard exit, allow partial completion

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--genome_fasta", required=True)
    parser.add_argument("--lost_tsv", required=True)
    parser.add_argument("--preserved_tsv", required=True)
    # The fa_prefix arg assumes standard naming from extract_se_motifs.py: prefix.exon.fa, prefix.3ss.fa etc.
    parser.add_argument("--lost_prefix", required=True, help="Prefix of lost FASTA files (e.g., .../lost)")
    parser.add_argument("--preserved_prefix", required=True, help="Prefix of preserved FASTA files")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--script_dir", required=True, help="Directory containing generate_constitutive_exons.py")
    parser.add_argument("--motif_db", help="Path to MEME Motif DB", required=False)
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # 1. Generate Constitutive Background
    const_bed = os.path.join(args.outdir, "constitutive_exons.bed")
    
    if not os.path.exists(const_bed):
        gen_script = os.path.join(args.script_dir, "generate_constitutive_exons.py")
        cmd_gen = f"python3 {gen_script} --gtf {args.gtf} --out_bed {const_bed} --exclude_files {args.lost_tsv} {args.preserved_tsv} --n_sample 5000"
        run_cmd(cmd_gen, "Generating Constitutive Exon Background")
    
    # 2. Extract Background Sequences
    regions = ['exon', '3ss', '5ss', 'intron_upstream', 'intron_downstream', '5ss_wide', '3ss_wide', 'entire_region']
    base_out = const_bed.replace(".bed", "")
    if base_out.endswith(".exon"): base_out = base_out.replace(".exon", "")
    
    const_fastas = {}
    for r in regions:
        bed = f"{base_out}.{r}.bed"
        fa = f"{base_out}.{r}.fa"
        const_fastas[r] = fa
        if not os.path.exists(fa) and os.path.exists(bed):
            cmd_fa = f"bedtools getfasta -s -fi {args.genome_fasta} -bed {bed} -fo {fa} -name"
            run_cmd(cmd_fa, f"Extracting Constitutive {r} sequences")

    # 3. Analyze Features (GC, Length) - Exon Only
    # Compare Lost Exons vs Constitutive, Preserved vs Constitutive
    lost_exon_fa = f"{args.lost_prefix}.exon.fa"
    pres_exon_fa = f"{args.preserved_prefix}.exon.fa"
    
    if os.path.exists(lost_exon_fa) and os.path.exists(const_fastas['exon']):
        compare_features(lost_exon_fa, const_fastas['exon'], "UFM1_dependent", "Constitutive", args.outdir, "Exon")
    if os.path.exists(pres_exon_fa) and os.path.exists(const_fastas['exon']):
        compare_features(pres_exon_fa, const_fastas['exon'], "UFM1_independent", "Constitutive", args.outdir, "Exon")

    # 4. AME Analysis
    if args.motif_db:
        for r in regions:
            # Skip wide variants in main loop if we handle mapping?
            # Or iterate all? The user wants plots for "3ss" and "5ss".
            # We should perform AME on "3ss_wide" but label it as "3ss"?
            # Or just iterate the original regions and map the filenames.
            
            if "_wide" in r: continue # Skip iterating explicitly if we map.
            
            print(f"[INFO] Running AME for Region: {r}")
            
            # Map valid AME regions to wide counterparts if needed
            r_use = r
            if r in ['3ss', '5ss']: r_use = f"{r}_wide"
                
            # Lost (Dependent)
            l_fa = f"{args.lost_prefix}.{r_use}.fa"
            if os.path.exists(l_fa) and os.path.exists(const_fastas[r_use]):
                run_ame_and_plot(l_fa, const_fastas[r_use], args.motif_db, "UFM1_dependent", "Constitutive", r, args.outdir, args.script_dir)
            
            # Preserved (Independent)
            p_fa = f"{args.preserved_prefix}.{r_use}.fa"
            if os.path.exists(p_fa) and os.path.exists(const_fastas[r_use]):
                run_ame_and_plot(p_fa, const_fastas[r_use], args.motif_db, "UFM1_independent", "Constitutive", r, args.outdir, args.script_dir)
                
        # Consolidated Plot (Dependent vs Independent) for each region
        # This mirrors the logic in analyze_ri_vs_constitutive
        for r in regions:
            ame_dir_dep = os.path.join(args.outdir, f"ame_{r}_UFM1_dependent_vs_Constitutive")
            ame_dir_ind = os.path.join(args.outdir, f"ame_{r}_UFM1_independent_vs_Constitutive")
            tsv_dep = os.path.join(ame_dir_dep, "ame.tsv")
            tsv_ind = os.path.join(ame_dir_ind, "ame.tsv")
            
            if os.path.exists(tsv_dep) and os.path.exists(tsv_ind):
                plot_pdf = os.path.join(args.outdir, f"motif_enrichment_comparison_{r}.pdf")
                plot_script = os.path.join(args.script_dir, "plot_ame_comparison.R")
                cmd_plot = f"Rscript {plot_script} --dep={tsv_dep} --indep={tsv_ind} --out={plot_pdf}"
                run_cmd(cmd_plot, f"Plotting Consolidated Comparison for {r}")

def get_gc_len(fasta_path):
    lens = []
    gcs = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        s = str(record.seq).upper()
        if len(s) == 0: continue
        lens.append(len(s))
        gcs.append((s.count("G") + s.count("C")) / len(s) * 100.0)
    return lens, gcs

def compare_features(f1, f2, l1, l2, outdir, feature_prefix):
    len1, gc1 = get_gc_len(f1)
    len2, gc2 = get_gc_len(f2)
    if not len1 or not len2: return

    # Stats
    try:
        _, p_len = mannwhitneyu(len1, len2)
        _, p_gc = mannwhitneyu(gc1, gc2)
    except:
        p_len, p_gc = 1.0, 1.0

    # Plot Length
    df_len = pd.DataFrame([{"Group": l1, "Length": x} for x in len1] + [{"Group": l2, "Length": x} for x in len2])
    plt.figure(figsize=(6, 4))
    sns.boxplot(data=df_len, x="Group", y="Length")
    plt.title(f"{feature_prefix} Length\np={p_len:.2e}")
    plt.savefig(os.path.join(outdir, f"Comparison_{feature_prefix}_Length_{l1}_vs_{l2}.png"))
    plt.close()

    # Plot GC
    df_gc = pd.DataFrame([{"Group": l1, "GC": x} for x in gc1] + [{"Group": l2, "GC": x} for x in gc2])
    plt.figure(figsize=(6, 4))
    sns.boxplot(data=df_gc, x="Group", y="GC")
    plt.title(f"{feature_prefix} GC Content\np={p_gc:.2e}")
    plt.savefig(os.path.join(outdir, f"Comparison_{feature_prefix}_GC_{l1}_vs_{l2}.png"))
    plt.close()

def run_ame_and_plot(target_fa, control_fa, motif_db, l1, l2, region, outdir, script_dir):
    out_ame = os.path.join(outdir, f"ame_{region}_{l1}_vs_{l2}")
    if not os.path.exists(out_ame): # Only run if needed
        cmd_ame = f"mamba run -n meme_env ame --evalue-report-threshold 1000 --control {control_fa} --oc {out_ame} {target_fa} {motif_db}"
        run_cmd(cmd_ame, f"AME {region}: {l1} vs {l2}")

if __name__ == "__main__":
    main()
