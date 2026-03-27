#!/usr/bin/env python3
"""
run_ese_analysis.py
Comprehensive ESE/RBP Motif Analysis Pipeline (Step 10B)

Runs the following analyses:
1. Sequence extraction (5'SS ±50bp, 3'SS ±50bp, full introns)
2. CAG and ESE consensus scanning
3. FIMO RBP scanning with CIS-BP database
4. AME statistical enrichment analysis
5. Generate comparison plots
"""

import os
import sys
import argparse
import subprocess
import shutil


def run_cmd(cmd, cwd=None):
    """Execute command and check for errors."""
    print(f"[EXEC] {cmd[:100]}...")
    result = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[WARN] Command returned {result.returncode}")
        if result.stderr:
            print(f"[STDERR] {result.stderr[:500]}")
    return result.returncode == 0


def main():
    parser = argparse.ArgumentParser(description="ESE/RBP Motif Analysis Pipeline")
    parser.add_argument("--events_rds", required=True, help="Path to UFM1_events_rich.rds")
    parser.add_argument("--genome_fasta", required=True, help="Path to genome FASTA")
    parser.add_argument("--gtf", required=True, help="Path to GTF annotation")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--motif_db", default="data/motifs/CisBP_Human_All.meme", help="CIS-BP motif database")
    parser.add_argument("--meme_env", default="meme_env", help="Conda environment for MEME suite")
    args = parser.parse_args()
    
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    
    src_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(src_dir)
    
    print("="*70)
    print("ESE/RBP MOTIF ANALYSIS PIPELINE")
    print("="*70)
    
    # === STEP 1: Sequence Extraction ===
    print("\n[STEP 1] Extracting ROI sequences...")
    extraction_script = os.path.join(src_dir, "ese_motif_extraction.R")
    run_cmd(f"Rscript {extraction_script} {args.events_rds} {args.genome_fasta} {args.gtf} {outdir}")
    
    # === STEP 2: Full Intron Extraction ===
    print("\n[STEP 2] Extracting full intron sequences...")
    full_intron_dir = os.path.join(outdir, "full_introns")
    os.makedirs(full_intron_dir, exist_ok=True)
    
    # Use R to extract full introns
    r_extract = f"""
    library(GenomicRanges)
    events <- readRDS('{args.events_rds}')
    ri <- events[events$EventType == 'RI']
    dep <- ri[ri$Group == 'UFM1_dependent']
    ind <- ri[ri$Group == 'UFM1_independent']
    
    write_bed <- function(gr, fp) {{
      chrom <- sub('^chr', '', as.character(seqnames(gr)))
      chrom <- ifelse(chrom == 'M', 'MT', chrom)
      bed_df <- data.frame(chrom, start(gr)-1, end(gr), paste0('intron_', seq_along(gr)), 0, as.character(strand(gr)))
      write.table(bed_df, fp, sep='\\t', quote=F, row.names=F, col.names=F)
    }}
    write_bed(dep, '{full_intron_dir}/UFM1_dependent_full.bed')
    write_bed(ind, '{full_intron_dir}/UFM1_independent_full.bed')
    """
    run_cmd(f'Rscript -e "{r_extract}"')
    run_cmd(f"bedtools getfasta -s -fi {args.genome_fasta} -bed {full_intron_dir}/UFM1_dependent_full.bed -fo {full_intron_dir}/UFM1_dependent_full.fa -name")
    run_cmd(f"bedtools getfasta -s -fi {args.genome_fasta} -bed {full_intron_dir}/UFM1_independent_full.bed -fo {full_intron_dir}/UFM1_independent_full.fa -name")
    
    # === STEP 3: ESE Consensus Scanning ===
    print("\n[STEP 3] Running ESE consensus scanning...")
    consensus_dir = os.path.join(outdir, "consensus_scan")
    consensus_script = os.path.join(src_dir, "ese_consensus_scanning.py")
    run_cmd(f"python3 {consensus_script} --extraction_dir {outdir} --outdir {consensus_dir}")
    
    # === STEP 4: FIMO RBP Scanning ===
    print("\n[STEP 4] Running FIMO RBP scanning...")
    fimo_dir = os.path.join(outdir, "fimo_results")
    os.makedirs(fimo_dir, exist_ok=True)
    
    for group in ['UFM1_dependent', 'UFM1_independent', 'Constitutive']:
        for roi in ['roi1_5ss', 'roi2_3ss']:
            fa_file = os.path.join(outdir, f"{group}.{roi}.fa")
            if os.path.exists(fa_file):
                fimo_out = os.path.join(fimo_dir, f"{group}_{roi.replace('.', '_')}")
                run_cmd(f"mamba run -n {args.meme_env} fimo --oc {fimo_out} --thresh 1e-3 --verbosity 1 {args.motif_db} {fa_file}")
    
    # === STEP 5: AME Statistical Analysis ===
    print("\n[STEP 5] Running AME enrichment analysis...")
    ame_dir = os.path.join(outdir, "ame_results")
    os.makedirs(ame_dir, exist_ok=True)
    
    # 5'SS comparisons
    for comparison, test, ctrl in [
        ("dep_vs_ctrl_5ss", "UFM1_dependent.roi1_5ss.fa", "Constitutive.roi1_5ss.fa"),
        ("ind_vs_ctrl_5ss", "UFM1_independent.roi1_5ss.fa", "Constitutive.roi1_5ss.fa"),
        ("dep_vs_ctrl_3ss", "UFM1_dependent.roi2_3ss.fa", "Constitutive.roi2_3ss.fa"),
        ("ind_vs_ctrl_3ss", "UFM1_independent.roi2_3ss.fa", "Constitutive.roi2_3ss.fa"),
    ]:
        test_path = os.path.join(outdir, test)
        ctrl_path = os.path.join(outdir, ctrl)
        if os.path.exists(test_path) and os.path.exists(ctrl_path):
            ame_out = os.path.join(ame_dir, comparison)
            run_cmd(f"mamba run -n {args.meme_env} ame --control {ctrl_path} --oc {ame_out} {test_path} {args.motif_db}")
    
    # Full intron AME
    run_cmd(f"mamba run -n {args.meme_env} ame --control {full_intron_dir}/UFM1_independent_full.fa --oc {ame_dir}/full_intron_dep_vs_ind {full_intron_dir}/UFM1_dependent_full.fa {args.motif_db}")
    
    # === STEP 6: Generate Plots ===
    print("\n[STEP 6] Generating visualization plots...")
    plots_dir = os.path.join(outdir, "plots")
    plot_script = os.path.join(src_dir, "plot_ese_enhanced.R")
    run_cmd(f"Rscript {plot_script} {outdir} {plots_dir}")
    
    # Generate FIMO-based scatter plots
    cag_plot_script = os.path.join(src_dir, "plot_cag_positional.py")
    if os.path.exists(cag_plot_script):
        run_cmd(f"python3 {cag_plot_script}")
    
    print("\n" + "="*70)
    print(f"ESE ANALYSIS COMPLETE")
    print(f"Results in: {outdir}")
    print("="*70)
    
    # Summary
    print("\nOutput files:")
    for subdir in ['consensus_scan', 'fimo_results', 'ame_results', 'plots', 'full_introns']:
        full_path = os.path.join(outdir, subdir)
        if os.path.exists(full_path):
            n_files = len([f for f in os.listdir(full_path) if not f.startswith('.')])
            print(f"  {subdir}/: {n_files} files")


if __name__ == "__main__":
    main()
