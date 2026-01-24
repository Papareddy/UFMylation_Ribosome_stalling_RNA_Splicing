#!/usr/bin/env python3
"""extract_3ss_sequences.py
Extracts +/- 20bp sequences around the 3' splice site (3'SS) for UFM1-dependent,
UFM1-independent, and constitutive introns.
"""

import pandas as pd
import argparse
import os
import subprocess
import sys

def run_cmd(cmd, description):
    print(f"[EXEC] {description}...")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed: {description}")
        sys.exit(1)

def get_col(row, base):
    if base in row: return row[base]
    if f"{base}.WT" in row: return row[f"{base}.WT"]
    if f"{base}.UFM" in row: return row[f"{base}.UFM"]
    raise KeyError(f"Column {base} not found")

def get_3ss_coords(row):
    chrom = get_col(row, "chr")
    if chrom.startswith("chr"): chrom = chrom.replace("chr", "")
    if chrom == "M": chrom = "MT"
    
    strand = get_col(row, "strand")
    u_ee = int(get_col(row, "upstreamEE"))
    d_es = int(get_col(row, "downstreamES"))
    
    # + strand: 3'SS is at downstreamES
    # - strand: 3'SS is at upstreamEE
    if strand == '+':
        pos = d_es
    else:
        pos = u_ee
    
    # Window +/- 20bp (41bp total: -20 to +20 relative to junction)
    return chrom, pos - 20, pos + 20, strand

def write_bed(df, out_path):
    with open(out_path, "w") as f:
        for idx, row in df.iterrows():
            c, s, e, strd = get_3ss_coords(row)
            f.write(f"{c}\t{s}\t{e}\tID={idx}\t0\t{strd}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lost", required=True)
    parser.add_argument("--preserved", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--genome_fasta", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--script_dir", required=True)
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # 1. Process Groups
    groups = {
        "UFM1_dependent": pd.read_csv(args.lost, sep="\t"),
        "UFM1_independent": pd.read_csv(args.preserved, sep="\t")
    }
    
    for name, df in groups.items():
        if "EventType" in df.columns:
            df = df[df["EventType"] == "RI"]
        
        bed = os.path.join(args.outdir, f"{name}.3ss_20bp.bed")
        fa = os.path.join(args.outdir, f"{name}.3ss_20bp.fa")
        
        write_bed(df, bed)
        run_cmd(f"bedtools getfasta -s -fi {args.genome_fasta} -bed {bed} -fo {fa} -name", f"Extracting {name} 3'SS sequences")

    # 2. Constitutive Background
    const_bed_base = os.path.join(args.outdir, "constitutive")
    gen_script = os.path.join(args.script_dir, "generate_constitutive_introns.py")
    
    # Run helper to get BEDs
    cmd_gen = f"python3 {gen_script} --gtf {args.gtf} --out_bed {const_bed_base}.bed --exclude_files {args.lost} {args.preserved} --n_sample 5000"
    run_cmd(cmd_gen, "Generating Constitutive Intron Background")
    
    # The generator produces .3ss.bed (which is MaxEnt sized)
    # We need to re-read that or trust its coordinates and expand/shrink.
    # Actually, let's just use the .3ss.bed it produced and remap to +/- 20bp.
    # rMATS 3'SS logic in generate_constitutive_introns.py:
    # strand +: [ES-20, ES+3]
    # strand -: [EE-3, EE+20]
    
    # It's easier to just parse the raw constitutive_introns.bed (full intron) and apply same logic.
    full_intron_bed = f"{const_bed_base}.intron.bed"
    if os.path.exists(full_intron_bed):
        const_df = pd.read_csv(full_intron_bed, sep="\t", header=None, names=["chr", "start", "end", "name", "score", "strand"])
        # In this BED, start/end ARE the intron.
        # strand +: 3'SS at 'end' (downstreamES)
        # strand -: 3'SS at 'start' (upstreamEE)
        
        bed_20 = os.path.join(args.outdir, "Constitutive.3ss_20bp.bed")
        fa_20 = os.path.join(args.outdir, "Constitutive.3ss_20bp.fa")
        
        with open(bed_20, "w") as f:
            for idx, row in const_df.iterrows():
                chrom = str(row["chr"])
                strand = row["strand"]
                pos = row["end"] if strand == "+" else row["start"]
                f.write(f"{chrom}\t{pos-20}\t{pos+20}\tID={idx}\t0\t{strand}\n")
        
        run_cmd(f"bedtools getfasta -s -fi {args.genome_fasta} -bed {bed_20} -fo {fa_20} -name", "Extracting Constitutive 3'SS sequences")

if __name__ == "__main__":
    main()
