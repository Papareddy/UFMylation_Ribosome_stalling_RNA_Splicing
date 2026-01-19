#!/usr/bin/env python3
"""extract_ri_motifs.py
Extracts sequences for Retained Intron (RI) events and performs motif/splice-signal analysis.
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
    # Try exact match, then .WT, then .UFM
    if base in row: return row[base]
    if f"{base}.WT" in row: return row[f"{base}.WT"]
    if f"{base}.UFM" in row: return row[f"{base}.UFM"]
    raise KeyError(f"Column {base} not found in row keys: {list(row.keys())}")

def parse_rmats_coords(row):
    """
    Returns (chrom, intron_start, intron_end, strand, ups_ee, down_es)
    RMATS RI: 
    upstreamEE = End of Upstream Exon
    downstreamES = Start of Downstream Exon
    Intron is enclosed between them.
    """
    chrom = get_col(row, "chr")
    
    # Normalize Chromosome names for Ensembl FASTA (1, 2, ..., MT)
    if chrom.startswith("chr"):
        chrom = chrom.replace("chr", "")
    if chrom == "M": chrom = "MT"
    
    strand = get_col(row, "strand")
    u_ee = int(get_col(row, "upstreamEE"))
    d_es = int(get_col(row, "downstreamES"))
    
    # Intron defined by rMATS
    return chrom, u_ee, d_es, strand, u_ee, d_es
    # Note: bedtools BED is 0-based, half-open [start, end).
    # rMATS is 0-based for starts, 1-based for ends usually?
    # Actually row["upstreamEE"] is 1-based end in rMATS usually.
    # Check header: upstreamEE.
    # If rMATS uses 0-based start, 1-based end:
    # Intron interval [upstreamEE, downstreamES).
    # Example: Exon [100, 200). Intron [200, 300). Exon [300, 400).
    # upstreamEE = 200. downstreamES = 300.
    # Correct.

def get_intervals(row):
    chrom, i_start, i_end, strand, u_ee, d_es = parse_rmats_coords(row)
    
    # Intervals to extract (0-based, half-open)
    # 5'SS and 3'SS for MaxEntScan
    # MaxEntScan 5'SS: 9mer (3 exon + 6 intron).
    # MaxEntScan 3'SS: 23mer (20 intron + 3 exon).
    # AME (Motif Analysis): Wider window (+/- 50bp)
    
    bed_entries = {}
    
    # Strand Logic
    if strand == '+':
        # 5'SS (Donor) at upstreamEE
        # Sequence: ...Exon [EE] Intron...
        # Want [EE-3, EE+6]
        d_start = u_ee - 3
        d_end = u_ee + 6
        
        # 3'SS (Acceptor) at downstreamES
        # Sequence: ...Intron [ES] Exon...
        # Want [ES-20, ES+3]
        a_start = d_es - 20
        a_end = d_es + 3
        
        bed_entries['5ss'] = (chrom, d_start, d_end, strand)
        bed_entries['3ss'] = (chrom, a_start, a_end, strand)
        
        # Wide (AME)
        bed_entries['5ss_wide'] = (chrom, u_ee - 50, u_ee + 50, strand)
        bed_entries['3ss_wide'] = (chrom, d_es - 50, d_es + 50, strand)
        
    else: # '-'
        # Transcript direction: High -> Low
        # Upstream Exon (Exon 1) is High Coord -> downstreamES logic?
        # User defined:
        # (-) Donor at downstreamES (Start of Downstream Exon, which is "Upstream" in transcript)
        # Sequence: ...Exon (High) | Intron (Low)...
        # Junction at downstreamES.
        # We want 3 bases of Exon (Start at d_es -> d_es+3) [RevComp'd becomes 5' end]
        # We want 6 bases of Intron (Start at d_es-6 -> d_es) [RevComp'd becomes 3' end]
        # Genomic Interval: [d_es - 6, d_es + 3]
        # RevComp([d_es-6, d_es+3]):
        #   [d_es, d_es+3] (Exon) becomes 5' part (GGR)
        #   [d_es-6, d_es] (Intron) becomes 3' part (AGT)
        # Correct.
        d_start = d_es - 6
        d_end = d_es + 3

        # (-) Acceptor at upstreamEE (End of Upstream Exon, which is "Downstream" in transcript)
        # Sequence: ...Intron | Exon...
        # Junction at upstreamEE.
        # Intron is High side? No, Exon 2 (Upstream Box, Low Coord) is "Downstream".
        # Intron is between them.
        # So Intron is High relative to Exon 2.
        # Junction at upstreamEE.
        # Genomic Interval: [upstreamEE - 3, upstreamEE + 20]
        # RevComp:
        #   [upstreamEE, upstreamEE+20] (Intron) -> 5' part (Py tract)
        #   [upstreamEE-3, upstreamEE] (Exon) -> 3' part (TAG/CAG)
        # Correct.
        a_start = u_ee - 3
        a_end = u_ee + 20
        
        bed_entries['5ss'] = (chrom, d_start, d_end, strand)
        bed_entries['3ss'] = (chrom, a_start, a_end, strand)

        # Wide (AME)
        # 5'SS (Donor) at downstreamES
        bed_entries['5ss_wide'] = (chrom, d_es - 50, d_es + 50, strand)
        # 3'SS (Acceptor) at upstreamEE
        bed_entries['3ss_wide'] = (chrom, u_ee - 50, u_ee + 50, strand)

    # Full Intron
    # [u_ee, d_es]
    bed_entries['intron'] = (chrom, u_ee, d_es, strand)
    
    return bed_entries

def write_beds(df, out_prefix):
    kinds = ['5ss', '3ss', 'intron', '5ss_wide', '3ss_wide']
    beds = {k: open(f"{out_prefix}.{k}.bed", "w") for k in kinds}
    
    for idx, row in df.iterrows():
        # Clean chromosome name if needed?
        # Assuming FASTA matches rMATS (e.g. chr1 vs 1).
        # Usually processed already.
        
        ivs = get_intervals(row)
        name = f"ID={idx}" # Use index or event ID
        
        for k, (c, s, e, strd) in ivs.items():
            # Formatting BED: chrom, start, end, name, score, strand
            beds[k].write(f"{c}\t{s}\t{e}\t{name}\t0\t{strd}\n")
            
    for f in beds.values(): f.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lost", required=True)
    parser.add_argument("--preserved", required=True)
    parser.add_argument("--genome_fasta", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # Load Data
    lost = pd.read_csv(args.lost, sep="\t")
    preserved = pd.read_csv(args.preserved, sep="\t")
    
    # Filter for RI
    if "EventType" in lost.columns:
        lost = lost[lost["EventType"] == "RI"]
    if "EventType" in preserved.columns:
        preserved = preserved[preserved["EventType"] == "RI"]
    
    # Process Groups
    groups = {"UFM1_dependent": lost, "UFM1_independent": preserved}
    
    for grp_name, df in groups.items():
        print(f"[INFO] Processing {grp_name} (n={len(df)})...")
        if df.empty: continue
        prefix = os.path.join(args.outdir, grp_name)
        
        # 1. Write BEDs
        write_beds(df, prefix)
        
        # 2. Extract Sequences
        for kind in ['5ss', '3ss', 'intron']:
            bed = f"{prefix}.{kind}.bed"
            fa = f"{prefix}.{kind}.fa"
            # getfasta -s for strand specificity
            cmd = f"bedtools getfasta -s -fi {args.genome_fasta} -bed {bed} -fo {fa} -name"
            run_cmd(cmd, f"Extracting {kind} sequences for {grp_name}")
            
        # 3. MaxEntScan Scores
        score_ss(f"{prefix}.5ss.fa", "5ss")
        score_ss(f"{prefix}.3ss.fa", "3ss")

    # 4. Motif Enrichment (AME)
    # Lost (Dependent) vs Preserved (Independent)
    if not lost.empty and not preserved.empty:
        print("[INFO] Running MEME/AME Enrichment (UFM1_dependent vs UFM1_independent)...")
        # Intron sequences
        lost_fa = os.path.join(args.outdir, "UFM1_dependent.intron.fa")
        pres_fa = os.path.join(args.outdir, "UFM1_independent.intron.fa")
        out_ame = os.path.join(args.outdir, "ame_results")
        
        # AME command
        # ame --control <control> <primary>
        # Check `ame` args.
        cmd = f"mamba run -n meme_env ame --control {pres_fa} --oc {out_ame} {lost_fa}" 
        # Note: AME uses a motif database.
        # If no database provided, AME can't score?
        # User said "De Novo Discovery ... DREME".
        # DREME finds motifs. AME enriches known motifs.
        # User said "Use MEME Suite (specifically AME) to find enriched motifs".
        # So I should run DREME.
        


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

if __name__ == "__main__":
    main()
