#!/usr/bin/env python3
"""plot_rna_map.py
Generates Positional Motif Enrichment (RNA Map) plots for specific RBPs.
"""

import pandas as pd
import argparse
import os
import subprocess
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import motifs
from Bio.Seq import Seq
import math

# Logger setup
def log(msg):
    with open("debug_rna_map.log", "a") as f:
        f.write(msg + "\n")

def run_cmd(cmd, description):
    log(f"[EXEC] {description}...")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        log(f"[ERROR] Failed: {description}")
        sys.exit(1)

def get_col(row, base):
    if base in row: return row[base]
    if f"{base}.WT" in row: return row[f"{base}.WT"]
    if f"{base}.UFM" in row: return row[f"{base}.UFM"]
    return None

def extract_coords(df, window=150):
    bed_5ss = []
    bed_3ss = []
    
    for idx, row in df.iterrows():
        chrom = get_col(row, "chr")
        if chrom is None: continue
        if str(chrom).startswith("chr"): chrom = str(chrom).replace("chr", "")
        if chrom == "M": chrom = "MT"
        
        strand = get_col(row, "strand")
        u_ee = int(get_col(row, "upstreamEE"))
        d_es = int(get_col(row, "downstreamES"))
        
        name = f"ID={idx}"
        
        if strand == '+':
            p5_start = u_ee - window
            p5_end = u_ee + window
            bed_5ss.append((chrom, p5_start, p5_end, name, 0, strand))
            
            p3_start = d_es - window
            p3_end = d_es + window
            bed_3ss.append((chrom, p3_start, p3_end, name, 0, strand))
            
        else: # '-'
            p5_start = d_es - window
            p5_end = d_es + window
            bed_5ss.append((chrom, p5_start, p5_end, name, 0, strand))
            
            p3_start = u_ee - window
            p3_end = u_ee + window
            bed_3ss.append((chrom, p3_start, p3_end, name, 0, strand))

    return bed_5ss, bed_3ss

def make_bed(data, outfile):
    with open(outfile, "w") as f:
        for x in data:
            f.write(f"{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]}\t{x[5]}\n")

class CustomMotif:
    def __init__(self, name, alt, matrix):
        self.name = name
        self.alt = alt
        self.pwm = matrix # List of dictionaries/arrays
        self.length = len(matrix)
        self.pssm = []
        bg = 0.25
        for row in self.pwm:
            pssm_row = {}
            for base, p in zip("ACGU", row):
                log_odds = math.log2((p + 1e-6) / bg)
                pssm_row[base] = log_odds
            self.pssm.append(pssm_row)

def parse_meme_manual(filepath):
    motifs_list = []
    current_name = None
    current_alt = None
    matrix = []
    in_matrix = False
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                if in_matrix and matrix:
                    motifs_list.append(CustomMotif(current_name, current_alt, matrix))
                    matrix = []
                    in_matrix = False
                continue
            
            if line.startswith("MOTIF"):
                if in_matrix and matrix:
                    motifs_list.append(CustomMotif(current_name, current_alt, matrix))
                    matrix = []
                    in_matrix = False
                
                parts = line.split()
                if len(parts) >= 3:
                    current_name = parts[1]
                    current_alt = parts[2]
                elif len(parts) == 2:
                    current_name = parts[1]
                    current_alt = ""
                else:
                    current_name = "Unknown"
            
            elif line.startswith("letter-probability matrix:"):
                in_matrix = True
                continue
                
            elif in_matrix:
                try:
                    probs = [float(x) for x in line.split()]
                    if len(probs) == 4:
                        matrix.append(probs)
                except ValueError:
                    pass
        
        if in_matrix and matrix:
            motifs_list.append(CustomMotif(current_name, current_alt, matrix))
    return motifs_list

def scan_sequences(fasta_file, motifs_list):
    log(f"[DEBUG] Scanning {fasta_file} for {len(motifs_list)} motifs...")
    
    seq_len = 0
    with open(fasta_file) as f:
        for line in f:
            if not line.startswith(">"):
                seq_len = len(line.strip())
                break
    if seq_len == 0: 
        log(f"[DEBUG] Empty or invalid FASTA: {fasta_file}")
        return None, 0
    
    total_counts = {m.name: np.zeros(seq_len) for m in motifs_list}
    n_seqs = 0
    hits_total = 0
    max_score_obs = -999.0
    
    from Bio import SeqIO
    for record in SeqIO.parse(fasta_file, "fasta"):
        n_seqs += 1
        seq = str(record.seq).upper().replace("T", "U")
        if len(seq) != seq_len: continue
        
        for m in motifs_list:
            w = m.length
            for i in range(seq_len - w + 1):
                subseq = seq[i : i+w]
                score = 0
                valid = True
                for j, base in enumerate(subseq):
                    if base not in "ACGU": 
                        valid = False; break
                    score += m.pssm[j][base]
                
                max_score_obs = max(max_score_obs, score)
                # Lower threshold to ensure hits? 2.0 is log2(4) = 4x.
                # Max score is usually > 10.
                if valid and score > 2.0: 
                     center = i + w // 2
                     if center < seq_len:
                         total_counts[m.name][center] += 1
                         hits_total += 1

    log(f"[DEBUG] Scanned {n_seqs} sequences. Max Score: {max_score_obs:.2f}. Total Hits: {hits_total}")
    
    if n_seqs > 0:
        density_maps = {k: v / n_seqs for k, v in total_counts.items()}
        return density_maps, seq_len
    return None, 0

def main():
    if os.path.exists("debug_rna_map.log"): os.remove("debug_rna_map.log")
    log("Starting RNA Map analysis...")
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--lost", required=True)
    parser.add_argument("--preserved", required=True)
    parser.add_argument("--constitutive", required=False)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--motif_db", required=True)
    parser.add_argument("--motifs", nargs="+", help="Motif names/IDs to plot")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--window", type=int, default=150, help="Extraction window (+/-)")
    parser.add_argument("--smooth", type=int, default=40, help="Smoothing window size (default: 40)")
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    log(f"Arguments: {args}")
    
    # 1. Load Data
    dfs = {}
    if args.lost: 
        log(f"Loading Lost (Dependent): {args.lost}")
        dfs["UFM1_dependent"] = pd.read_csv(args.lost, sep="\t")
    if args.preserved: 
        log(f"Loading Preserved (Independent): {args.preserved}")
        dfs["UFM1_independent"] = pd.read_csv(args.preserved, sep="\t")
    
    # 2. Extract BEDs
    beds_5ss = {}
    beds_3ss = {}
    
    for name, df in dfs.items():
        if "EventType" in df.columns: 
            log(f"Filtering {name} for RI events...")
            df = df[df["EventType"] == "RI"]
        
        log(f"Extracting coords for {name} (n={len(df)})...")
        b5, b3 = extract_coords(df, window=args.window)
        
        bed5_path = os.path.join(args.outdir, f"{name}.5ss.bed")
        bed3_path = os.path.join(args.outdir, f"{name}.3ss.bed")
        
        make_bed(b5, bed5_path)
        make_bed(b3, bed3_path)
        
        beds_5ss[name] = bed5_path
        beds_3ss[name] = bed3_path
        
    if args.constitutive:
        log(f"Loading Constitutive: {args.constitutive}")
        c_bed = []
        with open(args.constitutive) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 6: continue
                c_bed.append(parts)
        
        c5 = []; c3 = []
        for p in c_bed:
            chrom, s, e, nm, sc, strd = p[0], int(p[1]), int(p[2]), p[3], p[4], p[5]
            if strd == '+':
                c5.append([chrom, s-args.window, s+args.window, nm, 0, strd])
                c3.append([chrom, e-args.window, e+args.window, nm, 0, strd])
            else:
                c5.append([chrom, e-args.window, e+args.window, nm, 0, strd])
                c3.append([chrom, s-args.window, s+args.window, nm, 0, strd])
                
        make_bed(c5, os.path.join(args.outdir, "Constitutive.5ss.bed"))
        make_bed(c3, os.path.join(args.outdir, "Constitutive.3ss.bed"))
        beds_5ss["Constitutive"] = os.path.join(args.outdir, "Constitutive.5ss.bed")
        beds_3ss["Constitutive"] = os.path.join(args.outdir, "Constitutive.3ss.bed")

    # 3. Get FASTA
    fastas_5ss = {}
    fastas_3ss = {}
    
    for grp, bed in beds_5ss.items():
        fa = bed.replace(".bed", ".fa")
        if os.path.exists(fa): os.remove(fa)
        cmd = f"bedtools getfasta -s -fi {args.genome} -bed {bed} -fo {fa} -name"
        run_cmd(cmd, f"Extracting 5'SS sequences for {grp}")
        fastas_5ss[grp] = fa
        
    for grp, bed in beds_3ss.items():
        fa = bed.replace(".bed", ".fa")
        if os.path.exists(fa): os.remove(fa)
        cmd = f"bedtools getfasta -s -fi {args.genome} -bed {bed} -fo {fa} -name"
        run_cmd(cmd, f"Extracting 3'SS sequences for {grp}")
        fastas_3ss[grp] = fa

    # 4. Parse DB
    log(f"Loading motifs {args.motifs} from {args.motif_db}...")
    all_motifs = parse_meme_manual(args.motif_db)
    
    target_motifs = []
    for m in all_motifs:
        ids = [m.name, m.get("alt", "") if hasattr(m, "get") else m.alt] 
        
        found = False
        for req in args.motifs:
            req_u = req.upper()
            if (req_u == m.name.upper()) or (req_u == m.alt.upper()):
                found = True
                m.display_name = req 
                break
        if found:
            target_motifs.append(m)
            
    log(f"Found {len(target_motifs)} matching motifs: {[m.name for m in target_motifs]}")


    # 5. Scan and Plot
    os.makedirs(args.outdir, exist_ok=True)
    sns.set_style("whitegrid")
    
    for motif in target_motifs:
        dname = getattr(motif, "display_name", motif.name)
        log(f"Processing Motif: {motif.name} ({dname})...")
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Panel 1: 5'SS
        max_y = 0
        for grp, fa in fastas_5ss.items():
            density, length = scan_sequences(fa, [motif])
            if density:
                y = density[motif.name]
                s = pd.Series(y)
                # Use args.smooth for window
                y_smooth = s.rolling(window=args.smooth, center=True, min_periods=1).mean()
                
                x = np.arange(-length//2, length//2)
                if len(x) != len(y_smooth):
                    log(f"Warning: dimension mismatch {len(x)} vs {len(y_smooth)}")
                    continue
                    
                axes[0].plot(x, y_smooth, label=grp, linewidth=2)
                if not y_smooth.dropna().empty:
                    max_y = max(max_y, y_smooth.max())
        
        axes[0].set_title(f"{dname} ({motif.name}) - 5' Splice Site", fontsize=14)
        axes[0].set_xlabel("Dist. to 5'SS (bp)", fontsize=12)
        axes[0].set_ylabel("Motif Density", fontsize=12)
        axes[0].legend()
        axes[0].axvline(0, color='k', linestyle='--', alpha=0.5)
        # Add Exon/Intron labels
        axes[0].text(-args.window/2, max_y*1.05, "Exon", ha='center', fontweight='bold')
        axes[0].text(args.window/2, max_y*1.05, "Intron", ha='center', fontweight='bold')
        
        # Panel 2: 3'SS
        for grp, fa in fastas_3ss.items():
            density, length = scan_sequences(fa, [motif])
            if density:
                y = density[motif.name]
                s = pd.Series(y)
                # Use args.smooth for window
                y_smooth = s.rolling(window=args.smooth, center=True, min_periods=1).mean()
                x = np.arange(-length//2, length//2)
                
                axes[1].plot(x, y_smooth, label=grp, linewidth=2)
                if not y_smooth.dropna().empty:
                    max_y = max(max_y, y_smooth.max())
        
        axes[1].set_title(f"{dname} ({motif.name}) - 3' Splice Site", fontsize=14)
        axes[1].set_xlabel("Dist. to 3'SS (bp)", fontsize=12)
        axes[1].axvline(0, color='k', linestyle='--', alpha=0.5)
        axes[1].text(-args.window/2, max_y*1.05, "Intron", ha='center', fontweight='bold')
        axes[1].text(args.window/2, max_y*1.05, "Exon", ha='center', fontweight='bold')
        
        plt.suptitle(f"Positional Enrichment: {dname}", fontsize=16)
        plt.tight_layout()
        # Unique filename using ID
        out_plot = os.path.join(args.outdir, f"RNA_Map_{dname}_{motif.name}.png")
        plt.savefig(out_plot, dpi=300)
        plt.close()
        log(f"Saved {out_plot}")

if __name__ == "__main__":
    main()
