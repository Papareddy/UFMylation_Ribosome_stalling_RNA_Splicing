#!/usr/bin/env python3
"""generate_constitutive_introns.py
Generates a background set of constitutive introns from a GTF file,
excluding those found in rMATS RI events (Lost/Preserved).
"""

import os
import sys
import gzip
import random
import argparse
import pandas as pd
from collections import defaultdict

def parse_gtf_introns(gtf_file):
    """
    Parses GTF to find introns.
    Intron = Gap between Exons of the same transcript.
    Returns list of (chrom, start, end, strand, transcript_id)
    Note: GTF is 1-based.
    Bedtools/Internal logic: Use 0-based half-open [start, end).
    GTF Exon: [start, end] (inclusive 1-based).
    Intron: (exon1_end + 1, exon2_start - 1) inclusive 1-based?
    Wait. 0-based:
    Exon 1: [0, 100). (GTF: 1, 100).
    Exon 2: [200, 300). (GTF: 201, 300).
    Intron: [100, 200). (GTF: 101, 200).
    So Intron Start (0-based) = Exon 1 End (0-based).
    Intron End (0-based) = Exon 2 Start (0-based).
    Correct.
    """
    print(f"[INFO] Parsing GTF: {gtf_file}")
    transcripts = defaultdict(list)
    
    open_func = gzip.open if gtf_file.endswith(".gz") else open
    with open_func(gtf_file, "rt") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9: continue
            if parts[2] != "exon": continue
            
            # Extract transcript_id
            attr = parts[8]
            tid = None
            if 'transcript_id "' in attr:
                tid = attr.split('transcript_id "')[1].split('"')[0]
            elif 'transcript_id' in attr: # simplified
                 tid = attr.split('transcript_id')[1].strip().split(';')[0].strip()
            
            if not tid: continue
            
            chrom = parts[0]
            if chrom.startswith("chr"): chrom = chrom.replace("chr", "") # Normalize
            if chrom == "M": chrom = "MT"
            
            start = int(parts[3]) - 1 # 0-based
            end = int(parts[4])       # 0-based exclusive
            strand = parts[6]
            
            transcripts[tid].append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand
            })
            
    # Infer Introns
    introns = []
    for tid, exons in transcripts.items():
        if len(exons) < 2: continue
        # Sort by genomic coordinate
        exons.sort(key=lambda x: x["start"])
        
        # Iterate gaps
        # Exon i and i+1
        for i in range(len(exons) - 1):
            e1 = exons[i]
            e2 = exons[i+1]
            
            # Intron
            chrom = e1["chrom"]
            strand = e1["strand"]
            istart = e1["end"] # 0-based
            iend = e2["start"] # 0-based
            
            if iend > istart:
                introns.append((chrom, istart, iend, strand, tid))
                
    print(f"[INFO] Found {len(introns)} introns in GTF.")
    return introns

def get_col(row, base):
    # Same helper as extract_ri_motifs
    if base in row: return row[base]
    if f"{base}.WT" in row: return row[f"{base}.WT"]
    if f"{base}.UFM" in row: return row[f"{base}.UFM"]
    return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--exclude_files", nargs="+", help="rMATS TSV files to exclude")
    parser.add_argument("--n_sample", type=int, default=5000)
    parser.add_argument("--out_bed", required=True)
    args = parser.parse_args()
    
    # 1. Load GTF Introns
    all_introns = parse_gtf_introns(args.gtf)
    if not all_introns:
        print("[ERROR] No introns found.")
        sys.exit(1)
        
    # 2. Load Exclusion List (RI coordinates)
    excluded = set()
    if args.exclude_files:
        for fpath in args.exclude_files:
            if not os.path.exists(fpath): continue
            print(f"[INFO] Loading exclusions from {fpath}...")
            df = pd.read_csv(fpath, sep="\t")
            for _, row in df.iterrows():
                chrom = get_col(row, "chr")
                if chrom and chrom.startswith("chr"): chrom = chrom.replace("chr", "")
                if chrom == "M": chrom = "MT"
                
                # rMATS RI: Intron is [upstreamEE, downstreamES]
                # Assuming 0-based start, 0-based exclusive end (BED)?
                # Wait. upstreamEE (in rMATS) usually IS the start of the intron (0-based).
                # downstreamES is the end (0-based exclusive).
                try:
                    u_ee = int(get_col(row, "upstreamEE"))
                    d_es = int(get_col(row, "downstreamES"))
                    if chrom:
                        excluded.add((chrom, u_ee, d_es))
                except: pass
                
    print(f"[INFO] {len(excluded)} unique RI events to exclude.")
    
    # 3. Filter
    constitutive = []
    # Set lookup is faster
    # Keys: (chrom, start, end)
    # Are strand important for uniqueness? Intron is physical. 
    # Usually (chrom, start, end) defines the intron location.
    
    for iv in all_introns:
        chrom, start, end, strand, tid = iv
        key = (chrom, start, end)
        if key not in excluded:
            constitutive.append(iv)
            
    print(f"[INFO] {len(constitutive)} constitutive introns remaining.")
    
    # 4. Sample
    if len(constitutive) > args.n_sample:
        print(f"[INFO] Sampling {args.n_sample} introns...")
        sampled = random.sample(constitutive, args.n_sample)
    else:
        sampled = constitutive
        
    # 5. Write BED
    # BED: chrom, start, end, name, score, strand
    print(f"[INFO] Writing to {args.out_bed}...")
    
    # Define outputs
    base_out = args.out_bed.replace(".bed", "")
    if base_out.endswith(".intron"): base_out = base_out.replace(".intron", "")
    
    bed_intron = open(f"{base_out}.intron.bed", "w")
    bed_5ss = open(f"{base_out}.5ss.bed", "w")
    bed_3ss = open(f"{base_out}.3ss.bed", "w")
    
    for idx, (chrom, start, end, strand, tid) in enumerate(sampled):
        name = f"Const_{idx}_{tid}"
        
        # Intron
        bed_intron.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
        
        # Splice Sites
        if strand == "+":
            # 5'SS: [start-3, start+6]
            d5_s, d5_e = start - 3, start + 6
            # 3'SS: [end-20, end+3]
            d3_s, d3_e = end - 20, end + 3
        else:
            # 5'SS (Donor) at 'end': [end-6, end+3]
            d5_s, d5_e = end - 6, end + 3
            # 3'SS (Acceptor) at 'start': [start-3, start+20]
            d3_s, d3_e = start - 3, start + 20
            
        bed_5ss.write(f"{chrom}\t{d5_s}\t{d5_e}\t{name}\t0\t{strand}\n")
        bed_3ss.write(f"{chrom}\t{d3_s}\t{d3_e}\t{name}\t0\t{strand}\n")
        
    bed_intron.close()
    bed_5ss.close()
    bed_3ss.close()

if __name__ == "__main__":
    main()
