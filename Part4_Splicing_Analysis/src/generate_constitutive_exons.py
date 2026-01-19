#!/usr/bin/env python3
"""generate_constitutive_exons.py
Generates a background set of constitutive exons from a GTF file,
excluding those found in rMATS SE events (Lost/Preserved).
It produces BED files for:
- Exon
- 3'SS (Upstream of exon)
- 5'SS (Downstream of exon)
- Upstream Intron (Distal part)
- Downstream Intron (Proximal part)
"""

import os
import sys
import gzip
import random
import argparse
import pandas as pd
from collections import defaultdict

def parse_gtf_internal_exons(gtf_file):
    """
    Parses GTF to find INTERNAL exons (excluding first and last).
    Returns list of (chrom, start, end, strand, transcript_id)
    Start/End are 0-based, [start, end).
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
            
            # Filter for protein coding if possible (gene_type or transcript_type)
            # Simplified: just rely on structure for now, or filter by transcript_id presence in a list if provided.
            # Assuming 'retained_intron' or 'nonsense_mediated_decay' might be noisy, 
            # but usually 'constitutive' implies common to most isoforms. 
            # Here we just take internal exons of ANY transcript. 
            # The randomness helps smooth out noise.
            
            # Extract transcript_id
            attr = parts[8]
            tid = None
            if 'transcript_id "' in attr:
                tid = attr.split('transcript_id "')[1].split('"')[0]
            elif 'transcript_id' in attr:
                 tid = attr.split('transcript_id')[1].strip().split(';')[0].strip()
            
            if not tid: continue
            
            chrom = parts[0]
            if chrom.startswith("chr"): chrom = chrom.replace("chr", "")
            if chrom == "M": chrom = "MT"
            
            # GTF 1-based [start, end] -> BED 0-based [start-1, end)
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            
            transcripts[tid].append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand
            })
            
    # Identify Internal Exons
    internal_exons = []
    for tid, exons in transcripts.items():
        if len(exons) < 3: continue # Need at least 3 exons to have an internal one
        
        # Sort by coordinate
        exons.sort(key=lambda x: x["start"])
        
        # Identify internal indices based on strand
        # If (+): 0 is first, N is last. Internal: 1 to N-1.
        # If (-): 0 is last (genomic low), N is first (genomic high). Internal: 1 to N-1.
        # Regardless of strand, indices 1 to len-2 are "internal" physically.
        
        for i in range(1, len(exons) - 1):
            e = exons[i]
            # Need to capture flanking intron boundaries for context extraction?
            # Or just the exon coords and we infer flanks like extract_se_motifs.py?
            # extract_se_motifs uses cassette coords + offsets.
            # However, for constitutive exons, we need to know the intron boundaries 
            # to avoid going into the next exon if the intron is short.
            # But extract_se_motifs uses fixed offsets (e.g. +/- 150bp).
            # We will generate based on fixed offsets from the Exon boundaries to stay consistent.
            # We assume introns are long enough. If not, sequence will contain downstream exon, which is acceptable noise for background.
            
            internal_exons.append((e["chrom"], e["start"], e["end"], e["strand"], tid))
            
    print(f"[INFO] Found {len(internal_exons)} internal exons in GTF.")
    return internal_exons

def get_col(row, base):
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
    
    # 1. Load GTF Exons
    all_exons = parse_gtf_internal_exons(args.gtf)
    if not all_exons:
        print("[ERROR] No internal exons found.")
        sys.exit(1)
        
    # 2. Load Exclusion List (SE coordinates)
    # We exclude any exon that appears as a cassette exon in the rMATS inference.
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
                
                try:
                    # rMATS SE: exonStart_0base, exonEnd
                    estart = int(get_col(row, "exonStart_0base"))
                    eend = int(get_col(row, "exonEnd"))
                    if chrom:
                        excluded.add((chrom, estart, eend))
                except: pass
                
    print(f"[INFO] {len(excluded)} unique SE exons to exclude.")
    
    # 3. Filter
    constitutive = []
    for ex in all_exons:
        chrom, start, end, strand, tid = ex
        key = (chrom, start, end)
        if key not in excluded:
            constitutive.append(ex)
            
    print(f"[INFO] {len(constitutive)} constitutive exons remaining.")
    
    # 4. Sample
    if len(constitutive) > args.n_sample:
        print(f"[INFO] Sampling {args.n_sample} exons...")
        sampled = random.sample(constitutive, args.n_sample)
    else:
        sampled = constitutive
        
    # 5. Write BEDs
    print(f"[INFO] Writing to {args.out_bed}...")
    
    # Define outputs
    base_out = args.out_bed.replace(".bed", "")
    if base_out.endswith(".exon"): base_out = base_out.replace(".exon", "")
    
    files = {
        'exon': open(f"{base_out}.exon.bed", "w"),
        '3ss': open(f"{base_out}.3ss.bed", "w"),
        '5ss': open(f"{base_out}.5ss.bed", "w"),
        'intron_upstream': open(f"{base_out}.intron_upstream.bed", "w"),
        'intron_downstream': open(f"{base_out}.intron_downstream.bed", "w"),
        '5ss_wide': open(f"{base_out}.5ss_wide.bed", "w"),
        '3ss_wide': open(f"{base_out}.3ss_wide.bed", "w")
    }
    
    for idx, (chrom, c_start, c_end, strand, tid) in enumerate(sampled):
        name = f"ConstEx_{idx}_{tid}"
        
        # Mirroring extract_se_motifs.py logic for coordinates
        
        # Exon
        files['exon'].write(f"{chrom}\t{c_start}\t{c_end}\t{name}\t0\t{strand}\n")
        
        if strand == '+':
            # 3'SS (Acceptor): [start-20, start+3]
            files['3ss'].write(f"{chrom}\t{c_start-20}\t{c_start+3}\t{name}\t0\t{strand}\n")
            files['3ss_wide'].write(f"{chrom}\t{c_start-50}\t{c_start+50}\t{name}\t0\t{strand}\n")
            
            # 5'SS (Donor): [end-3, end+6]
            files['5ss'].write(f"{chrom}\t{c_end-3}\t{c_end+6}\t{name}\t0\t{strand}\n")
            files['5ss_wide'].write(f"{chrom}\t{c_end-50}\t{c_end+50}\t{name}\t0\t{strand}\n")
            
            # Upstream Intron (Distal): [end of upstream intron - 150, end of upstream intron]
            # End of upstream intron = c_start.
            # Range: [c_start-170, c_start-20] (Assuming -20 is splice site)
            files['intron_upstream'].write(f"{chrom}\t{c_start-170}\t{c_start-20}\t{name}\t0\t{strand}\n")
            
            # Downstream Intron (Proximal): [end, end+150]
            # Start of downstream intron = c_end.
            # Range: [c_end+6, c_end+156] (Assuming +6 is splice site)
            files['intron_downstream'].write(f"{chrom}\t{c_end+6}\t{c_end+156}\t{name}\t0\t{strand}\n")
            
        else: # (-) Strand
            # 3'SS (Acceptor) is at c_end (Genomic High)
            # [c_end-3, c_end+20]
            files['3ss'].write(f"{chrom}\t{c_end-3}\t{c_end+20}\t{name}\t0\t{strand}\n")
            files['3ss_wide'].write(f"{chrom}\t{c_end-50}\t{c_end+50}\t{name}\t0\t{strand}\n")
            
            # 5'SS (Donor) is at c_start (Genomic Low)
            # [c_start-6, c_start+3]
            files['5ss'].write(f"{chrom}\t{c_start-6}\t{c_start+3}\t{name}\t0\t{strand}\n")
            files['5ss_wide'].write(f"{chrom}\t{c_start-50}\t{c_start+50}\t{name}\t0\t{strand}\n")
            
            # Upstream Intron (Transcriptional Upstream = Genomic High)
            # Starts at c_end.
            # Range: [c_end+20, c_end+170]
            files['intron_upstream'].write(f"{chrom}\t{c_end+20}\t{c_end+170}\t{name}\t0\t{strand}\n")
            
            # Downstream Intron (Transcriptional Downstream = Genomic Low)
            # Ends at c_start.
            # Range: [c_start-156, c_start-6]
            files['intron_downstream'].write(f"{chrom}\t{c_start-156}\t{c_start-6}\t{name}\t0\t{strand}\n")
            
    for f in files.values(): f.close()

if __name__ == "__main__":
    main()
