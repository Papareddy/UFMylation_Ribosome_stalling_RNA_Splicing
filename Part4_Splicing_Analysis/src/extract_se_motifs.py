#!/usr/bin/env python3
"""extract_se_motifs.py
Extracts sequences for Skipped Exon (SE) events and performs motif/splice-signal analysis options.
Focuses on the Cassette Exon and its immediate splice sites.
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
    # rMATS SE specific columns
    # exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE
    raise KeyError(f"Column {base} not found in row keys: {list(row.keys())}")

def parse_se_coords(row):
    """
    Returns (chrom, strand, cass_start, cass_end, ups_ee, down_es)
    """
    chrom = get_col(row, "chr")
    if chrom.startswith("chr"): chrom = chrom.replace("chr", "")
    if chrom == "M": chrom = "MT"
    
    strand = get_col(row, "strand")
    
    # Cassette Exon
    cass_start = int(get_col(row, "exonStart_0base")) # 0-based
    cass_end = int(get_col(row, "exonEnd"))     # 1-based (or 0-based end? rMATS usually 0-based start, 1-based end for intervals)
    
    # Flanking Exons (we need them to define introns)
    ups_ee = int(get_col(row, "upstreamEE")) # End of upstream exon
    down_es = int(get_col(row, "downstreamES")) # Start of downstream exon
    
    return chrom, strand, cass_start, cass_end, ups_ee, down_es

def get_intervals(row):
    chrom, strand, c_start, c_end, u_ee, d_es = parse_se_coords(row)
    
    bed_entries = {}
    
    # Analyze the Cassette Exon Splice Sites (The ones that are "Skipped")
    # Upstream Intron: [u_ee, c_start]
    # Downstream Intron: [c_end, d_es]
    
    if strand == '+':
        # 3'SS (Acceptor) of Cassette Exon (at c_start)
        # Sequence: ...Intron | Exon...
        # [c_start - 20, c_start + 3]
        a_start = c_start - 20
        a_end = c_start + 3
        bed_entries['3ss'] = (chrom, a_start, a_end, strand)
        
        # Wide (AME)
        bed_entries['3ss_wide'] = (chrom, c_start - 50, c_start + 50, strand)
        
        # 5'SS (Donor) of Cassette Exon (at c_end)
        # Sequence: ...Exon | Intron...
        # [c_end - 3, c_end + 6]
        d_start = c_end - 3
        d_end = c_end + 6
        bed_entries['5ss'] = (chrom, d_start, d_end, strand)
        
        # Wide (AME)
        bed_entries['5ss_wide'] = (chrom, c_end - 50, c_end + 50, strand)
        
    else: # '-'
        # Transcript direction: High -> Low
        
        # 3'SS (Acceptor) of Cassette Exon (at c_end)
        # Sequence: ...Intron (High) | Exon (Low)... 
        # Junction at c_end.
        # Intron is Downstream in genomic coords (High coords).
        # Interval: [c_end - 3, c_end + 20] -> RevComp
        # genomic c_end is Exon Start (High).
        # genomic c_end+1 is Intron Start?
        # Actually in rMATS SE:
        # (+) u_ee < c_start < c_end < d_es
        # (-) u_ee < c_start < c_end < d_es (coords are always genomic increasing)
        # But 'upstream' logically is High coords.
        # rMATS 'upstreamEE' for (-) strand is actually the Downstream Exon in Transcription?
        # WARNING: rMATS variable names are genomic or logical?
        # Usually rMATS sorts by genomic coordinate:
        # upstream = lowest coord.
        # downstream = highest coord.
        # So for (-):
        # 'downstream' is logical Upstream (Donor).
        # 'upstream' is logical Downstream (Acceptor).
        # Wait, if strand is (-), Transcription goes d_es -> d_ee ... c_end -> c_start ... u_ee -> u_es
        
        # Let's verify standard rMATS behavior.
        # If I assume 'upstream' variable means 'genomic upstream' (low coord):
        # Then for (-):
        # Logical Upstream Exon is 'downstream' set.
        # Logical Intron 1 (upstream of SE) is between 'downstream' exon and 'cassette'.
        # That is interval [c_end, d_es].
        # Logical 3'SS (Acceptor) of SE is at c_end?
        # Yes.
        
        # MaxEntScan 3'SS (Acceptor) for (-):
        # We need [c_end - 3, c_end + 20]
        # RevComp converts:
        #   [c_end, c_end+20] (Intron) -> 5' part (Py tract)
        #   [c_end-3, c_end] (Exon) -> 3' part (TAG)
        # Correct.
        a_start = c_end - 3
        a_end = c_end + 20
        bed_entries['3ss'] = (chrom, a_start, a_end, strand)
        
        # Wide (AME)
        bed_entries['3ss_wide'] = (chrom, c_end - 50, c_end + 50, strand)
        
        # 5'SS (Donor) of Cassette Exon (at c_start)
        # Sequence: ...Exon | Intron (Low)...
        # Junction at c_start.
        # Interval: [c_start - 6, c_start + 3]
        # RevComp:
        #   [c_start, c_start+3] (Exon) -> 5' part (GGR)
        #   [c_start-6, c_start] (Intron) -> 3' part (AGT)
        # Correct.
        # 5'SS (Donor) of Cassette Exon (at c_start)
        # Sequence: ...Exon | Intron (Low)...
        # Interval: [c_start - 6, c_start + 3]
        d_start = c_start - 6
        d_end = c_start + 3
        bed_entries['5ss'] = (chrom, d_start, d_end, strand)
        
        # Wide (AME)
        bed_entries['5ss_wide'] = (chrom, c_start - 50, c_start + 50, strand)
        
        # Flanking Introns (150bp)
        # Upstream Intron (genomic higher than c_end): [c_end + 20, c_end + 170] -> RevComp (Intron downstream of cassette in transcription)
        # Downstream Intron (genomic lower than c_start): [c_start - 156, c_start - 6] -> RevComp (Intron upstream of cassette in transcription)
        
        # Warning: For (-) strand:
        # Transcript starts high, goes low.
        # Upstream Intron is GENOMIC HIGH (relative to exon).
        # Interval: [c_end + 6, c_end + 156]?
        #   c_end is Exon Start (High).
        #   Intron starts at c_end + 1.
        #   3'SS region is [c_end - 3, c_end + 20]. (Intron part is c_end to c_end+20).
        #   So Upstream Flank should be further into intron: [c_end + 20, c_end + 170].
        #   RevComp'd, this becomes the distal part of the upstream intron.
        # Downstream Intron is GENOMIC LOW.
        #   c_start is Exon End (Low).
        #   Intron starts at c_start - 1.
        #   5'SS region is [c_start - 6, c_start + 3]. (Intron part is c_start-6 to c_start).
        #   So Downstream Flank should be [c_start - 156, c_start - 6].
        
        bed_entries['intron_upstream'] = (chrom, c_end + 20, c_end + 170, strand)
        bed_entries['intron_downstream'] = (chrom, c_start - 156, c_start - 6, strand)

    # Cassette Exon Body
    bed_entries['exon'] = (chrom, c_start, c_end, strand)
    
    # For (+) strand, calculate flanks too
    if strand == '+':
        # Upstream Intron (Transcriptionally Upstream): [c_start - 170, c_start - 20]
        #   End of Upstream Intron is c_start.
        #   3'SS is [c_start - 20, c_start + 3].
        #   So flank is [c_start - 170, c_start - 20].
        bed_entries['intron_upstream'] = (chrom, c_start - 170, c_start - 20, strand)
        
        # Downstream Intron (Transcriptionally Downstream): [c_end + 6, c_end + 156]
        #   Start of Downstream Intron is c_end.
        #   5'SS is [c_end - 3, c_end + 6].
        #   So flank is [c_end + 6, c_end + 156].
        bed_entries['intron_downstream'] = (chrom, c_end + 6, c_end + 156, strand)
    
    return bed_entries

def write_beds(df, out_prefix):
    keys = ['5ss', '3ss', 'exon', 'intron_upstream', 'intron_downstream', '5ss_wide', '3ss_wide']
    beds = {k: open(f"{out_prefix}.{k}.bed", "w") for k in keys}
    
    for idx, row in df.iterrows():
        try:
            ivs = get_intervals(row)
            name = f"ID={idx}"
            for k, (c, s, e, strd) in ivs.items():
                if k in beds:
                    beds[k].write(f"{c}\t{s}\t{e}\t{name}\t0\t{strd}\n")
        except Exception as e:
            # Skip malformed lines
            continue
            
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
    
    # Filter for SE
    if "EventType" in lost.columns:
        lost = lost[lost["EventType"] == "SE"]
    if "EventType" in preserved.columns:
        preserved = preserved[preserved["EventType"] == "SE"]
    
    groups = {"UFM1_dependent": lost, "UFM1_independent": preserved}
    
    for grp_name, df in groups.items():
        print(f"[INFO] Processing {grp_name} (n={len(df)})...")
        if df.empty: continue
        prefix = os.path.join(args.outdir, grp_name)
        
        # 1. Write BEDs
        write_beds(df, prefix)
        
        # 2. Extract Sequences
        for kind in ['5ss', '3ss', 'exon', 'intron_upstream', 'intron_downstream']:
            bed = f"{prefix}.{kind}.bed"
            fa = f"{prefix}.{kind}.fa"
            # getfasta -s for strand specificity
            cmd = f"bedtools getfasta -s -fi {args.genome_fasta} -bed {bed} -fo {fa} -name"
            run_cmd(cmd, f"Extracting {kind} sequences for {grp_name}")
            
    # Note: MaxEnt scoring and AME will be orchestrated by a separate analyzer script
    # similar to analyze_ri_vs_constitutive.py, but updated for SE.

if __name__ == "__main__":
    main()
