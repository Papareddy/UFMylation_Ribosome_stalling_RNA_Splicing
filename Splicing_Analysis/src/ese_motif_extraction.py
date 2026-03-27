#!/usr/bin/env python3
"""
ese_motif_extraction.py
Step A: Sequence Extraction & FASTA Generation for ESE/RBP motif analysis.

Extracts ROI sequences around 5'SS, 3'SS, and Branch Point regions for:
- UFM1-Dependent Introns
- UFM1-Independent Introns
- Constitutive Introns (Control)
"""

import pandas as pd
import argparse
import os
import subprocess
import sys
from pybedtools import BedTool
import pybedtools

def run_cmd(cmd, description):
    """Execute a shell command."""
    print(f"[EXEC] {description}...")
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed: {description}")
        sys.exit(1)


def get_col(row, base):
    """Flexible column getter for rMATS/GRanges-style DataFrames."""
    # Handle common column name variations
    aliases = {
        'chr': ['chr', 'seqnames', 'chrom', 'chromosome'],
        'upstreamEE': ['upstreamEE', 'riExonStart_0base'],
        'downstreamES': ['downstreamES', 'riExonEnd'],
    }
    
    # If base has known aliases, try them
    if base in aliases:
        for alias in aliases[base]:
            if alias in row.index:
                return row[alias]
    
    # Direct match
    if base in row.index: return row[base]
    if f"{base}.WT" in row.index: return row[f"{base}.WT"]
    if f"{base}.UFM" in row.index: return row[f"{base}.UFM"]
    
    raise KeyError(f"Column {base} not found in row keys: {list(row.index)}")



def extract_roi_beds(df, out_prefix, group_name):
    """
    Extract ROI BED files for a given group of RI events.
    
    ROI_1: 5'SS ±30bp (centered on Exon|Intron junction)
    ROI_2: 3'SS ±30bp (centered on Intron|Exon junction)
    ROI_3: Branch Point -15 to -45bp upstream of 3'SS
    """
    beds = {
        'roi1_5ss': [],      # ±30bp around 5'SS
        'roi2_3ss': [],      # ±30bp around 3'SS
        'roi3_bp': [],       # Branch point region
    }
    
    for idx, row in df.iterrows():
        chrom = get_col(row, "chr")
        if chrom.startswith("chr"):
            chrom = chrom.replace("chr", "")
        if chrom == "M": chrom = "MT"
        
        strand = get_col(row, "strand")
        u_ee = int(get_col(row, "upstreamEE"))  # Donor junction
        d_es = int(get_col(row, "downstreamES"))  # Acceptor junction
        
        name = f"{group_name}_{idx}"
        
        if strand == '+':
            # 5'SS (Donor) at upstreamEE: Exon|Intron junction
            # ROI_1: ±30bp -> [u_ee - 30, u_ee + 30]
            roi1_start, roi1_end = u_ee - 30, u_ee + 30
            
            # 3'SS (Acceptor) at downstreamES: Intron|Exon junction
            # ROI_2: ±30bp -> [d_es - 30, d_es + 30]
            roi2_start, roi2_end = d_es - 30, d_es + 30
            
            # Branch Point: -15 to -45 upstream of 3'SS
            # [d_es - 45, d_es - 15]
            roi3_start, roi3_end = d_es - 45, d_es - 15
            
        else:  # Minus strand
            # For minus strand, coordinates are reversed in terms of biological meaning
            # 5'SS (Donor) is at downstreamES (genomic position)
            roi1_start, roi1_end = d_es - 30, d_es + 30
            
            # 3'SS (Acceptor) is at upstreamEE (genomic position)
            roi2_start, roi2_end = u_ee - 30, u_ee + 30
            
            # Branch Point: -15 to -45 upstream of 3'SS (which is upstreamEE on minus)
            # [u_ee + 15, u_ee + 45]
            roi3_start, roi3_end = u_ee + 15, u_ee + 45
        
        # Ensure valid coordinates
        roi1_start = max(0, roi1_start)
        roi2_start = max(0, roi2_start)
        roi3_start = max(0, roi3_start)
        
        beds['roi1_5ss'].append((chrom, roi1_start, roi1_end, name, 0, strand))
        beds['roi2_3ss'].append((chrom, roi2_start, roi2_end, name, 0, strand))
        beds['roi3_bp'].append((chrom, roi3_start, roi3_end, name, 0, strand))
    
    # Write BED files
    for roi_name, entries in beds.items():
        bed_path = f"{out_prefix}.{roi_name}.bed"
        with open(bed_path, 'w') as f:
            for entry in entries:
                f.write("\t".join(map(str, entry)) + "\n")
        print(f"[INFO] Wrote {len(entries)} entries to {bed_path}")
    
    return beds


def extract_fasta(bed_file, genome_fasta, output_fasta):
    """Extract sequences from genome using bedtools getfasta."""
    cmd = f"bedtools getfasta -s -fi {genome_fasta} -bed {bed_file} -fo {output_fasta} -name"
    run_cmd(cmd, f"Extracting sequences for {os.path.basename(bed_file)}")


def load_constitutive_introns(gtf_path, ri_chroms, n_sample=500):
    """
    Load constitutive introns from GTF as control group.
    Samples introns from protein-coding transcripts.
    """
    print(f"[INFO] Loading constitutive introns from GTF...")
    
    # Parse GTF for exons to derive introns
    exons = []
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'exon': continue
            
            chrom = parts[0].replace('chr', '')
            if chrom == 'M': chrom = 'MT'
            if chrom not in ri_chroms: continue  # Only consider chromosomes with RI events
            
            start, end = int(parts[3]) - 1, int(parts[4])  # 0-based
            strand = parts[6]
            
            # Extract transcript ID
            attrs = parts[8]
            tx_id = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    tx_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                    break
            
            if tx_id:
                exons.append({'chrom': chrom, 'start': start, 'end': end, 
                              'strand': strand, 'transcript_id': tx_id})
    
    if not exons:
        print("[WARN] No exons found in GTF.")
        return pd.DataFrame()
    
    exon_df = pd.DataFrame(exons)
    
    # Derive introns per transcript
    introns = []
    for tx_id, tx_exons in exon_df.groupby('transcript_id'):
        tx_exons = tx_exons.sort_values('start')
        for i in range(len(tx_exons) - 1):
            intron_start = tx_exons.iloc[i]['end']
            intron_end = tx_exons.iloc[i+1]['start']
            if intron_end > intron_start + 50:  # Minimum intron length
                introns.append({
                    'chr': tx_exons.iloc[i]['chrom'],
                    'upstreamEE': intron_start,
                    'downstreamES': intron_end,
                    'strand': tx_exons.iloc[i]['strand']
                })
    
    if not introns:
        print("[WARN] No introns derived from GTF.")
        return pd.DataFrame()
    
    intron_df = pd.DataFrame(introns)
    
    # Sample if too many
    if len(intron_df) > n_sample:
        intron_df = intron_df.sample(n=n_sample, random_state=42)
    
    print(f"[INFO] Loaded {len(intron_df)} constitutive introns as control.")
    return intron_df


def main():
    parser = argparse.ArgumentParser(description="ESE Motif Extraction (Step A)")
    parser.add_argument("--events_rds", required=True, help="Path to UFM1_events_rich.rds file")
    parser.add_argument("--gtf", required=True, help="Path to GTF annotation file")
    parser.add_argument("--genome_fasta", required=True, help="Path to genome FASTA file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--species", default="mouse", choices=["mouse", "human"])
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    # Load RDS using rpy2 or pre-converted TSV
    # For simplicity, convert RDS to TSV first via R
    tsv_path = os.path.join(args.outdir, "events_temp.tsv")
    r_cmd = f"""
    library(GenomicRanges)
    events <- readRDS("{args.events_rds}")
    df <- as.data.frame(events)
    write.table(df, "{tsv_path}", sep="\\t", quote=FALSE, row.names=FALSE)
    """
    
    print("[INFO] Converting RDS to TSV...")
    subprocess.run(["Rscript", "-e", r_cmd], check=True)
    
    # Load events
    events = pd.read_csv(tsv_path, sep="\t")
    
    # Filter for RI events only
    if "EventType" in events.columns:
        events = events[events["EventType"] == "RI"]
    
    print(f"[INFO] Loaded {len(events)} RI events.")
    
    # Split by group
    dependent = events[events["Group"] == "UFM1_dependent"]
    independent = events[events["Group"] == "UFM1_independent"]
    
    print(f"[INFO] UFM1-Dependent: {len(dependent)}, UFM1-Independent: {len(independent)}")
    
    # Get chromosomes from RI events for control filtering
    ri_chroms = set(events['chr'].str.replace('chr', '').replace('M', 'MT').unique()) if 'chr' in events.columns else set()
    
    # Load constitutive introns as control
    constitutive = load_constitutive_introns(args.gtf, ri_chroms)
    
    # Process each group
    groups = {
        "UFM1_dependent": dependent,
        "UFM1_independent": independent,
        "Constitutive": constitutive
    }
    
    for group_name, df in groups.items():
        if df.empty:
            print(f"[WARN] Skipping {group_name} (no events).")
            continue
        
        prefix = os.path.join(args.outdir, group_name)
        
        # Step 1: Extract ROI BED files
        extract_roi_beds(df, prefix, group_name)
        
        # Step 2: Extract FASTA sequences
        for roi in ['roi1_5ss', 'roi2_3ss', 'roi3_bp']:
            bed_file = f"{prefix}.{roi}.bed"
            fa_file = f"{prefix}.{roi}.fa"
            if os.path.exists(bed_file):
                extract_fasta(bed_file, args.genome_fasta, fa_file)
    
    print("[INFO] Step A complete. FASTA files generated.")
    print(f"[INFO] Output directory: {args.outdir}")


if __name__ == "__main__":
    main()
