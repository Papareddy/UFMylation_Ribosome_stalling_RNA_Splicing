#!/usr/bin/env python3
"""generate_genomic_SE_background.py
Scans a GTF, identifies all Exons, and classifies them using the same logic as SE analysis.
This provides the 'Genome Background' for Skipped Exon position distribution.
"""

import sys
import argparse
import pandas as pd
from collections import defaultdict
import os
import gzip
import re
from dataclasses import dataclass
from typing import List, Optional, Dict

@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int
    strand: str = "."
    def length(self) -> int: return max(0, self.end - self.start)
    def overlaps(self, other: "Interval") -> bool: return self.chrom == other.chrom and self.end > other.start and other.end > self.start

@dataclass
class TranscriptModel:
    gene_id_base: str; transcript_id: str; chrom: str; strand: str
    exons: List[Interval]; cds: List[Interval]; start_codons: List[Interval]; stop_codons: List[Interval]
    def sort_all(self) -> None:
        self.exons.sort(key=lambda x: (x.start, x.end)); self.cds.sort(key=lambda x: (x.start, x.end))
        self.start_codons.sort(key=lambda x: (x.start, x.end)); self.stop_codons.sort(key=lambda x: (x.start, x.end))
    def overlaps_any(self, iv: Interval, lst: List[Interval]) -> bool: return any(iv.overlaps(x) for x in lst)
    def cds_span(self) -> Optional[Interval]:
        if not self.cds: return None
        return Interval(self.chrom, min(c.start for c in self.cds), max(c.end for c in self.cds), self.strand)

def classify_interval_relative_to_transcript(iv: Interval, tm: TranscriptModel) -> str:
    # 1. Start/Stop Codon Overlap (Highest Priority)
    if tm.overlaps_any(iv, tm.start_codons): return "Start_Codon"
    if tm.overlaps_any(iv, tm.stop_codons): return "Stop_Codon"
    
    # 2. CDS Overlap
    if tm.overlaps_any(iv, tm.cds): return "CDS"
    
    # 3. UTR Logic
    c_span = tm.cds_span()
    if not c_span:
        if tm.overlaps_any(iv, tm.exons): return "NonCoding"
        return "Intronic" 
        
    # Check UTRs based on CDS span
    if tm.strand == '+':
        if iv.end <= c_span.start: return "5UTR"
        if iv.start >= c_span.end: return "3UTR"
    else: # '-'
        if iv.start >= c_span.end: return "5UTR"
        if iv.end <= c_span.start: return "3UTR"
        
    # If strictly inside CDS span (between start and stop) but overlaps nothing?
    # For Exons, they usually overlap 'exon' (themselves) or 'cds'.
    # If it's an exon but not CDS overlap, it's UTR.
    # The UTR logic above covers it.
    return "CDS"

_gtf_attr_re = re.compile(r'(\S+) \"([^\"]*)\";')
def parse_gtf_attributes(attr: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for m in _gtf_attr_re.finditer(attr): d[m.group(1)] = m.group(2)
    return d

def open_text_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def load_transcript_models(gtf_path: str) -> Dict[str, TranscriptModel]:
    by_tx: Dict[str, TranscriptModel] = {}
    print(f"[INFO] Loading GTF: {gtf_path}")
    with open_text_maybe_gzip(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith('#'): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9: continue
            chrom, _, feature, start_1, end_1, _, strand, _, attrs = parts
            chrom = chrom.replace("chr", "")
            if feature not in {'exon', 'CDS', 'start_codon', 'stop_codon'}: continue
            a = parse_gtf_attributes(attrs)
            gene_id, transcript_id = a.get('gene_id'), a.get('transcript_id')
            if not gene_id or not transcript_id: continue
            gene_id_base = gene_id.split('.')[0]
            start0, end0 = int(start_1) - 1, int(end_1)
            iv = Interval(chrom=chrom, start=start0, end=end0, strand=strand)
            
            if transcript_id not in by_tx:
                by_tx[transcript_id] = TranscriptModel(gene_id_base, transcript_id, chrom, strand, [], [], [], [])
            tm = by_tx[transcript_id]
            if feature == 'exon': tm.exons.append(iv)
            elif feature == 'CDS': tm.cds.append(iv)
            elif feature == 'start_codon': tm.start_codons.append(iv)
            elif feature == 'stop_codon': tm.stop_codons.append(iv)
    
    for tm in by_tx.values(): tm.sort_all()
    print(f"[INFO] Loaded {len(by_tx)} transcripts.")
    return by_tx

def main():
    parser = argparse.ArgumentParser(description="Generate Genome-wide Exon Distribution stats (for SE background).")
    parser.add_argument("--gtf", required=True, help="Path to GTF file")
    parser.add_argument("--out_tsv", required=True, help="Output TSV for background counts")
    args = parser.parse_args()

    tm_dict = load_transcript_models(args.gtf)
    
    # Logic:
    # 1. Collect all unique exons in the genome.
    # 2. For each unique exon, classify it against ALL transcripts that contain it or overlap it?
    #    Actually, exons are defined BY the transcript.
    #    An exon in Isoform A is UTR, in Isoform B is CDS.
    #    Wait, rMATS SE events are usually "Cassette Exons".
    #    The "Background" should be "All Exons that COULD be skipped".
    #    But "All Exons" is a good enough proxy for "Genomic Feature Distribution".
    #    So we iterate all unique (chrom, start, end, strand) intervals that are exons.
    
    exon_to_tms = defaultdict(list)
    
    print("[INFO] Extracting Exons...")
    count_exons = 0
    for tid, tm in tm_dict.items():
        for ex in tm.exons:
            # Key by interval
            exon_to_tms[ex].append(tm)
            count_exons += 1
            
    print(f"[INFO] Found {len(exon_to_tms)} unique exons (from {count_exons} total transcript-exons).")
    
    counts = defaultdict(int)
    
    print("[INFO] Classifying Exons...")
    for i, (iv, tms) in enumerate(exon_to_tms.items()):
        if i % 10000 == 0: print(f"... processed {i} unique exons")
        
        hits = set()
        for tm in tms:
            # Ensure the transcript actually contains this exon exactly?
            # Or just classification logic?
            # 'exon_to_tms' ensures tm has this exon.
            cls = classify_interval_relative_to_transcript(iv, tm)
            hits.add(cls)
        
        # Apply Priority (Same as splicing_functional_impat.py)
        final_cls = "Unknown"
        if "Start_Codon" in hits: final_cls = "Start_Codon"
        elif "Stop_Codon" in hits: final_cls = "Stop_Codon"
        elif "CDS" in hits: final_cls = "CDS"
        elif "3UTR" in hits: final_cls = "3UTR"
        elif "5UTR" in hits: final_cls = "5UTR"
        elif "NonCoding" in hits: final_cls = "NonCoding"
        elif "Intronic" in hits: final_cls = "Intronic" # Unlikely for an exon unless mapping error
        
        counts[final_cls] += 1
        
    print("[INFO] Writing counts...")
    with open(args.out_tsv, "w") as fh:
        fh.write("RI_Position\tcount\n")  # Header matches RI_Position for plot compatibility
        for cls, n in counts.items():
            fh.write(f"{cls}\t{n}\n")
    
    print(f"[DONE] Saved to {args.out_tsv}")

if __name__ == "__main__":
    main()
