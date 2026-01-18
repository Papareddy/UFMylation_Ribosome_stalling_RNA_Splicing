#!/usr/bin/env python3
"""splicing_functional_impact.py
Robust functional annotation of pre-filtered rMATS events.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple


import pandas as pd
import numpy as np

import warnings

# Silence pandas FutureWarnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning, message=r".*Mean of empty slice.*",)


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int
    strand: str = "."
    def length(self) -> int: return max(0, self.end - self.start)
    def overlaps(self, other: "Interval") -> bool: return self.chrom == other.chrom and self.end > other.start and other.end > self.start
    def overlap_len(self, other: "Interval") -> int:
        if not self.overlaps(other): return 0
        return min(self.end, other.end) - max(self.start, other.start)

def safe_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")

def _split_event_id(event_id: str) -> List[str]:
    if event_id is None: return []
    return str(event_id).strip().split("|")

def row_from_event_id(event_type: str, event_id: str) -> Dict[str, object]:
    parts = _split_event_id(event_id)
    if len(parts) < 2: raise ValueError(f"Bad event_id (too few fields): {event_id}")
    chrom, strand = parts[0], parts[1]
    chrom = chrom.replace("chr", "")
    def as_int(x: str) -> int: return int(float(x))
    d: Dict[str, object] = {"chr": chrom, "strand": strand}
    et = str(event_type)
    if et == "SE":
        if len(parts) != 8: raise ValueError(f"SE event_id must have 8 fields, got {len(parts)}: {event_id}")
        d.update({"exonStart_0base": as_int(parts[2]), "exonEnd": as_int(parts[3]), "upstreamES": as_int(parts[4]), "upstreamEE": as_int(parts[5]), "downstreamES": as_int(parts[6]), "downstreamEE": as_int(parts[7])})
    elif et == "RI":
        if len(parts) != 8: raise ValueError(f"RI event_id must have 8 fields, got {len(parts)}: {event_id}")
        d.update({"riExonStart_0base": as_int(parts[2]), "riExonEnd": as_int(parts[3]), "upstreamES": as_int(parts[4]), "upstreamEE": as_int(parts[5]), "downstreamES": as_int(parts[6]), "downstreamEE": as_int(parts[7])})
    elif et in {"A3SS", "A5SS"}:
        if len(parts) != 8: raise ValueError(f"{et} event_id must have 8 fields, got {len(parts)}: {event_id}")
        d.update({"longExonStart_0base": as_int(parts[2]), "longExonEnd": as_int(parts[3]), "shortES": as_int(parts[4]), "shortEE": as_int(parts[5]), "flankingES": as_int(parts[6]), "flankingEE": as_int(parts[7])})
    elif et == "MXE":
        if len(parts) != 10: raise ValueError(f"MXE event_id must have 10 fields, got {len(parts)}: {event_id}")
        d.update({"1stExonStart_0base": as_int(parts[2]), "1stExonEnd": as_int(parts[3]), "2ndExonStart_0base": as_int(parts[4]), "2ndExonEnd": as_int(parts[5]), "upstreamES": as_int(parts[6]), "upstreamEE": as_int(parts[7]), "downstreamES": as_int(parts[8]), "downstreamEE": as_int(parts[9])})
    else: raise ValueError(f"Unsupported event_type for event_id parsing: {event_type}")
    return d

def load_merged_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if "EventType" in df.columns and "event_type" not in df.columns: df = df.rename(columns={"EventType": "event_type"})
    if "event_id" in df.columns and "ID" not in df.columns: df = df.rename(columns={"event_id": "ID"})
    if "GeneID" in df.columns:
        df["GeneID"] = df["GeneID"].astype(str).str.replace('"', "", regex=False).str.split(".", n=1).str[0]
    if "geneSymbol" in df.columns:
        df["geneSymbol"] = df["geneSymbol"].astype(str).str.replace('"', "", regex=False)
    for c in ("FDR_num.WT", "dPSI_num.WT", "FDR_num.UFM", "dPSI_num.UFM"):
        if c in df.columns: df[c] = safe_numeric(df[c])
    return df

_gtf_attr_re = re.compile(r'(\S+) \"([^\"]*)\";')
def parse_gtf_attributes(attr: str) -> Dict[str, str]:
    d: Dict[str, str] = {}
    for m in _gtf_attr_re.finditer(attr): d[m.group(1)] = m.group(2)
    return d

def open_text_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

@dataclass
class TranscriptModel:
    gene_id_base: str; transcript_id: str; chrom: str; strand: str
    exons: List[Interval]; cds: List[Interval]; start_codons: List[Interval]; stop_codons: List[Interval]
    def sort_all(self) -> None:
        self.exons.sort(key=lambda x: (x.start, x.end)); self.cds.sort(key=lambda x: (x.start, x.end))
        self.start_codons.sort(key=lambda x: (x.start, x.end)); self.stop_codons.sort(key=lambda x: (x.start, x.end))
    def transcript_span(self) -> Optional[Interval]:
        if not self.exons: return None
        return Interval(self.chrom, min(e.start for e in self.exons), max(e.end for e in self.exons), self.strand)
    def has_exon_exact(self, exon_iv: Interval) -> bool: return any(e.chrom == exon_iv.chrom and e.start == exon_iv.start and e.end == exon_iv.end for e in self.exons)
    def overlaps_any(self, iv: Interval, lst: List[Interval]) -> bool: return any(iv.overlaps(x) for x in lst)
    def cds_overlap_bp(self, iv: Interval) -> int: return int(sum(iv.overlap_len(c) for c in self.cds if c.chrom == iv.chrom))
    def exon_junctions(self) -> List[int]:
        if len(self.exons) < 2: return []
        ex = sorted(self.exons, key=lambda x: (x.start, x.end))
        return [ex[i].end for i in range(len(ex) - 1)]
    def final_junction_pos(self) -> Optional[int]:
        j = self.exon_junctions()
        if not j: return None
        return j[-1] if self.strand != '-' else j[0]
    
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
        # Non-coding transcript? Treat as "NonCoding" or check exons
        # If it overlaps exons but no CDS -> "NonCodingExon"
        if tm.overlaps_any(iv, tm.exons): return "NonCoding"
        return "Intronic" # or Intergenic
        
    # Check UTRs based on CDS span
    if tm.strand == '+':
        if iv.end <= c_span.start: return "5UTR"
        if iv.start >= c_span.end: return "3UTR"
    else: # '-'
        if iv.start >= c_span.end: return "5UTR"
        if iv.end <= c_span.start: return "3UTR"
        
    # If we are here, it didn't strictly fall outside CDS but didn't overlap CDS?
    # This implies it's in the INTROMS of the CDS region (but 'overlaps_any' checked CDS intervals).
    # IF the interval is inside the genomic span of the CDS but not overlapping CDS exons:
    # It must be an intronic region within the CDS range.
    # But usually RI *is* the intron.
    # If it overlaps CDS exons, we caught it.
    # If it is purely intronic relative to this transcript, return "Intronic"
    return "Intronic"

def load_transcript_models(gtf_path: str) -> Dict[str, List[TranscriptModel]]:
    by_gene_tx: Dict[str, Dict[str, TranscriptModel]] = defaultdict(dict)
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
            if transcript_id not in by_gene_tx[gene_id_base]:
                by_gene_tx[gene_id_base][transcript_id] = TranscriptModel(gene_id_base, transcript_id, chrom, strand, [], [], [], [])
            tm = by_gene_tx[gene_id_base][transcript_id]
            if feature == 'exon': tm.exons.append(iv)
            elif feature == 'CDS': tm.cds.append(iv)
            elif feature == 'start_codon': tm.start_codons.append(iv)
            elif feature == 'stop_codon': tm.stop_codons.append(iv)
    out: Dict[str, List[TranscriptModel]] = {}
    for g, txd in by_gene_tx.items():
        lst = list(txd.values())
        for tm in lst: tm.sort_all()
        out[g] = lst
    return out

def load_cds_intervals_by_gene(gtf_path: str) -> Dict[str, List[Interval]]:
    cds_by_gene: Dict[str, List[Interval]] = {}
    tx = load_transcript_models(gtf_path)
    for g, lst in tx.items():
        cds_all: List[Interval] = []
        for tm in lst: cds_all.extend(tm.cds)
        cds_by_gene[g] = sorted(cds_all, key=lambda x: (x.chrom, x.start, x.end))
    return cds_by_gene

def event_effective_exclusion_intervals(row: pd.Series, event_type: str) -> List[Interval]:
    _, exc = event_intervals(row, event_type)
    return [] if event_type in {'SE', 'RI'} else exc

def annotate_transcript_aware(df: pd.DataFrame, event_type: str, tx_by_gene: Dict[str, List[TranscriptModel]]) -> pd.DataFrame:
    out = df.copy()
    results = []
    for _, row in out.iterrows():
        gene_id = str(row.get('GeneID', '')).replace('"', '').split('.')[0]
        res = { 
            'advanced_n_transcripts_matched': 0, 
            'advanced_any_protein_pair': False, 
            'advanced_transcript_ids': "", 
            'inclusion_transcript_ids': "",
            'exclusion_transcript_ids': "",
            'advanced_overlaps_start_codon_any': False, 
            'advanced_overlaps_stop_codon_any': False, 
            'advanced_frameshift_any': False, 
            'advanced_frameshift_frac': 0.0, 
            'advanced_nmd_proxy_any': False,
            'RI_Position': 'NA'
        }
        try:
            inc_iv, _ = event_intervals(row, event_type)
            exc_eff = event_effective_exclusion_intervals(row, event_type)
        except Exception: 
            results.append(res); continue
        
        tx_list = tx_by_gene.get(gene_id, [])
        if not tx_list:
            results.append(res); continue
            
        matched_inc: List[TranscriptModel] = []
        matched_exc: List[TranscriptModel] = []
        
        for tm in tx_list:
            span = tm.transcript_span()
            # Basic overlap check with the event region
            if tm.chrom != inc_iv.chrom or tm.strand != inc_iv.strand or span is None or not span.overlaps(inc_iv): 
                continue
                
            # Inclusion Isoform Logic
            # SE/A3SS/A5SS/MXE: Must contain the inclusion exon(s) exactly
            # RI: "Inclusion" means retaining the intron (so riExon), "Exclusion" means splicing it out.
            #     RI logic in original code: if tm.has_exon_exact(inc_iv) -> "matched" (for whatever purpose).
            #     Wait, for RI, the 'inc_iv' is the RI exon (intron retention). 
            #     If a transcript has this exon, it IS the inclusion isoform.
            #     If a transcript has the upstream and downstream exons spliced together (skipping the RI part), it's exclusion.
            
            is_inclusion = False
            is_exclusion = False

            if event_type == 'SE':
                if tm.has_exon_exact(inc_iv):
                    is_inclusion = True
                # Exclusion: overlaps the genomic range but skips this specific exon.
                # Specifically for SE, check if it joins upstream and downstream directly?
                # Simpler proxy: overlaps the range but DOES NOT have the exon. 
                # Better: Check splice junctions. But let's stick to exon presence for now as proxy.
                elif not tm.has_exon_exact(inc_iv):
                     # Likely exclusion if it covers the region. We checked span.overlaps(inc_iv) above.
                     is_exclusion = True
            
            elif event_type == 'RI':
                 # inc_iv is the retained intron (as an exon).
                 if tm.has_exon_exact(inc_iv):
                     is_inclusion = True
                 else:
                     is_exclusion = True
            
            else:
                # A3SS, A5SS, MXE
                # Inclusion: has the "long" or "1st" exon
                if tm.has_exon_exact(inc_iv):
                    is_inclusion = True
                else:
                    is_exclusion = True

            if is_inclusion:
                matched_inc.append(tm)
            if is_exclusion:
                matched_exc.append(tm)

        # "Matched" usually meant inclusion-containing in old logic? 
        # Original: if (event_type == 'SE' and tm.has_exon_exact(inc_iv))... matched.append
        # So original 'matched' was strictly Inclusion isoforms.
        # We will keep 'matched' as inclusion for backward compatibility of 'advanced_transcript_ids' 
        # but populating new fields properly.
        
        # Original logic re-implementation for 'matched' list (Inclusion only)
        matched = matched_inc 

        res['advanced_n_transcripts_matched'] = len(matched)
        
        # Populate IDs
        def get_ids(tms):
            return ';'.join(sorted(set(str(t.transcript_id).split('.')[0] for t in tms if t.transcript_id)))
            
        if matched:
            res['advanced_transcript_ids'] = get_ids(matched)
            res['inclusion_transcript_ids'] = get_ids(matched)
            
            res['advanced_overlaps_start_codon_any'] = any(tm.overlaps_any(inc_iv, tm.start_codons) for tm in matched)
            res['advanced_overlaps_stop_codon_any'] = any(tm.overlaps_any(inc_iv, tm.stop_codons) for tm in matched)
            
            fs_flags, pair_flags, nmd_flags = [], [], []
            for tm in matched:
                inc_ov = tm.cds_overlap_bp(inc_iv)
                exc_ov = int(sum(tm.cds_overlap_bp(e) for e in exc_eff)) if exc_eff else 0
                d = inc_ov - exc_ov
                fs = (d != 0) and (abs(d) % 3 != 0)
                fs_flags.append(fs)
                
                # Pair logic: checking if this event makes sense in this transcript
                pair_flags.append(bool(inc_ov > 0 and exc_ov > 0) if event_type in {'A3SS', 'A5SS', 'MXE'} else bool(inc_ov > 0))
                
                fj = tm.final_junction_pos()
                if fj is not None:
                    upstream = (inc_iv.end + 55) < fj if tm.strand == '-' else (inc_iv.start + 55) < fj
                    nmd_flags.append(fs and upstream)
            
            res['advanced_any_protein_pair'] = any(pair_flags)
            res['advanced_frameshift_any'] = any(fs_flags)
            res['advanced_frameshift_frac'] = float(sum(1 for x in fs_flags if x) / max(1, len(fs_flags)))
            res['advanced_nmd_proxy_any'] = any(nmd_flags)
        
        if matched_exc:
            res['exclusion_transcript_ids'] = get_ids(matched_exc)
            
        # --- Positional Analysis (RI and SE) ---
        target_iv = None
        if event_type == 'RI':
            # Derive Intron Interval (upstreamEE to downstreamES)
            try:
                if all(k in row for k in ("upstreamEE", "downstreamES", "chr", "strand")):
                    uEE = int(row["upstreamEE"])
                    dES = int(row["downstreamES"])
                    if dES > uEE:
                        chrom = str(row["chr"]).replace("chr","")
                        strand = str(row["strand"])
                        target_iv = Interval(chrom, uEE, dES, strand)
            except Exception: pass
        elif event_type == 'SE':
            target_iv = inc_iv

        if target_iv:
            hits = set()
            for tm in tx_list:
                 cls = classify_interval_relative_to_transcript(target_iv, tm)
                 hits.add(cls)
            
            # Apply Priority
            pos_class = "Unknown"
            if "Start_Codon" in hits: pos_class = "Start_Codon"
            elif "Stop_Codon" in hits: pos_class = "Stop_Codon"
            elif "CDS" in hits: pos_class = "CDS"
            elif "3UTR" in hits: pos_class = "3UTR"
            elif "5UTR" in hits: pos_class = "5UTR"
            elif "NonCoding" in hits: pos_class = "NonCoding"
            elif "Intronic" in hits: pos_class = "Intronic"
            
            res['RI_Position'] = pos_class # Keep column name for compatibility
            res['Event_Position'] = pos_class

        results.append(res)
    return pd.concat([out.reset_index(drop=True), pd.DataFrame(results)], axis=1)

def iv(chrom: str, start0: int, end1: int, strand: str) -> Interval:
    return Interval(chrom=chrom, start=int(start0), end=int(end1), strand=strand)

def event_intervals(row: pd.Series, event_type: str) -> Tuple[Interval, List[Interval]]:
    chrom, strand = str(row["chr"]), str(row["strand"])
    if event_type == "SE":
        return iv(chrom, row["exonStart_0base"], row["exonEnd"], strand), [iv(chrom, row["upstreamES"], row["upstreamEE"], strand), iv(chrom, row["downstreamES"], row["downstreamEE"], strand)]
    if event_type == "RI":
        return iv(chrom, row["riExonStart_0base"], row["riExonEnd"], strand), [iv(chrom, row["upstreamES"], row["upstreamEE"], strand), iv(chrom, row["downstreamES"], row["downstreamEE"], strand)]
    if event_type in {"A3SS", "A5SS"}:
        return iv(chrom, row["longExonStart_0base"], row["longExonEnd"], strand), [iv(chrom, row["shortES"], row["shortEE"], strand)]
    if event_type == "MXE":
        c1s, c1e, c2s, c2e = ("1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd")
        return iv(chrom, row[c1s], row[c1e], strand), [iv(chrom, row[c2s], row[c2e], strand)]
    raise ValueError(f"Unsupported event_type: {event_type}")

def annotate_with_cds(df: pd.DataFrame, event_type: str, cds_by_gene: Dict[str, List[Interval]]) -> pd.DataFrame:
    out = df.copy()
    results = []
    for _, row in out.iterrows():
        res = {'inc_len_bp': 0, 'exc_len_bp': 0, 'inc_cds_overlap_bp': 0, 'exc_cds_overlap_bp': 0, 'inc_touches_cds': False, 'exc_touches_cds': False, 'inc_cds_overlap_mod3': None, 'exc_cds_overlap_mod3': None, 'delta_cds_overlap_bp': None, 'delta_cds_overlap_mod3': None, 'frameshift_candidate_inc': None, 'frameshift_candidate_delta': None}
        gene_id = str(row.get("GeneID", "")).replace('"', "").split(".")[0]
        try:
            inc_iv, exc_ivs = event_intervals(row, event_type)
            res['inc_len_bp'] = inc_iv.length()
            exc_ivs_eff = [] if event_type in {"SE", "RI"} else exc_ivs
            res['exc_len_bp'] = int(sum(iv_.length() for iv_ in exc_ivs_eff))
            cds_list = cds_by_gene.get(gene_id, [])
            if cds_list:
                inc_ov = int(sum(inc_iv.overlap_len(c) for c in cds_list if inc_iv.chrom == c.chrom))
                res['inc_cds_overlap_bp'] = inc_ov
                res['inc_touches_cds'] = inc_ov > 0
                if res['inc_touches_cds']:
                    m = int(inc_ov % 3)
                    res['inc_cds_overlap_mod3'] = m
                    res['frameshift_candidate_inc'] = bool(m != 0)
                
                exc_ov = 0
                if exc_ivs_eff:
                    exc_ov = int(sum(ex_iv.overlap_len(c) for ex_iv in exc_ivs_eff for c in cds_list if ex_iv.chrom == c.chrom))
                res['exc_cds_overlap_bp'] = exc_ov
                res['exc_touches_cds'] = exc_ov > 0
                if res['exc_touches_cds']:
                    res['exc_cds_overlap_mod3'] = int(exc_ov % 3)

                d = inc_ov - exc_ov
                res['delta_cds_overlap_bp'] = d
                if d != 0:
                    md = int(abs(d) % 3)
                    res['delta_cds_overlap_mod3'] = md
                    res['frameshift_candidate_delta'] = bool(md != 0)
                else:
                    res['delta_cds_overlap_mod3'] = 0
                    res['frameshift_candidate_delta'] = False
        except Exception: pass
        results.append(res)
    out = pd.concat([out.reset_index(drop=True), pd.DataFrame(results)], axis=1)
    out["cds_overlap_bp"] = out["inc_cds_overlap_bp"]
    out["touches_cds"] = out["inc_touches_cds"]
    out["cds_overlap_mod3"] = out["inc_cds_overlap_mod3"]
    out["frameshift_candidate"] = out["frameshift_candidate_inc"]
    return out

def summarize(df_all: pd.DataFrame, outdir: str):
    # This function can be simplified if not needed, or kept for summary stats
    pass

def main() -> None:
    p = argparse.ArgumentParser(description="Functional annotation of pre-filtered rMATS events.")
    p.add_argument("--lost_table", required=True, help="Merged TSV for 'lost' events.")
    p.add_argument("--preserved_table", required=True, help="Merged TSV for 'preserved' events.")
    p.add_argument("--gtf", required=True, help="GENCODE GTF (.gtf or .gtf.gz).")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--advanced", action="store_true", help="Add transcript-aware annotations.")
    p.add_argument("--event_types", nargs="+", default=["SE", "RI", "A3SS", "A5SS", "MXE"], help="Event types to process.")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if args.advanced:
        print(f"[INFO] Loading transcript models from GTF: {args.gtf}")
        tx_by_gene = load_transcript_models(args.gtf)
        cds_by_gene = {g: [c for tm in tms for c in tm.cds] for g, tms in tx_by_gene.items()}
        print(f"[INFO] Loaded CDS for {len(cds_by_gene):,} genes and transcripts for {len(tx_by_gene):,} genes")
    else:
        print(f"[INFO] Loading CDS from GTF: {args.gtf}")
        cds_by_gene = load_cds_intervals_by_gene(args.gtf)
        print(f"[INFO] Loaded CDS for {len(cds_by_gene):,} genes")

    dataset_tables = {"UFM1_dependent": args.lost_table, "UFM1_independent": args.preserved_table}
    all_datasets_rows: List[pd.DataFrame] = []

    for dataset_name, table_path in dataset_tables.items():
        print(f"[INFO] Loading merged table ({dataset_name}): {table_path}")
        df0 = load_merged_table(table_path)
        if "event_type" not in df0.columns: raise SystemExit("Merged table must contain EventType column.")
        
        df0 = df0[df0["event_type"].isin(args.event_types)].copy()
        if df0.empty:
            print(f"[WARN] No events of specified types found in {table_path}. Skipping.")
            continue
            
        if "ID" not in df0.columns: df0["ID"] = df0["event_id"]
        if "FDR" not in df0.columns: df0["FDR"] = df0.get("FDR_num.WT", pd.NA)
        if "IncLevelDifference" not in df0.columns: df0["IncLevelDifference"] = df0.get("dPSI_num.WT", pd.NA)

        coord_dicts = [row_from_event_id(r["event_type"], r["ID"]) for _, r in df0.iterrows()]
        df0 = pd.concat([df0.reset_index(drop=True), pd.DataFrame(coord_dicts)], axis=1)

        all_rows = []
        for ev in sorted(df0["event_type"].unique()):
            df_ev = df0[df0["event_type"] == ev].copy()
            print(f"[INFO] {dataset_name} / {ev}: {len(df_ev):,} events from merged_table")
            df_ann = annotate_with_cds(df_ev, ev, cds_by_gene)
            if args.advanced:
                df_ann = annotate_transcript_aware(df_ann, ev, tx_by_gene)
            df_ann["dataset"] = dataset_name
            all_rows.append(df_ann)
        
        if all_rows:
            df_dataset = pd.concat(all_rows, ignore_index=True)
            all_datasets_rows.append(df_dataset)
            ds_out = os.path.join(args.outdir, dataset_name)
            os.makedirs(ds_out, exist_ok=True)
            out_sig_annot = os.path.join(ds_out, "rmats_sig_annotated.tsv")
            df_dataset.to_csv(out_sig_annot, sep="\t", index=False)

    if all_datasets_rows:
        df_all = pd.concat(all_datasets_rows, ignore_index=True)
        # Final summary stats
        n = len(df_all)
        n_inc_cds = int(df_all["inc_touches_cds"].sum())
        n_fs_delta = int(df_all["frameshift_candidate_delta"].fillna(False).astype(bool).sum())
        print(f"[DONE] Wrote {n:,} total events")
        print(f"       inc touches CDS: {n_inc_cds:,} ({(n_inc_cds/n*100):.1f}%)")
        print(f"       frameshift_candidate_delta (heuristic): {n_fs_delta:,} ({(n_fs_delta/n*100):.1f}%)")
        if args.advanced:
            n_adv_fs = int(df_all.get('advanced_frameshift_any', pd.Series([False]*n)).fillna(False).astype(bool).sum())
            print(f"       advanced frameshift_any: {n_adv_fs:,} ({(n_adv_fs/n*100):.1f}%)")

if __name__ == "__main__":
    main()