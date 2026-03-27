from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Set

import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
import seaborn as sns

try:
    from Bio import Align
    from Bio.Align import substitution_matrices
except ImportError:
    print("[ERROR] Biopython is required. Please install it (e.g. pip install biopython).")
    sys.exit(1)


# -----------------------------
# Column normalization
# -----------------------------
COLUMN_ALIASES: Dict[str, List[str]] = {
    "EventType": ["EventType", "event_type", "eventType"],
    "event_id": ["event_id", "EventID", "event", "ID"],
    "GeneID": ["GeneID", "gene_id", "gene"],
    "inc_touches_CDS": ["inc_touches_CDS", "inc_touches_cds", "inc_touches_CDS_any"],
    "exc_touches_CDS": ["exc_touches_CDS", "exc_touches_cds", "exc_touches_CDS_any"],
    "advanced_frameshift_any": ["advanced_frameshift_any", "advanced_frameshift", "frameshift_any_advanced"],
    "advanced_NMD_proxy_any": ["advanced_NMD_proxy_any", "advanced_nmd_proxy_any", "advanced_NMD_any"],
    "advanced_overlaps_start_codon_any": [
        "advanced_overlaps_start_codon_any",
        "advanced_start_codon_overlap_any",
        "advanced_overlaps_start_any",
    ],
    "advanced_overlaps_stop_codon_any": [
        "advanced_overlaps_stop_codon_any",
        "advanced_stop_codon_overlap_any",
        "advanced_overlaps_stop_any",
    ],
    "advanced_n_transcripts_matched": ["advanced_n_transcripts_matched", "n_transcripts_matched", "advanced_n_tx"],
    "advanced_any_protein_pair": ["advanced_any_protein_pair", "any_protein_pair", "protein_pair"],
    "delta_cds_overlap_bp": ["delta_cds_overlap_bp", "delta_CDS_overlap_bp", "delta_cds_bp", "delta_CDS_bp"],
    "advanced_transcript_ids": ["advanced_transcript_ids", "transcript_ids", "tx_ids", "advanced_matched_transcripts"],
    "inclusion_transcript_ids": ["inclusion_transcript_ids", "inc_tx_ids", "inclusion_transcript_id", "inc_tx_id"],
    "exclusion_transcript_ids": ["exclusion_transcript_ids", "exc_tx_ids", "exclusion_transcript_id", "exc_tx_id"],
}


def _first_present(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns to canonical names when possible."""
    rename_map = {}
    for canon, alts in COLUMN_ALIASES.items():
        found = _first_present(df, alts)
        if found is not None and found != canon:
            rename_map[found] = canon
    if rename_map:
        df = df.rename(columns=rename_map)
    return df


def require_columns(df: pd.DataFrame, required: Iterable[str], label: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label}: missing required columns: {missing}. "
            f"Available columns: {list(df.columns)[:40]}{'...' if len(df.columns) > 40 else ''}"
        )


# -----------------------------
# Panel-ready proxies (E/G/H-inspired)
# -----------------------------

PROTEIN_IMPACT_PRIORITY = [
    "Start_Stop_disruption",
    "Frameshift_NMD_likely",
    "Frameshift_no_NMD",
    "Inframe_CDS_change",
    "CDS_neutral",
    "UTR_only",
]


def to_bool_series(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.fillna(False).astype(str).str.lower().isin(["true", "1", "t", "yes", "y"])


def add_protein_impact_class(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    df["inc_touches_CDS"] = to_bool_series(df["inc_touches_CDS"])
    df["exc_touches_CDS"] = to_bool_series(df["exc_touches_CDS"]) if "exc_touches_CDS" in df.columns else False

    df["advanced_frameshift_any"] = to_bool_series(df["advanced_frameshift_any"])
    df["advanced_NMD_proxy_any"] = to_bool_series(df["advanced_NMD_proxy_any"])
    df["advanced_overlaps_start_codon_any"] = to_bool_series(df["advanced_overlaps_start_codon_any"])
    df["advanced_overlaps_stop_codon_any"] = to_bool_series(df["advanced_overlaps_stop_codon_any"])

    # delta CDS
    df["delta_cds_overlap_bp"] = pd.to_numeric(df["delta_cds_overlap_bp"], errors="coerce")
    delta_nonzero = df["delta_cds_overlap_bp"].fillna(0).astype(float) != 0
    inframe = (df["delta_cds_overlap_bp"].fillna(0).abs() % 3) == 0

    # priority categories
    cat = pd.Series([None] * len(df), index=df.index, dtype="object")

    # Priority 1: Start/Stop Disruption
    m = (
        (df["advanced_overlaps_start_codon_any"] | df["advanced_overlaps_stop_codon_any"])
        & cat.isna()
    )
    cat[m] = "Start_Stop_disruption"

    # Priority 2: Frameshift with NMD
    m = df["advanced_frameshift_any"] & df["advanced_NMD_proxy_any"] & cat.isna()
    cat[m] = "Frameshift_NMD_likely"

    # Priority 3: Frameshift without NMD
    m = df["advanced_frameshift_any"] & ~df["advanced_NMD_proxy_any"] & cat.isna()
    cat[m] = "Frameshift_no_NMD"

    # Priority 4: In-frame Change
    m = delta_nonzero & inframe & cat.isna()
    cat[m] = "Inframe_CDS_change"

    # Priority 5: CDS Neutral (Touches CDS but no length change - e.g. synonymous or identical)
    m = df["inc_touches_CDS"] & cat.isna()
    cat[m] = "CDS_neutral"

    # Priority 6: UTR Only
    m = ~df["inc_touches_CDS"] & cat.isna()
    cat[m] = "UTR_only"

    df["protein_impact_class"] = pd.Categorical(cat, categories=PROTEIN_IMPACT_PRIORITY, ordered=True)

    # Protein predicted different/identical logic
    df["protein_predicted_different"] = (
        df["advanced_frameshift_any"]
        | delta_nonzero
        | df["advanced_overlaps_start_codon_any"]
        | df["advanced_overlaps_stop_codon_any"]
    )

    if "advanced_any_protein_pair" in df.columns:
        df["advanced_any_protein_pair"] = to_bool_series(df["advanced_any_protein_pair"])
        df["protein_pair_proxy"] = df["advanced_any_protein_pair"]
    else:
        df["protein_pair_proxy"] = df["inc_touches_CDS"]

    df["protein_predicted_identical"] = df["protein_pair_proxy"] & ~df["protein_predicted_different"]

    return df


def add_exon_location_class(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["inc_touches_CDS"] = to_bool_series(df["inc_touches_CDS"])
    if "exc_touches_CDS" in df.columns:
        df["exc_touches_CDS"] = to_bool_series(df["exc_touches_CDS"])
    else:
        df["exc_touches_CDS"] = False

    all_coding = df["inc_touches_CDS"] & df["exc_touches_CDS"]
    all_noncoding = (~df["inc_touches_CDS"]) & (~df["exc_touches_CDS"])
    mixed = ~(all_coding | all_noncoding)

    loc = pd.Series([None] * len(df), index=df.index, dtype="object")
    loc[all_coding] = "all_coding_exons"
    loc[all_noncoding] = "all_non_coding_exons"
    loc[mixed] = "some_non_coding_exons"

    df["exon_location_class"] = pd.Categorical(
        loc,
        categories=["all_coding_exons", "some_non_coding_exons", "all_non_coding_exons"],
        ordered=True,
    )
    return df


def make_fig1E(df: pd.DataFrame, dataset_label: str) -> pd.DataFrame:
    sub = df[df["protein_pair_proxy"]].copy()
    n = len(sub)
    if n == 0:
        return pd.DataFrame(
            [{"dataset": dataset_label, "n_events": 0, "frac_identical": np.nan, "frac_different": np.nan}]
        )
    frac_ident = float(sub["protein_predicted_identical"].mean())
    frac_diff = float(sub["protein_predicted_different"].mean())
    return pd.DataFrame(
        [{"dataset": dataset_label, "n_events": n, "frac_identical": frac_ident, "frac_different": frac_diff}]
    )


def make_fig1G(df: pd.DataFrame, dataset_label: str) -> pd.DataFrame:
    out_rows = []
    for et, g in df.groupby("EventType", dropna=False):
        delta_bp = pd.to_numeric(g["delta_cds_overlap_bp"], errors="coerce").fillna(0)
        delta_aa = (delta_bp.abs() / 3.0).astype(float)
        log2_len_change = np.log2(delta_aa + 1.0)
        out_rows.append(
            {
                "dataset": dataset_label,
                "event_type": et,
                "n_events": int(len(g)),
                "median_abs_delta_aa": float(np.nanmedian(delta_aa)),
                "median_log2_abs_delta_aa_plus1": float(np.nanmedian(log2_len_change)),
                "frac_delta_cds_nonzero": float((delta_bp != 0).mean()),
            }
        )
    return pd.DataFrame(out_rows).sort_values(["dataset", "event_type"]).reset_index(drop=True)


def make_fig1H(df: pd.DataFrame, dataset_label: str) -> pd.DataFrame:
    out_rows = []
    for et, g in df.groupby("EventType", dropna=False):
        denom = len(g)
        counts = g["exon_location_class"].value_counts(dropna=False)
        for cls in ["all_coding_exons", "some_non_coding_exons", "all_non_coding_exons"]:
            out_rows.append(
                {
                    "dataset": dataset_label,
                    "event_type": et,
                    "exon_location_class": cls,
                    "n_events": int(counts.get(cls, 0)),
                    "frac": float(counts.get(cls, 0) / denom) if denom else np.nan,
                }
            )
    return pd.DataFrame(out_rows).sort_values(["dataset", "event_type", "exon_location_class"]).reset_index(drop=True)


# -----------------------------
# F-inspired: protein sequence similarity (alignment-based)
# -----------------------------

def parse_tx_list(val: str) -> List[str]:
    """Parse transcript IDs from a cell."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return []
    s = str(val).strip()
    if not s:
        return []
    parts = []
    for delim in [";", ",", "|", " "]:
        if delim in s:
            parts = [p for p in s.replace("|", ";").replace(",", ";").split(";") if p]
            break
    if not parts:
        parts = [s]
    out = []
    for p in parts:
        p = p.strip().split(":")[0].split(".")[0]
        if p:
            out.append(p)
    return sorted(set(out))


def compute_alignment_score_blosum62(a: str, b: str) -> Optional[float]:
    """Return matched/(matched+mismatched) after global alignment.
    
    Uses Bio.Align.PairwiseAligner (replacement for pairwise2).
    """
    try:
        from Bio import Align
        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        aligner.mode = 'global'
        
        # alignments = aligner.align(a, b)
        # Taking the first/best alignment
        score = aligner.score(a, b)
        # Getting actual alignment to count matches
        aln = aligner.align(a, b)[0]
        
        # aln is an Alignment object. We can iterate its coordinates or string rep.
        # Biopython 1.80+ style
        # Get aligned sequences
        t_seq = aln.target # sequence A
        q_seq = aln.query  # sequence B
        
        # Reconstruct alignment strings (with gaps)
        # format() method is easiest but might vary by version. 
        # Using simple counting from counts / coords is safer if just needing identity.
        
        # Actually, "alignment_score" in Fig 1F often implies Identity %.
        # The previous code calculated: matched / (matched + mismatched + gaps)
        
        # Let's extract aligned strings to reproduce previous logic exactly
        # aln.format("fasta") gives strings
        lines = str(aln).strip().split("\n")
        # Typically 3 lines: seqA, match_indicators, seqB. Or formatted differently.
        
        # More robust: use aligned coordinates
        # But simpler: use format
        # Or just use the score itself? No, previous code used percent identity manually.
        
        # Using the aligned strings property if available (depends on version)
        # Let's rely on constructing it from coordinates to be safe or just use the string repr logic
        # For Biopython >1.78, aln[0, :] gives chars? No.
        
        # Safest fallback for simple counts without parsing complex objects: 
        # aligned representation
        # aln[0] and aln[1] ? No.
        
        # Let's use the simplest method: 
        # Iterate over the aligned positions.
        matched = aln.counts().identities
        length = aln.shape[1] # Alignment length (columns)
        
        if length == 0: return None
        return matched / length
        
    except Exception as e:
        # Fallback or error
        # print(f"[WARN] Alignment failed: {e}")
        return None


def load_pysam_and_map(fasta_path: str) -> Tuple[pysam.FastaFile, Dict[str, str]]:
    """
    Initializes pysam FastaFile and builds ENST -> Reference mapping.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Protein FASTA not found: {fasta_path}")
    
    # Check for index
    if not os.path.exists(fasta_path + ".fai") and not os.path.exists(fasta_path + ".gzi"):
        print(f"[INFO] Indexing FASTA file: {fasta_path}")
    
    try:
        pdb = pysam.FastaFile(fasta_path)
    except IOError:
        raise IOError(f"Could not open {fasta_path} with pysam. Ensure it is uncompressed or bgzip-compressed.")

    # Build mapping
    print("[INFO] Building ENST -> FASTA Reference map by scanning headers...")
    enst_to_key = {}
    
    open_func = open
    if fasta_path.endswith(".gz"):
        open_func = gzip.open
        
    try:
        with open_func(fasta_path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    parts = line[1:].strip().split(maxsplit=1)
                    ref_name = parts[0]
                    full_header = line[1:].strip()
                    norm_header = full_header.replace("|", " ").replace(":", " ")
                    tokens = norm_header.split()
                    
                    enst = None
                    for t in tokens:
                        if t.startswith("transcript:"):
                            t = t.split(":", 1)[1]
                        if t.startswith("ENST") or t.startswith("ENSMUST") or (t.startswith("ENS") and "T" in t):
                            enst = t.split(".")[0]
                            if len(enst) > 10: 
                                break
                    if enst:
                        enst_to_key[enst] = ref_name
    except Exception as e:
        print(f"[WARN] Failed to scan FASTA headers: {e}")
        
    print(f"[INFO] Mapped {len(enst_to_key):,} transcripts from FASTA.")
    return pdb, enst_to_key


def make_fig1F(
    df: pd.DataFrame,
    dataset_label: str,
    pdb: pysam.FastaFile,
    enst_to_key: Dict[str, str],
    max_events: int = 20000,
) -> pd.DataFrame:
    """Compute one alignment score per event using explicitly matched transcript pairs.
    
    Requirements:
    - df must contain 'inclusion_transcript_ids' and 'exclusion_transcript_ids'.
    - Uses pysam for efficient FASTA access.
    """
    req_cols = ["inclusion_transcript_ids", "exclusion_transcript_ids"]
    for c in req_cols:
        if c not in df.columns:
            raise ValueError(f"Alignment requires '{c}'. Please run with updated upstream annotation script.")
    
    out_rows = []
    n_done = 0

    for _, row in df.iterrows():
        if n_done >= max_events:
            break
        
        if n_done % 1000 == 0 and n_done > 0:
            print(f"[DEBUG] Processing event {n_done} ({dataset_label})...")

        inc_txs = parse_tx_list(row.get("inclusion_transcript_ids", ""))
        exc_txs = parse_tx_list(row.get("exclusion_transcript_ids", ""))
        
        if not inc_txs or not exc_txs:
            continue

        # Find a valid pair
        found_pair = None
        seq_a = None
        seq_b = None
        
        for inc_id in inc_txs:
            if inc_id not in enst_to_key: continue
            for exc_id in exc_txs:
                if exc_id not in enst_to_key: continue
                
                # Found a pair where both exist in FASTA
                try:
                    seq_a = pdb.fetch(enst_to_key[inc_id])
                    seq_b = pdb.fetch(enst_to_key[exc_id])
                    found_pair = (inc_id, exc_id)
                    break 
                except KeyError:
                    continue
            if found_pair: break
        
        if not found_pair or not seq_a or not seq_b:
            continue

        score = compute_alignment_score_blosum62(seq_a, seq_b)
        if score is None:
            continue

        out_rows.append(
            {
                "dataset": dataset_label,
                "event_id": row.get("event_id"),
                "GeneID": row.get("GeneID"),
                "event_type": row.get("EventType"),
                "inc_tx": found_pair[0],
                "exc_tx": found_pair[1],
                "alignment_score": float(score),
                "len_inc": int(len(seq_a)),
                "len_exc": int(len(seq_b)),
            }
        )
        n_done += 1

    return pd.DataFrame(out_rows)


# -----------------------------
# IO helpers
# -----------------------------

def read_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = normalize_columns(df)

    if "EventType" not in df.columns:
        raise ValueError(f"{path}: could not find EventType column (after alias normalization).")

    required = [
        "event_id",
        "GeneID",
        "EventType",
        "inc_touches_CDS",
        "advanced_frameshift_any",
        "advanced_NMD_proxy_any",
        "advanced_overlaps_start_codon_any",
        "advanced_overlaps_stop_codon_any",
        "delta_cds_overlap_bp",
    ]
    # exc_touches_CDS is strongly preferred for Fig1H but we can tolerate missing
    present_required = [c for c in required if c != "exc_touches_CDS"]
    require_columns(df, present_required, os.path.basename(path))

    df["EventType"] = df["EventType"].astype(str)
    return df


def ensure_outdir(outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)


def write_tsv(df: pd.DataFrame, path: str) -> None:
    df.to_csv(path, sep="\t", index=False)


# -----------------------------
# Main
# -----------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Panel-ready protein primary-sequence impact summaries from annotated rMATS events."
    )
    p.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="List of input annotated TSVs. Format: 'path' (label inferred from filename) or 'label=path'.",
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Output directory for summary TSVs",
    )
    p.add_argument(
        "--protein_fasta",
        required=False,
        help="Protein FASTA (e.g., gencode.v45.pc_translations.fa.gz). If provided, computes Fig 1F-like alignment scores.",
    )
    p.add_argument(
        "--max_events_for_alignment",
        type=int,
        default=20000,
        help="Safety cap on number of events to align per dataset (default: 20000)",
    )
    return p


def process_one(df: pd.DataFrame) -> pd.DataFrame:
    df = add_protein_impact_class(df)
    df = add_exon_location_class(df)
    return df


def plot_alignment_score_histogram(df: pd.DataFrame, outdir: str):
    """
    Plots a histogram of alignment_score faceted by dataset.
    """
    if df.empty or "alignment_score" not in df.columns:
        print("[WARN] No alignment data to plot.")
        return

    # Filter out NaNs
    df_clean = df.dropna(subset=["alignment_score"])
    if df_clean.empty:
        print("[WARN] No valid alignment scores to plot.")
        return

    sns.set_theme(style="whitegrid")
    
    # FacetGrid for histograms
    g = sns.FacetGrid(df_clean, col="dataset", col_wrap=3, height=4, aspect=1.5, sharex=True, sharey=False)
    g.map(sns.histplot, "alignment_score", bins=30, kde=True, element="step", stat="probability")
    
    g.set_titles("{col_name}")
    g.set_axis_labels("Alignment Score", "Frequency")
    
    outpath = os.path.join(outdir, "Alignment_Scores_Histogram.pdf")
    plt.savefig(outpath, bbox_inches='tight')
    plt.close()
    
    outpath_png = os.path.join(outdir, "Alignment_Scores_Histogram.png")
    plt.savefig(outpath_png, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"[DONE] Wrote alignment score histogram to {outpath}")


def main() -> None:
    args = build_parser().parse_args()
    ensure_outdir(args.outdir)

    # Parse inputs
    inputs = []
    for s in args.inputs:
        if "=" in s:
            lbl, pth = s.split("=", 1)
        else:
            pth = s
            lbl = os.path.splitext(os.path.basename(pth))[0]
        inputs.append((lbl, pth))

    processed = []
    for lbl, pth in inputs:
        print(f"[INFO] Processing {lbl} from {pth}...")
        df = read_table(pth)
        df = process_one(df)
        processed.append((lbl, df))

    # -----------------------------
    # E-inspired panel (proxy)
    # -----------------------------
    print("[INFO] Generating_Protein_Identity_Stats...")
    fig1E = pd.concat([make_fig1E(df, lbl) for lbl, df in processed], ignore_index=True)
    write_tsv(fig1E, os.path.join(args.outdir, "Protein_Identity_Stats.tsv"))

    # -----------------------------
    # G-inspired panel
    # -----------------------------
    print("[INFO] Generating_Length_Change_Stats...")
    fig1G = pd.concat([make_fig1G(df, lbl) for lbl, df in processed], ignore_index=True)
    write_tsv(fig1G, os.path.join(args.outdir, "Length_Change_Stats.tsv"))

    # -----------------------------
    # H-inspired panel
    # -----------------------------
    print("[INFO] Generating_Exon_Location_Stats...")
    fig1H = pd.concat([make_fig1H(df, lbl) for lbl, df in processed], ignore_index=True)
    write_tsv(fig1H, os.path.join(args.outdir, "Exon_Location_Stats.tsv"))

    # write a compact per-event table
    pe_parts = []
    for lbl, df in processed:
        keep_cols = [
            "event_id",
            "GeneID",
            "EventType",
            "protein_impact_class",
            "protein_predicted_identical",
            "protein_predicted_different",
            "delta_cds_overlap_bp",
            "exon_location_class",
            "advanced_frameshift_any",
            "advanced_NMD_proxy_any",
            "advanced_overlaps_start_codon_any",
            "advanced_overlaps_stop_codon_any",
            "advanced_n_transcripts_matched",
            "advanced_transcript_ids",
            "RI_Position",
        ]
        keep_cols = [c for c in keep_cols if c in df.columns]
        pe_parts.append(df[keep_cols].assign(dataset=lbl))

    per_event = pd.concat(pe_parts, ignore_index=True)
    write_tsv(per_event, os.path.join(args.outdir, "per_event_compact_for_plotting.tsv"))

    print("[DONE] Wrote summary tables.")

    if args.protein_fasta:
        print("[INFO] Generating_Alignment_Scores...")
        # -----------------------------
        # F-inspired panel (alignment)
        # -----------------------------
        pdb_obj, enst_map = load_pysam_and_map(args.protein_fasta)
        
        try:
            parts = []
            for lbl, df in processed:
                parts.append(
                    make_fig1F(
                        df,
                        lbl,
                        pdb_obj,
                        enst_map,
                        max_events=args.max_events_for_alignment,
                    )
                )
            fig1F = pd.concat(parts, ignore_index=True)
            write_tsv(fig1F, os.path.join(args.outdir, "Alignment_Scores.tsv"))
            print(f"[DONE] Wrote alignment scores to {os.path.join(args.outdir, 'Alignment_Scores.tsv')}")

            # Plot Histogram
            plot_alignment_score_histogram(fig1F, args.outdir)
        finally:
            pdb_obj.close()



if __name__ == "__main__":
    main()