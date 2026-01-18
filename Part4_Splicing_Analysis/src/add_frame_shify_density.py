#!/usr/bin/env python3
"""add_frame_shift_density.py

Compute relative CDS positions and plot density curves (matplotlib) along CDS (TSS→TES).

Key features
- Split and grouping options:
  - --event_types SE A3SS A5SS MXE RI
  - --pool (pool the listed event types)
  - --split_by_dataset (lost vs preserved)
- --overlay_datasets (overlay lost vs preserved within the same plot)
- --one_plot_per_class (one plot per protein_impact_class)
- Plot multiple protein_impact_class curves (or single via --plot_classes)
- Normalization modes for comparing lost vs preserved:
  - --normalize within  : each curve integrates to 1 (shape-only; default)
  - --normalize pooled  : scales each curve by group fraction (area reflects abundance)
  - --normalize ratio   : plots group density divided by pooled density (enrichment vs baseline)
- Bin-based statistics (Strategy 3):
  - --bin_stats writes per-subset tables comparing lost vs preserved per class using binned distributions

- CDS landmarks:
  - Start codon (CDS start) and stop codon (CDS end) are drawn as vertical dotted lines at x=0 and x=1.
  - Quartile landmarks are drawn at x=0.25, 0.5, 0.75.

Inputs
- per_event_compact_for_plotting.tsv (or '-' for stdin)
- GENCODE v45 GTF (hg38), .gtf or .gtf.gz

Outputs
- <out_prefix>*.png and <out_prefix>*.pdf
- Optional: <out_prefix>*.density.tsv and <out_prefix>*.per_event_with_relpos.tsv
- Optional bin stats TSVs: <out_prefix>*.binstats.tsv
"""

import argparse
import gzip
import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


@dataclass
class CDSModel:
    chrom: str
    strand: str
    intervals: List[Tuple[int, int]]  # 1-based inclusive
    cds_len: int


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_gtf_attributes(attr: str) -> Dict[str, str]:
    out = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr):
        out[m.group(1)] = m.group(2)
    return out


def _norm_chrom(chrom: str) -> str:
    if chrom is None:
        return ""
    c = str(chrom)
    if c.startswith("chr"):
        c = c[3:]
    return c


def load_cds_from_gtf(gtf_path: str) -> Tuple[Dict[str, CDSModel], Dict[str, CDSModel]]:
    """Return transcript and gene-level CDS models."""
    tx_to: Dict[str, CDSModel] = {}
    gene_to_intervals: Dict[str, Dict[str, List[Tuple[int, int]]]] = {}

    with open_text(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _frame, attr = parts
            if feature != "CDS":
                continue
            a = parse_gtf_attributes(attr)
            tx = a.get("transcript_id")
            gene = a.get("gene_id")
            if not tx and not gene:
                continue

            chrom_n = _norm_chrom(chrom)
            s = int(start)
            e = int(end)

            if tx:
                tx0 = tx.split(".")[0]
                if tx0 not in tx_to:
                    tx_to[tx0] = CDSModel(chrom=chrom_n, strand=strand, intervals=[], cds_len=0)
                tx_to[tx0].intervals.append((s, e))

            if gene:
                gene0 = gene.split(".")[0]
                if gene0 not in gene_to_intervals:
                    gene_to_intervals[gene0] = {"chrom": chrom_n, "strand": strand, "intervals": []}
                gene_to_intervals[gene0]["intervals"].append((s, e))

    for _tx, m in tx_to.items():
        m.intervals.sort(key=lambda x: x[0])
        m.cds_len = sum((e - s + 1) for s, e in m.intervals)

    gene_to: Dict[str, CDSModel] = {}
    for gene, rec in gene_to_intervals.items():
        ivs = sorted(rec["intervals"], key=lambda x: x[0])
        merged: List[Tuple[int, int]] = []
        for s, e in ivs:
            if not merged or s > merged[-1][1] + 1:
                merged.append((s, e))
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        cds_len = sum((e - s + 1) for s, e in merged)
        gene_to[gene] = CDSModel(chrom=rec["chrom"], strand=rec["strand"], intervals=merged, cds_len=cds_len)

    return tx_to, gene_to


def first_enst_from_list(s: str) -> Optional[str]:
    if s is None:
        return None
    s = str(s).strip()
    if not s:
        return None
    tx = s.split(";")[0].strip()
    if not tx:
        return None
    return tx.split(".")[0]


def parse_event_region(event_id: str) -> Optional[Tuple[str, int, int]]:
    """Parse event_id.

    Supports your pipe-delimited format:
      chr10|+|133168210|133168312|...
    Strategy: use [min(coords), max(coords)] as a coarse event region.
    """
    if event_id is None:
        return None
    s = str(event_id).strip()
    if not s:
        return None

    if "|" in s:
        parts = s.split("|")
        if len(parts) < 4:
            return None
        chrom = _norm_chrom(parts[0])
        nums = []
        for x in parts[2:]:
            try:
                nums.append(int(x))
            except Exception:
                pass
        if len(nums) < 2:
            return None
        return chrom, min(nums), max(nums)

    chrom_m = re.match(r"^(chr[^:]+|[^:]+):", s)
    if not chrom_m:
        return None
    chrom = _norm_chrom(chrom_m.group(1))
    blocks = re.findall(r"(\d+)-(\d+)", s)
    if not blocks:
        return None
    coords = []
    for a, b in blocks:
        coords.append(int(a))
        coords.append(int(b))
    return chrom, min(coords), max(coords)


def cds_position_for_genomic_point(cds: CDSModel, gpos: int) -> Optional[int]:
    if cds.cds_len <= 0:
        return None

    ivs = cds.intervals
    if cds.strand == "-":
        ivs = sorted(ivs, key=lambda x: x[0], reverse=True)

    cum = 0
    for s, e in ivs:
        if cds.strand == "+":
            if s <= gpos <= e:
                return cum + (gpos - s)
            cum += (e - s + 1)
        else:
            if s <= gpos <= e:
                return cum + (e - gpos)
            cum += (e - s + 1)
    return None


def midpoint_overlap_with_cds(chrom: str, start: int, end: int, cds: CDSModel) -> Optional[int]:
    if cds.chrom != _norm_chrom(chrom):
        return None
    overlaps = []
    for s, e in cds.intervals:
        os = max(start, s)
        oe = min(end, e)
        if os <= oe:
            overlaps.append((os, oe))
    if not overlaps:
        return None
    overlaps.sort(key=lambda x: (x[1] - x[0]), reverse=True)
    os, oe = overlaps[0]
    return (os + oe) // 2


def kde_gaussian(x: np.ndarray, grid: np.ndarray, bw: float) -> np.ndarray:
    if x.size == 0:
        return np.zeros_like(grid)
    x = x.reshape(-1, 1)
    g = grid.reshape(1, -1)
    z = (g - x) / bw
    dens = np.exp(-0.5 * z * z).sum(axis=0) / (x.shape[0] * bw * np.sqrt(2 * np.pi))
    return dens


def binned_fraction(x: np.ndarray, nbins: int) -> np.ndarray:
    """Return per-bin fractions over [0,1]."""
    if x.size == 0:
        return np.zeros(nbins, dtype=float)
    edges = np.linspace(0.0, 1.0, nbins + 1)
    # include 1.0 in last bin
    x2 = np.clip(x, 0.0, 1.0 - 1e-12)
    counts, _ = np.histogram(x2, bins=edges)
    tot = counts.sum()
    return counts / tot if tot > 0 else np.zeros(nbins, dtype=float)


def _logcomb(n: int, k: int) -> float:
    """log( n choose k ) using lgamma for stability."""
    from math import lgamma
    if k < 0 or k > n:
        return float("-inf")
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def fisher_exact_twosided(a: int, b: int, c: int, d: int):
    """Two-sided Fisher exact test for a 2x2 table [[a,b],[c,d]].
    Returns (odds_ratio, p_value).
    """
    from math import exp

    denom = b * c
    if denom == 0:
        odds = float("inf") if (a * d) > 0 else float("nan")
    else:
        odds = (a * d) / denom

    r1 = a + b
    r2 = c + d
    c1 = a + c
    n = r1 + r2

    lo = max(0, c1 - r2)
    hi = min(c1, r1)

    def log_p(x: int) -> float:
        # Hypergeometric pmf in log-space:
        return _logcomb(c1, x) + _logcomb(n - c1, r1 - x) - _logcomb(n, r1)

    logp_obs = log_p(a)
    p_obs = exp(logp_obs)

    p_two = 0.0
    for x in range(lo, hi + 1):
        px = exp(log_p(x))
        if px <= p_obs + 1e-15:
            p_two += px

    p_two = min(1.0, max(0.0, p_two))
    return odds, p_two


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR adjustment. Returns q-values in original order."""
    m = len(pvals)
    if m == 0:
        return []
    idx = np.argsort(pvals)
    q = np.empty(m, dtype=float)
    prev = 1.0
    for rank, i in enumerate(idx[::-1], start=1):
        k = m - rank + 1
        val = (m / k) * pvals[i]
        prev = min(prev, val)
        q[i] = prev
    return np.clip(q, 0.0, 1.0).tolist()


def star_from_q(q: float) -> str:
    if not np.isfinite(q):
        return "ns"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def binned_counts(x: np.ndarray, nbins: int):
    edges = np.linspace(0.0, 1.0, nbins + 1)
    if x.size == 0:
        return np.zeros(nbins, dtype=int), edges
    x2 = np.clip(x, 0.0, 1.0 - 1e-12)
    counts, _ = np.histogram(x2, bins=edges)
    return counts.astype(int), edges


def guess_contrast_pairs(present_datasets):
    pairs = []
    # Standard
    if "UFM1_dependent" in present_datasets and "UFM1_independent" in present_datasets:
        pairs.append(("UFM1_dependent", "UFM1_independent"))
    # Directional
    if "UFM1_dependent_dPSI_positive" in present_datasets and "UFM1_independent_dPSI_positive" in present_datasets:
        pairs.append(("UFM1_dependent_dPSI_positive", "UFM1_independent_dPSI_positive"))
    if "UFM1_dependent_dPSI_negative" in present_datasets and "UFM1_independent_dPSI_negative" in present_datasets:
        pairs.append(("UFM1_dependent_dPSI_negative", "UFM1_independent_dPSI_negative"))
    # Legacy (just in case)
    if "lost" in present_datasets and "preserved" in present_datasets:
        pairs.append(("lost", "preserved"))
    return pairs

def compute_bin_fisher(sub_df: pd.DataFrame, cls: str, nbins: int, ds1: str, ds2: str) -> Optional[pd.DataFrame]:
    """Per-bin Fisher (ds1 vs ds2) for one protein_impact_class."""
    unique_ds = sub_df["dataset"].unique().tolist()
    if (ds1 not in unique_ds) or (ds2 not in unique_ds):
        return None

    x_l = sub_df[(sub_df["dataset"] == ds1) & (sub_df["protein_impact_class"] == cls)]["relative_cds_position"].to_numpy(dtype=float)
    x_p = sub_df[(sub_df["dataset"] == ds2) & (sub_df["protein_impact_class"] == cls)]["relative_cds_position"].to_numpy(dtype=float)

    c_l, edges = binned_counts(x_l, nbins)
    c_p, _ = binned_counts(x_p, nbins)

    NL = int(c_l.sum())
    NP = int(c_p.sum())
    if NL == 0 or NP == 0:
        return None

    rows = []
    pvals = []
    odds_list = []
    centers = (edges[:-1] + edges[1:]) / 2.0

    for i in range(nbins):
        a = int(c_l[i])          # ds1 in bin
        b = NL - a               # ds1 out of bin
        c = int(c_p[i])          # ds2 in bin
        d = NP - c               # ds2 out of bin

        odds, p = fisher_exact_twosided(a, b, c, d)
        pvals.append(float(p))
        odds_list.append(float(odds))

        rows.append({
            "bin_index": i + 1,
            "bin_start": float(edges[i]),
            "bin_end": float(edges[i + 1]),
            "bin_center": float(centers[i]),
            f"n_{ds1}_in_bin": a,
            f"n_{ds1}_out_bin": b,
            f"n_{ds2}_in_bin": c,
            f"n_{ds2}_out_bin": d,
            f"odds_ratio_{ds1}_vs_{ds2}": float(odds),
            "p_value": float(p),
        })

    qvals = bh_fdr(pvals)

    for i in range(nbins):
        q = float(qvals[i])
        rows[i]["q_value"] = q
        rows[i]["star"] = star_from_q(q)
        odds = odds_list[i]
        if np.isfinite(odds):
            rows[i]["enriched_dataset"] = ds1 if odds > 1.0 else ds2
        else:
            rows[i]["enriched_dataset"] = "NA"

    return pd.DataFrame(rows)


def perm_pvalue_stat(
    x: np.ndarray,
    y: np.ndarray,
    nperm: int,
    stat_fn,
    rng: np.random.Generator,
    verbose: bool = False,
    label: str = "",
) -> float:
    """Two-sided permutation p-value for a scalar stat_fn(x,y)."""
    if x.size == 0 or y.size == 0:
        return float("nan")
    obs = float(stat_fn(x, y))
    if verbose:
        tag = f" [{label}]" if label else ""
        print(f"[INFO] Permutation test{tag}: running nperm={nperm}")
    z = np.concatenate([x, y])
    n_x = x.size
    more = 0
    for i in range(nperm):
        if verbose and (i > 0) and (i % 1000 == 0):
            tag = f" [{label}]" if label else ""
            print(f"[INFO] Permutation test{tag}: completed {i}/{nperm}")
        rng.shuffle(z)
        x_p = z[:n_x]
        y_p = z[n_x:]
        s = float(stat_fn(x_p, y_p))
        if abs(s) >= abs(obs):
            more += 1
    if verbose:
        tag = f" [{label}]" if label else ""
        print(f"[INFO] Permutation test{tag}: completed {nperm}/{nperm}")
    return (more + 1) / (nperm + 1)


def main():
    ap = argparse.ArgumentParser(description="Compute relative CDS positions and plot density curves (matplotlib).")
    ap.add_argument("--per_event", required=True, help="per_event_compact_for_plotting.tsv or '-' for stdin")
    ap.add_argument("--gtf", required=True, help="GENCODE v45 GTF (.gtf or .gtf.gz)")
    ap.add_argument("--out_prefix", required=True, help="Output prefix (dir/prefix)")
    ap.add_argument("--bw", type=float, default=0.08, help="KDE bandwidth on [0,1] scale")
    ap.add_argument("--grid_n", type=int, default=800, help="Number of x grid points")
    ap.add_argument(
        "--normalize",
        choices=["within", "pooled", "ratio", "log2ratio"],
        default="within",
        help=(
            "Normalization for comparing lost vs preserved: within|pooled|ratio|log2ratio. "
            "In ratio/log2ratio modes, f_pooled is computed from ALL events in the current subset "
            "(lost + preserved together)."
        ),
    )

    ap.add_argument("--write_per_event", action="store_true", help="Write per-event TSV with relpos")
    ap.add_argument("--write_density_tsv", action="store_true", help="Write x-grid + densities TSV")

    ap.add_argument("--split_by_dataset", action="store_true",
                    help="Separate plots/TSVs for dataset groups (lost, preserved), if present.")
    ap.add_argument("--dataset_order", default="UFM1_dependent,UFM1_independent",
                    help="Comma-separated dataset order (default: UFM1_dependent,UFM1_independent).")
    ap.add_argument(
        "--overlay_datasets",
        action="store_true",
        help="Overlay lost vs preserved curves within the same plot (requires 'dataset' column). This is the correct mode to compare positional differences between lost and preserved."
    )
    ap.add_argument(
        "--one_plot_per_class",
        action="store_true",
        help="If set, generate a separate plot for each protein_impact_class in --plot_classes (lost vs preserved overlaid if --overlay_datasets)."
    )

    ap.add_argument("--plot_classes",
                    default="CDS_neutral,Frameshift_NMD_likely,Frameshift_no_NMD,Inframe_CDS_change,Start_Stop_disruption,UTR_only",
                    help="Comma-separated protein_impact_class values to include as curves.")

    ap.add_argument("--event_types", nargs="*", default=None,
                    help="Space-separated EventType values to plot (e.g. SE A3SS A5SS MXE RI).")
    ap.add_argument("--pool", action="store_true",
                    help="If set together with --event_types, pool all listed event types into a single plot. If not set, produce one plot per listed event type.")

    ap.add_argument("--bin_stats", action="store_true",
                    help="If set, write bin-based statistics comparing lost vs preserved per class.")
    ap.add_argument("--nbins", type=int, default=10, help="Number of bins over [0,1] for bin-based stats (default: 10)")
    ap.add_argument("--nperm", type=int, default=5000, help="Permutation count for bin-based p-values (default: 5000)")
    ap.add_argument("--seed", type=int, default=1, help="RNG seed for permutation tests")
    ap.add_argument("--verbose", action="store_true", help="Verbose logging (e.g., bin_stats permutation progress).")
    ap.add_argument(
        "--annotate_binstats",
        action="store_true",
        help="If set with --bin_stats, add per-class bin-stat summary (Δmedian and permutation p-value) onto the plot."
    )
    ap.add_argument(
        "--bin_fisher",
        action="store_true",
        help=(
            "Bin-wise Fisher exact tests (lost vs preserved) over --nbins bins. "
            "Applies BH-FDR across bins and (for log2ratio + --overlay_datasets) annotates stars per bin."
        ),
    )
    ap.add_argument(
        "--write_bin_fisher_tsv",
        action="store_true",
        help="If set with --bin_fisher, write per-bin Fisher/FDR TSV (<out_prefix>.bin_fisher.tsv).",
    )

    args = ap.parse_args()

    if args.per_event == "-":
        df = pd.read_csv(sys.stdin, sep="\t", dtype=str)
    else:
        df = pd.read_csv(args.per_event, sep="\t", dtype=str)

    required = {"event_id", "GeneID", "protein_impact_class"}
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns in per_event file: {missing}")

    if "dataset" not in df.columns:
        raise SystemExit("Missing required column 'dataset' (expected values like lost/preserved).")

    if "EventType" not in df.columns:
        raise SystemExit("Missing required column 'EventType' (SE/RI/A3SS/A5SS/MXE).")

    tx_col = "advanced_transcript_ids" if "advanced_transcript_ids" in df.columns else None

    print(f"[INFO] Loading CDS from GTF: {args.gtf}")
    tx_cds_map, gene_cds_map = load_cds_from_gtf(args.gtf)
    print(f"[INFO] Loaded CDS models for {len(tx_cds_map):,} transcripts and {len(gene_cds_map):,} genes")

    relpos = []
    label_used = []

    for _, r in df.iterrows():
        ev = r.get("event_id", "")
        gene = str(r.get("GeneID", "")).split(".")[0]

        reg = parse_event_region(ev)
        if reg is None:
            relpos.append(np.nan)
            label_used.append("")
            continue

        chrom, st, en = reg

        cds = None
        label = ""
        if tx_col:
            tx = first_enst_from_list(r.get(tx_col, ""))
            if tx and tx in tx_cds_map:
                cds = tx_cds_map[tx]
                label = tx
        if cds is None and gene and gene in gene_cds_map:
            cds = gene_cds_map[gene]
            label = gene

        if cds is None:
            relpos.append(np.nan)
            label_used.append(label)
            continue

        mid = midpoint_overlap_with_cds(chrom, st, en, cds)
        if mid is None:
            relpos.append(np.nan)
            label_used.append(label)
            continue

        cpos = cds_position_for_genomic_point(cds, mid)
        if cpos is None or cds.cds_len <= 0:
            relpos.append(np.nan)
            label_used.append(label)
            continue

        rp = cpos / float(cds.cds_len)
        rp = min(1.0, max(0.0, rp))
        relpos.append(rp)
        label_used.append(label)

    df["relative_cds_position"] = pd.to_numeric(relpos, errors="coerce")
    df["relpos_model_used"] = label_used

    n_total = len(df)
    n_rel = int(np.isfinite(df["relative_cds_position"]).sum())
    print(f"[INFO] relative_cds_position computed for {n_rel}/{n_total} events")

    # Labels for legend/TSV columns
    label_map = {
        "Frameshift_NMD_likely": "NMD_likely",
        "Frameshift_no_NMD": "non_NMD",
        "CDS_neutral": "CDS_neutral",
        "Inframe_CDS_change": "Inframe_CDS_change",
        "Start_Stop_disruption": "Start_Stop_disruption",
        "UTR_only": "UTR_only",
    }

    classes_to_plot = [c.strip() for c in str(args.plot_classes).split(",") if c.strip()]
    classes_to_plot = [c for c in classes_to_plot if c in label_map]
    if not classes_to_plot:
        raise SystemExit("No valid classes to plot. Check --plot_classes.")

    dataset_order = [x.strip() for x in str(args.dataset_order).split(",") if x.strip()]

    base = df.copy()
    base = base[np.isfinite(base["relative_cds_position"])].copy()

    # user event types filter
    user_event_types = None
    if args.event_types is not None and len(args.event_types) > 0:
        user_event_types = [str(x).strip() for x in args.event_types if str(x).strip()]
        present = set(base["EventType"].dropna().unique().tolist())
        user_event_types = [et for et in user_event_types if et in present]
        if len(user_event_types) == 0:
            user_event_types = None

    # Decide event-group iteration: if --event_types given, either pool or split
    def _iter_event_groups():
        if user_event_types is None:
            yield ("ALL", None)  # no filtering
            return
        if args.pool:
            tag = "_".join(user_event_types)
            yield (f"pooled_{tag}", user_event_types)
        else:
            for et in user_event_types:
                yield (f"event_{et}", [et])

    # Decide dataset iteration
    def _iter_datasets():
        # If overlay is active, we must pass the full dataset to _run_one (ds=None)
        # so that it can plot both groups on the same axes.
        if args.overlay_datasets:
            yield None
            return
        if not args.split_by_dataset:
            yield None
            return
        present = set(base["dataset"].dropna().unique().tolist())
        ordered = [d for d in dataset_order if d in present]
        rest = sorted([d for d in present if d not in ordered])
        for d in ordered + rest:
            yield d

    rng = np.random.default_rng(args.seed)

    def _run_one(
        sub_df: pd.DataFrame,
        out_prefix: str,
        title: Optional[str] = None,
        classes_this_plot: Optional[List[str]] = None,
    ) -> None:
        sub_df = sub_df.copy()
        sub_df = sub_df[np.isfinite(sub_df["relative_cds_position"])].copy()
        classes_local = classes_to_plot if (classes_this_plot is None) else classes_this_plot
        sub_df = sub_df[sub_df["protein_impact_class"].isin(set(classes_local))].copy()
        if sub_df.empty:
            print(f"[WARN] No rows to plot for: {out_prefix}")
            return

        grid = np.linspace(0.0, 1.0, args.grid_n)

        # pooled baseline for ratio mode
        pooled_density = None
        if args.normalize in {"ratio", "log2ratio"}:
            # f_pooled is computed from ALL events in the current subset
            pooled_x = sub_df["relative_cds_position"].to_numpy(dtype=float)
            pooled_density = kde_gaussian(pooled_x, grid, bw=args.bw)
            pooled_density = np.maximum(pooled_density, 1e-12)

        # Compute per-(dataset,class) densities
        dens = {}  # (dataset, class) -> density
        counts = {}  # (dataset, class) -> n

        ds_levels = ["ALL"]
        present_ds = set(sub_df["dataset"].dropna().unique().tolist())
        
        if args.overlay_datasets:
            ds_levels = [d for d in dataset_order if d in present_ds]
            if not ds_levels:
                ds_levels = sorted(list(present_ds))
        elif args.split_by_dataset:
            ds_levels = [d for d in dataset_order if d in present_ds]
            if not ds_levels:
                ds_levels = sorted(list(present_ds))

        for ds in ds_levels:
            ds_df = sub_df if ds == "ALL" else sub_df[sub_df["dataset"] == ds]
            for cls in classes_local:
                x = ds_df.loc[ds_df["protein_impact_class"] == cls, "relative_cds_position"].to_numpy(dtype=float)
                counts[(ds, cls)] = int(x.size)
                dens[(ds, cls)] = kde_gaussian(x, grid, bw=args.bw)

        # normalization
        if args.normalize == "pooled":
            # scale by group fraction within sub_df (per ds)
            for ds in ds_levels:
                tot = sum(counts[(ds, c)] for c in classes_local)
                tot = max(tot, 1)
                for cls in classes_local:
                    frac = counts[(ds, cls)] / tot
                    dens[(ds, cls)] = dens[(ds, cls)] * frac

        if args.normalize == "ratio":
            for ds in ds_levels:
                for cls in classes_local:
                    dens[(ds, cls)] = dens[(ds, cls)] / pooled_density

        if args.normalize == "log2ratio":
            eps = 1e-6
            for ds in ds_levels:
                for cls in classes_local:
                    dens[(ds, cls)] = np.log2((dens[(ds, cls)] + eps) / (pooled_density + eps))

        # Plot
        plt.figure(figsize=(6.4, 3.9))
        ax = plt.gca()
        for cls in classes_local:
            if args.overlay_datasets or args.split_by_dataset:
                for ds in ds_levels:
                    lab = f"{label_map[cls]} ({ds})"
                    ax.plot(grid, dens[(ds, cls)], linewidth=2.0, label=lab)
            else:
                ax.plot(grid, dens[("ALL", cls)], linewidth=2.0, label=label_map[cls])

        ax.set_xlim(0, 1)
        ax.set_xlabel("")
        ylabel = "Density"
        if args.normalize == "pooled":
            ylabel = "Pooled probability density"
        elif args.normalize == "ratio":
            ylabel = "Relative density vs pooled"
        elif args.normalize == "log2ratio":
            ylabel = r"$\log_2\left(\frac{f_{dataset}(x)}{f_{pooled}(x)}\right)$"
        ax.set_ylabel(ylabel)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["TSS", "TES"])

        # CDS landmarks
        ax.axvline(0.0, color="black", linestyle=":", linewidth=1.2)
        ax.axvline(1.0, color="black", linestyle=":", linewidth=1.2)
        for q in (0.25, 0.5, 0.75):
            ax.axvline(q, color="gray", linestyle="--", linewidth=0.7)

        ax.text(0.0, 1.02, "Start (AUG)", transform=ax.get_xaxis_transform(), ha="left", va="bottom", fontsize=8)
        ax.text(1.0, 1.02, "Stop", transform=ax.get_xaxis_transform(), ha="right", va="bottom", fontsize=8)

        if title:
            ax.set_title(title)
        ax.legend(frameon=False, fontsize=8)

        if args.normalize == "log2ratio":
            ax.axhline(0.0, color="gray", linewidth=1.0, linestyle="--")

        # Identify pairs for stats
        possible_pairs = guess_contrast_pairs(present_ds)
        # Filter pairs: both members must be in ds_levels (which is derived from dataset_order)
        contrast_pairs = [p for p in possible_pairs if p[0] in ds_levels and p[1] in ds_levels]

        # Bin-wise Fisher
        bin_fisher_tables: List[pd.DataFrame] = []
        if args.bin_fisher:
            for cls in classes_local:
                for (ds1, ds2) in contrast_pairs:
                    t = compute_bin_fisher(sub_df, cls, args.nbins, ds1, ds2)
                    if t is None:
                        continue
                    t.insert(0, "comparison", f"{ds1}_vs_{ds2}")
                    t.insert(0, "label", label_map.get(cls, cls))
                    t.insert(0, "class", cls)
                    bin_fisher_tables.append(t)

            if args.write_bin_fisher_tsv and bin_fisher_tables:
                out_bf = f"{out_prefix}.bin_fisher.tsv"
                pd.concat(bin_fisher_tables, ignore_index=True).to_csv(out_bf, sep="\t", index=False)
                print(f"[DONE] Wrote: {out_bf}")

        # Plot stars
        if args.normalize == "log2ratio" and args.overlay_datasets and bin_fisher_tables:
            y0, y1 = ax.get_ylim()
            yr = y1 - y0
            # If multiple pairs, might need offsetting. For now, just plot first? 
            # Or simplified: top/bottom for first pair only to avoid clutter.
            # Assuming typically one pair is relevant per plot.
            
            # Group by comparison
            bf_df = pd.concat(bin_fisher_tables, ignore_index=True)
            comparisons = bf_df["comparison"].unique()
            
            # Heuristic offsets
            offsets = [(-0.04, 0.06), (-0.08, 0.10), (-0.12, 0.14)] # top, bottom relative factors
            
            for i, comp in enumerate(comparisons):
                if i >= len(offsets): break # limited space
                off_t, off_b = offsets[i]
                
                top_y = y1 + off_t * yr
                bot_y = y0 + off_b * yr
                
                sub_t = bf_df[bf_df["comparison"] == comp]
                
                for _, rr in sub_t.iterrows():
                    star = rr.get("star", "ns")
                    if star == "ns": continue
                    xc = float(rr["bin_center"])
                    enr = str(rr.get("enriched_dataset", "NA"))
                    
                    # Parse comp to know ds names
                    # comp format: "A_vs_B"
                    # enriched_dataset is one of them
                    parts = comp.split("_vs_")
                    ds1 = parts[0]
                    # if enr == ds1 -> top, else bottom
                    
                    if enr == ds1:
                        yy, va = top_y, "top"
                    else:
                        yy, va = bot_y, "bottom"
                    
                    ax.text(xc, yy, star, ha="center", va=va, fontsize=10)

        # Optional: annotate bin-stats on the plot
        if args.bin_stats and args.annotate_binstats:
            lines = ["Bin-stats (Δmedian, p):"]
            
            for (ds1, ds2) in contrast_pairs:
                for cls in classes_local:
                    x_1 = sub_df[(sub_df["dataset"] == ds1) & (sub_df["protein_impact_class"] == cls)]["relative_cds_position"].to_numpy(dtype=float)
                    x_2 = sub_df[(sub_df["dataset"] == ds2) & (sub_df["protein_impact_class"] == cls)]["relative_cds_position"].to_numpy(dtype=float)

                    if x_1.size == 0 or x_2.size == 0: continue

                    # medians
                    med_1 = float(np.median(x_1))
                    med_2 = float(np.median(x_2))
                    d_med = med_1 - med_2

                    def _stat_bins(a, b):
                        return float(np.sum(np.abs(binned_fraction(a, args.nbins) - binned_fraction(b, args.nbins))))

                    p_perm = perm_pvalue_stat(
                        x_1.copy(), x_2.copy(), args.nperm, _stat_bins, rng, verbose=args.verbose, label=f"{label_map[cls]} {ds1}v{ds2}"
                    )
                    
                    try:
                        p_txt = f"{p_perm:.3g}" if np.isfinite(p_perm) else "NA"
                    except: p_txt = "NA"
                    
                    try:
                        d_txt = f"{d_med:+.3f}" if np.isfinite(d_med) else "NA"
                    except: d_txt = "NA"
                    
                    star = star_from_q(p_perm)
                    s_str = f" {star}" if star != "ns" else ""
                    
                    # Tag with ds1 vs ds2 if multiple pairs
                    tag = f" ({ds1} vs {ds2})" if len(contrast_pairs) > 1 else ""
                    lines.append(f"{label_map[cls]}{tag}: {d_txt}, p={p_txt}{s_str} ({x_1.size},{x_2.size})")

            if len(lines) > 1:
                ax.text(
                    0.02, 0.98, "\n".join(lines),
                    transform=ax.transAxes, ha="left", va="top", fontsize=7.5, family="monospace"
                )

        plt.tight_layout()
        out_png = f"{out_prefix}.png"
        out_pdf = f"{out_prefix}.pdf"
        plt.savefig(out_png, dpi=200)
        plt.savefig(out_pdf)
        plt.close()
        print(f"[DONE] Wrote: {out_png}")
        print(f"[DONE] Wrote: {out_pdf}")

        if args.write_density_tsv:
            out_dt = f"{out_prefix}.density.tsv"
            out = pd.DataFrame({"x": grid})
            if args.overlay_datasets or args.split_by_dataset:
                for cls in classes_local:
                    for ds in ds_levels:
                        out[f"{label_map[cls]}__{ds}"] = dens[(ds, cls)]
            else:
                for cls in classes_local:
                    out[label_map[cls]] = dens[("ALL", cls)]
            out.to_csv(out_dt, sep="\t", index=False)
            print(f"[DONE] Wrote: {out_dt}")

        if args.write_per_event:
            out_pe = f"{out_prefix}.per_event_with_relpos.tsv"
            sub_df.to_csv(out_pe, sep="\t", index=False)
            print(f"[DONE] Wrote: {out_pe}")

        if args.bin_stats:
            # Write full stats table
            rows = []
            for (ds1, ds2) in contrast_pairs:
                for cls in classes_local:
                    x_1 = sub_df[(sub_df["dataset"] == ds1) & (sub_df["protein_impact_class"] == cls)]["relative_cds_position"].to_numpy(dtype=float)
                    x_2 = sub_df[(sub_df["dataset"] == ds2) & (sub_df["protein_impact_class"] == cls)]["relative_cds_position"].to_numpy(dtype=float)
                    
                    n_1 = int(x_1.size)
                    n_2 = int(x_2.size)
                    if n_1 == 0 or n_2 == 0: continue

                    f_1 = binned_fraction(x_1, args.nbins)
                    f_2 = binned_fraction(x_2, args.nbins)
                    l1 = float(np.sum(np.abs(f_1 - f_2)))

                    def _stat_bins(a, b):
                        return float(np.sum(np.abs(binned_fraction(a, args.nbins) - binned_fraction(b, args.nbins))))

                    p_perm = perm_pvalue_stat(
                         x_1.copy(), x_2.copy(), args.nperm, _stat_bins, rng, verbose=args.verbose, label=f"{label_map[cls]} {ds1}v{ds2}"
                    )
                    
                    med_1 = float(np.median(x_1))
                    med_2 = float(np.median(x_2))
                    d_med = med_1 - med_2
                    
                    rows.append({
                        "class": cls,
                        "label": label_map[cls],
                        "comparison": f"{ds1}_vs_{ds2}",
                        f"n_{ds1}": n_1,
                        f"n_{ds2}": n_2,
                        f"median_{ds1}": med_1,
                        f"median_{ds2}": med_2,
                        "delta_median": d_med,
                        "L1_bin_distance": l1,
                        "perm_pvalue_L1": p_perm,
                        "star": star_from_q(p_perm) # Approximate (p vs q? usually we use p here or do FDR across classes?)
                        # Keeping raw p-star for now as requested
                    })
                    
            if rows:
                out_bs = f"{out_prefix}.binstats.tsv"
                pd.DataFrame(rows).to_csv(out_bs, sep="\t", index=False)
                print(f"[DONE] Wrote: {out_bs}")

    # Run jobs
    n_jobs = 0
    for (etag, et_list) in _iter_event_groups():
        for ds in _iter_datasets():
            sub_df = base
            suffix_parts = []
            title_parts = ["CDS positional density"]

            if et_list is not None:
                sub_df = sub_df[sub_df["EventType"].isin(set(et_list))]
                suffix_parts.append(etag)
                title_parts.append(f"EventTypes={','.join(et_list)}" if args.pool else f"EventType={et_list[0]}")

            if ds is not None:
                sub_df = sub_df[sub_df["dataset"] == ds]
                suffix_parts.append(f"dataset_{ds}")
                title_parts.append(f"dataset={ds}")

            suffix = "__".join(suffix_parts)
            out_prefix = args.out_prefix if not suffix else f"{args.out_prefix}__{suffix}"
            title = " | ".join(title_parts) if suffix_parts else None

            # If --overlay_datasets is set, do NOT split jobs by dataset; comparison happens within the same plot.
            if args.overlay_datasets and ds is not None:
                continue

            # If we already split by dataset at job level, disable split_by_dataset inside plot
            # (so you get one curve set, not duplicated labels).
            old_split = args.split_by_dataset
            if ds is not None:
                args.split_by_dataset = False

            if args.one_plot_per_class:
                for cls in classes_to_plot:
                    cls_out = f"{out_prefix}__class_{label_map[cls]}"
                    cls_title = title
                    if cls_title:
                        cls_title = f"{cls_title} | class={label_map[cls]}"
                    else:
                        cls_title = f"class={label_map[cls]}"
                    _run_one(sub_df, cls_out, title=cls_title, classes_this_plot=[cls])
            else:
                _run_one(sub_df, out_prefix, title=title)

            args.split_by_dataset = old_split

            n_jobs += 1

    print(f"[INFO] Completed {n_jobs} plot job(s).")


if __name__ == "__main__":
    main()
