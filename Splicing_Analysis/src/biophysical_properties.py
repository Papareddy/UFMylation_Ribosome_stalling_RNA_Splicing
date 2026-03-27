#!/usr/bin/env python3
import argparse, os, re
import pandas as pd

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

def clean_seq(s: str) -> str:
    s = s.upper().replace("*", "")
    s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s)
    return s

def net_charge_ph7(seq: str) -> float:
    pa = ProteinAnalysis(seq)
    if hasattr(pa, "charge_at_pH"):
        return float(pa.charge_at_pH(7.0))
    return float("nan")

def aliphatic_index(seq: str) -> float:
    pa = ProteinAnalysis(seq)
    comp = pa.get_amino_acids_percent()
    return 100.0 * (comp.get("A", 0) + 2.9 * comp.get("V", 0) + 3.9 * (comp.get("I", 0) + comp.get("L", 0)))

def disorder_proxy(seq: str) -> float:
    dis = set("PESQKAG")
    return sum(1 for aa in seq if aa in dis) / len(seq)

def compute_props(fasta_path: str, group: str) -> pd.DataFrame:
    rows = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = clean_seq(str(rec.seq))
        if len(seq) < 5:
            continue

        pa = ProteinAnalysis(seq)

        try:
            gravy = pa.gravy()
        except Exception:
            gravy = float("nan")
        try:
            ii = pa.instability_index()
        except Exception:
            ii = float("nan")
        try:
            arom = pa.aromaticity()
        except Exception:
            arom = float("nan")
        try:
            pi = pa.isoelectric_point()
        except Exception:
            pi = float("nan")
        try:
            mw = pa.molecular_weight()
        except Exception:
            mw = float("nan")

        rows.append({
            "group": group,
            "id": rec.id,
            "len": len(seq),
            "mw": mw,
            "pI": pi,
            "gravy": gravy,
            "aromaticity": arom,
            "instability_index": ii,
            "aliphatic_index": aliphatic_index(seq),
            "net_charge_pH7": net_charge_ph7(seq),
            "disorder_proxy": disorder_proxy(seq),
            "source_fasta": os.path.basename(fasta_path),
        })

    return pd.DataFrame(rows)

def discover_fastas(indir: str, include_genome: bool = True):
    """
    Discover deduplicated bestORF AA FASTAs.
    - foreground: lost/preserved (and optional favored/unfavoured)
    - control: genome-wide exon background from cache_bg_exons (if present)
    """
    pats = [
        r"^(lost|preserved)\.bestORF\.aa\.dedup\.fa$",
        r"^(lost|preserved)\.(favored|unfavoured)\.bestORF\.aa\.dedup\.fa$", # Corrected regex (removed space)
    ]

    out = []
    
    print(f"[INFO] Discovering fastas in directory: {indir}")
    files_in_indir = os.listdir(indir)
    print(f"[INFO] Files found: {files_in_indir}")

    # foreground (lost/preserved)
    for fn in files_in_indir:
        for p in pats:
            if re.match(p, fn):
                group = fn.replace(".bestORF.aa.dedup.fa", "")
                print(f"[INFO] Matched {fn} with pattern {p}. Group: {group}")
                out.append((group, os.path.join(indir, fn)))
                break

    # genome control from cache_bg_exons/
    if include_genome:
        bg_dir = os.path.join(indir, "cache_bg_exons")
        if os.path.isdir(bg_dir):
            bg_files = [f for f in os.listdir(bg_dir) if f.endswith(".exons.bestORF.aa.dedup.fa") and f.startswith("bgexons__")]
            for fn in sorted(bg_files):
                out.append(("genome", os.path.join(bg_dir, fn)))

    # deduplicate identical paths (just in case)
    seen = set()
    out2 = []
    for g, p in out:
        if p not in seen:
            out2.append((g, p))
            seen.add(p)

    return sorted(out2, key=lambda x: x[0])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--indir", required=True, help="Directory containing lost/preserved *.bestORF.aa.dedup.fa (and cache_bg_exons/ for genome control)")
    ap.add_argument("-o", "--out_tsv", default="biophys_properties.tsv")
    ap.add_argument("--no_genome", action="store_true", help="Do NOT include genome background from indir/cache_bg_exons/")
    args = ap.parse_args()

    fastas = discover_fastas(args.indir, include_genome=(not args.no_genome))
    if not fastas:
        raise SystemExit(f"No dedup bestORF fastas found in: {args.indir}")

    dfs = []
    for group, fp in fastas:
        if not os.path.exists(fp):
            continue
        print(f"[INFO] {group}: {fp}")
        df = compute_props(fp, group)
        dfs.append(df)

    if not dfs:
        raise SystemExit("No sequences parsed from discovered FASTAs (empty outputs).")

    all_df = pd.concat(dfs, ignore_index=True)
    all_df.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[OK] Wrote: {args.out_tsv} (n={len(all_df)})")

if __name__ == "__main__":
    main()
