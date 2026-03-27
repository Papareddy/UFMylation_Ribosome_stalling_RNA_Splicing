import pandas as pd
from Bio import SeqIO
import argparse
import os
from scipy.stats import mannwhitneyu

def get_stats(fa_path):
    lens = []
    gcs = []
    if not os.path.exists(fa_path): return [], []
    for r in SeqIO.parse(fa_path, "fasta"):
        s = str(r.seq).upper()
        if len(s) == 0: continue
        lens.append(len(s))
        gc = (s.count("G") + s.count("C")) / len(s)
        gcs.append(gc)
    return lens, gcs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lost", required=True)
    parser.add_argument("--preserved", required=True)
    args = parser.parse_args()
    
    l_len, l_gc = get_stats(args.lost)
    p_len, p_gc = get_stats(args.preserved)
    
    print(f"Stats for {os.path.basename(args.lost)} vs {os.path.basename(args.preserved)}")
    
    # Length
    u, p = mannwhitneyu(l_len, p_len)
    print(f"Length: Lost={sum(l_len)/len(l_len):.1f}bp, Pres={sum(p_len)/len(p_len):.1f}bp, p={p:.2e}")
    
    # GC
    u, p = mannwhitneyu(l_gc, p_gc)
    print(f"GC: Lost={sum(l_gc)/len(l_gc):.1%}, Pres={sum(p_gc)/len(p_gc):.1%}, p={p:.2e}")

if __name__ == "__main__":
    main()
