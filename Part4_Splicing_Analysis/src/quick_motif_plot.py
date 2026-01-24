
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ame_tsv", required=True)
    parser.add_argument("--out_pdf", required=True)
    args = parser.parse_args()

    # Read AME TSV
    # Skip header comments (lines starting with #)
    try:
        df = pd.read_csv(args.ame_tsv, sep="\t", comment="#")
    except Exception as e:
        print(f"Error reading TSV: {e}")
        return

    if df.empty:
        print("No motifs found in results.")
        return

    # Filter for reasonable significance if needed, or just plot top N
    # Let's plot top 10 by p-value
    df = df.sort_values("p-value").head(20)

    # transform p-value
    df["neg_log_p"] = -np.log10(df["p-value"])
    
    # Use motif_alt_ID if available and readable, else motif_ID
    df["Motif"] = df.apply(lambda x: x["motif_alt_ID"] if pd.notna(x["motif_alt_ID"]) and x["motif_alt_ID"] != "." else x["motif_ID"], axis=1)

    plt.figure(figsize=(8, 6))
    sns.barplot(data=df, x="neg_log_p", y="Motif", palette="viridis")
    plt.xlabel("-log10(p-value)")
    plt.title("Motif Enrichment: Dependent vs Independent Introns")
    plt.tight_layout()
    plt.savefig(args.out_pdf)
    print(f"Plot saved to {args.out_pdf}")

if __name__ == "__main__":
    main()
