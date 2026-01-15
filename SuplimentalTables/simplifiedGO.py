import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from gprofiler import GProfiler
import sys
import argparse

# Step 0: Parse arguments and show help if needed
parser = argparse.ArgumentParser(description="Perform GO enrichment analysis from input file with Entrez ID and Cluster.")
parser.add_argument("input_file", help="Input tab-separated file with two columns: 'GeneId' (Entrez ID) and 'Cluster'.")
args = parser.parse_args()

input_file = args.input_file
base_name = os.path.splitext(os.path.basename(input_file))[0]
output_dir = f"{base_name}_results"
os.makedirs(output_dir, exist_ok=True)
print(f"\n=== GO Enrichment Pipeline ===")
print(f"Input file: {input_file}")
print("Expecting two columns:")
print(" - GeneId: Entrez ID")
print(" - Cluster: cluster/group label")
print(f"Results will be saved in: {output_dir}\n")

# Step 1: Load data with cluster assignments
print(f"Loading data from '{input_file}'...")
data = pd.read_csv(input_file, sep='\s+')
print("Data loaded successfully. Saving a preview to the output directory...")
data.head().to_csv(os.path.join(output_dir, "data_preview.tsv"), sep='\t', index=False)
print("Data preview saved as 'data_preview.tsv'.")

# Step 2: Summarize cluster sizes
cluster_counts = data['Cluster'].value_counts().sort_index()
cluster_counts.to_csv(os.path.join(output_dir, "cluster_sizes.tsv"), sep='\t')
print("Cluster sizes saved as 'cluster_sizes.tsv'.")

# Plot cluster sizes
plt.figure(figsize=(10, 5))
sns.barplot(x=cluster_counts.index, y=cluster_counts.values, palette="tab10")
plt.xticks(rotation=45, ha='right')
plt.xlabel("Cluster")
plt.ylabel("Number of Genes")
plt.title("Number of Genes per Cluster")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "cluster_sizes.png"), dpi=300)
plt.show()
print("Cluster sizes plot saved as 'cluster_sizes.png'.")

# Step 3: GO term enrichment
print("Performing GO term enrichment using g:Profiler...")
gp = GProfiler(return_dataframe=True)
cluster_results = []

for cluster in data['Cluster'].unique():
    print(f"Analyzing cluster '{cluster}'...")
    genes_in_cluster = data[data['Cluster'] == cluster]['GeneId'].tolist()
    if not genes_in_cluster:
        continue
    result = gp.profile(
        organism='athaliana',
        query=genes_in_cluster,
        significance_threshold_method='fdr'
    )
    if not result.empty:
        result['Cluster'] = cluster
        cluster_results.append(result)
        print(f"Significant GO terms found for cluster '{cluster}'.")
    else:
        print(f"No significant GO terms found for cluster '{cluster}'.")

if cluster_results:
    go_df = pd.concat(cluster_results)
    go_df['-log10(p_value)'] = -np.log10(go_df['p_value'])
    go_df.to_csv(os.path.join(output_dir, "go_terms.tsv"), sep='\t', index=False)
    print("GO term enrichment results saved as 'go_terms.tsv'.")

    # Filter top GO terms and plot heatmap
    top_go_terms = go_df.groupby('Cluster').head(10)
    heatmap_data = top_go_terms.pivot(index="name", columns="Cluster", values="-log10(p_value)").fillna(0)
    plt.figure(figsize=(10, 8))
    sns.heatmap(heatmap_data, cmap="coolwarm", annot=False, cbar_kws={'label': '-log10(p-value)'})
    plt.title("GO Term Enrichment Heatmap")
    plt.xlabel("Cluster")
    plt.ylabel("GO Term Name")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "go_term_heatmap.png"), dpi=300)
    plt.show()
    print("GO Term Enrichment Heatmap saved as 'go_term_heatmap.png'.")

    # Split by GO categories and save separately
    for category, label in {
        "GO:BP": "bp",
        "GO:MF": "mf",
        "GO:CC": "cc"
    }.items():
        go_df_cat = go_df[go_df['source'] == category]
        if not go_df_cat.empty:
            out_file = os.path.join(output_dir, f"go_terms_{label}.tsv")
            go_df_cat.to_csv(out_file, sep='\t', index=False)
            print(f"{category} GO terms saved as '{out_file}'")

            # Plot category-specific heatmap
            top_cat = go_df_cat.groupby('Cluster').head(10)
            heatmap_data_cat = top_cat.pivot(index="name", columns="Cluster", values="-log10(p_value)").fillna(0)
            plt.figure(figsize=(10, 8))
            sns.heatmap(heatmap_data_cat, cmap="coolwarm", annot=False, cbar_kws={'label': '-log10(p-value)'})
            plt.title(f"{category} GO Term Enrichment Heatmap")
            plt.xlabel("Cluster")
            plt.ylabel("GO Term Name")
            plt.tight_layout()
            out_heatmap = os.path.join(output_dir, f"go_term_heatmap_{label}.png")
            plt.savefig(out_heatmap, dpi=300)
            plt.close()
            print(f"{category} heatmap saved as '{out_heatmap}'")
else:
    print("No significant GO terms found for any cluster.")
