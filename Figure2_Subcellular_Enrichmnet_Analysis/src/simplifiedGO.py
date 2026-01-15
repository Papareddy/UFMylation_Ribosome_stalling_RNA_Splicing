import argparse
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from gprofiler import GProfiler

def run_go_analysis(gene_list, organism):
    """
    Runs GO enrichment analysis using gprofiler.
    """
    gp = GProfiler(return_dataframe=True)
    go_results = gp.profile(
        organism=organism,
        query=gene_list,
        sources=['GO:BP', 'GO:MF', 'GO:CC'],
        significance_threshold_method='fdr'
    )
    return go_results

def plot_go_terms(df, output_file, title, cluster_order):
    """
    Generates and saves a dot plot for GO terms for a single category.
    """
    if df.empty:
        print("DataFrame for plotting is empty. Skipping plot generation.")
        return

    # Filter for top N terms per cluster
    top_n = 10
    df_plot = df.groupby('Cluster').apply(lambda x: x.nsmallest(top_n, 'p_value')).reset_index(drop=True)

    if df_plot.empty:
        print("No terms left to plot after filtering.")
        return

    df_plot['short_name'] = df_plot['name'].str.wrap(60)
    df_plot['Cluster'] = pd.Categorical(df_plot['Cluster'], categories=cluster_order, ordered=True)
    
    # Sort y-axis terms by name for consistency
    df_plot = df_plot.sort_values('short_name')

    # Dynamically adjust height
    height = max(5, len(df_plot['short_name'].unique()) * 0.25)

    plt.figure(figsize=(10, height))
    ax = sns.scatterplot(
        data=df_plot,
        x="Cluster",
        y="short_name",
        hue="Fold_Enrichment",
        size="-log10(p_value)",
        palette="viridis_r",
        sizes=(50, 500)
    )

    ax.set_xlabel("Cluster")
    ax.set_ylabel("GO Term")
    ax.set_title(title)
    plt.xticks(rotation=45, ha='right')
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout(rect=[0, 0, 0.85, 1]) # Adjust layout to make space for legend

    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Plot saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Perform GO term analysis for all clusters and save results separated by GO category.")
    parser.add_argument("--input_file", default="Figure2_Subcellular_Enrichmnet_Analysis/data/Microsome_ANS_vsDMSO_Col0_ufm1_Limma_with_log2counts.tsv", help="Path to the input TSV file.")
    parser.add_argument("--output_dir", default="Figure2_Subcellular_Enrichmnet_Analysis/results", help="Directory to save results.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Read data
    data = pd.read_csv(args.input_file, sep='\t', engine='python')

    # 2. Get clusters in a specific order
    clusters_to_analyze = data['cluster'].dropna().unique()
    desired_order = ["Col0_Up", "Col0_Down", "ufm1_Up", "ufm1_Down"]
    clusters_in_order = [c for c in desired_order if c in clusters_to_analyze]

    all_go_results = []

    # 3. Loop through each cluster and perform GO analysis
    for cluster_name in clusters_in_order:
        print(f"--- Analyzing cluster: {cluster_name} ---")
        genes = data[data['cluster'] == cluster_name]['gene_id'].dropna().tolist()
        if not genes:
            continue
        go_results = run_go_analysis(genes, 'athaliana')
        if go_results is not None and not go_results.empty:
            go_results['Cluster'] = cluster_name
            all_go_results.append(go_results)

    if not all_go_results:
        print("No GO results found. Exiting.")
        return

    # 4. Combine all results and add calculated columns
    final_go_df = pd.concat(all_go_results, ignore_index=True)
    final_go_df['-log10(p_value)'] = -np.log10(final_go_df['p_value'])
    final_go_df['Fold_Enrichment'] = (final_go_df['intersection_size'] / final_go_df['query_size']) / (final_go_df['term_size'] / final_go_df['effective_domain_size'])

    # 5. Save results and plots SEPARATELY for each GO category
    for go_category, label in {"GO:BP": "Biological Process", "GO:MF": "Molecular Function", "GO:CC": "Cellular Component"}.items():
        print(f"--- Processing Category: {label} ---")
        
        category_df = final_go_df[final_go_df['source'] == go_category]
        
        if category_df.empty:
            print(f"No terms found for {label}.")
            continue

        # Save category-specific TSV
        output_tsv = os.path.join(args.output_dir, f"go_enrichment_{go_category.split(':')[1]}.tsv")
        category_df.to_csv(output_tsv, sep='\t', index=False)
        print(f"Saved results for {label} to {output_tsv}")

        # Generate category-specific plot
        plot_file = os.path.join(args.output_dir, f"go_plot_{go_category.split(':')[1]}.png")
        plot_go_terms(category_df, plot_file, f"GO Enrichment: {label}", clusters_in_order)


if __name__ == "__main__":
    main()