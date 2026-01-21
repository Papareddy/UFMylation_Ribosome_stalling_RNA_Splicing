#!/usr/bin/env python3
"""
Step 15B: Expression Boxplots by Gene Set

Compare WT (DMSO) vs UFM1 (DOX) expression levels across gene sets:
1. Background
2. All UFM1-dependent genes
3. All UFM1-independent genes
4. UFM1-dependent with SRSF3/PCBP2 motif (not overlapping with independent)

Separate plots for:
- Treatment: ANS, CTRL
- Fraction: NUCL, CYTO
"""

import argparse
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def parse_args():
    parser = argparse.ArgumentParser(description="Step 15B: Expression Boxplots by Gene Set")
    parser.add_argument("--dependent_tsv", required=True, help="Path to UFM1_dependent.tsv from Step 1")
    parser.add_argument("--independent_tsv", required=True, help="Path to UFM1_independent.tsv from Step 1")
    parser.add_argument("--motif_targets", required=True, help="Path to UFM1_dependent_SRSF3_PCBP2_targets.tsv from Step 15")
    parser.add_argument("--expression_file", required=True, help="Path to gene expression TPM file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    return parser.parse_args()


def load_gene_sets(dependent_tsv, independent_tsv, motif_targets_file):
    """Load gene sets for comparison."""
    
    # All UFM1-dependent genes
    dep_df = pd.read_csv(dependent_tsv, sep='\t')
    all_dep_genes = set(dep_df['GeneID'].str.replace(r'\.\d+$', '', regex=True).unique())
    print(f"[INFO] All UFM1-dependent genes: {len(all_dep_genes)}")
    
    # All UFM1-independent genes
    indep_df = pd.read_csv(independent_tsv, sep='\t')
    all_indep_genes = set(indep_df['GeneID'].str.replace(r'\.\d+$', '', regex=True).unique())
    print(f"[INFO] All UFM1-independent genes: {len(all_indep_genes)}")
    
    # UFM1-dependent with SRSF3/PCBP2 motif (not overlapping with independent)
    motif_df = pd.read_csv(motif_targets_file, sep='\t')
    # Filter for Dependent_Specific group (not Overlap_Independent)
    if 'group' in motif_df.columns:
        motif_dep_specific = set(motif_df[motif_df['group'] == 'Dependent_Specific']['gene_id'].unique())
    else:
        # Older format without group column - just get genes not in independent
        motif_genes = set(motif_df['gene_id'].unique())
        motif_dep_specific = motif_genes - all_indep_genes
    print(f"[INFO] UFM1-dependent with motif (not overlapping): {len(motif_dep_specific)}")
    
    return {
        'all_dep': all_dep_genes,
        'all_indep': all_indep_genes,
        'motif_dep_specific': motif_dep_specific
    }


def load_expression_data(expr_file):
    """Load and organize expression data by condition."""
    print(f"[INFO] Loading expression data: {expr_file}")
    df = pd.read_csv(expr_file, sep='\t')
    
    # Clean gene IDs
    gene_col = df.columns[0]
    df[gene_col] = df[gene_col].str.replace(r'\.\d+$', '', regex=True)
    df = df.rename(columns={gene_col: 'gene_id'})
    
    # Organize columns by condition
    # Pattern: [AAVS|UFM]N_FRAC_TREAT_repN
    # We need: UFM samples, FRAC (NUCL/CYTO), TREAT (DMSO/DOX), ANS/CTRL
    
    conditions = {}
    for col in df.columns[1:]:
        if 'AAVS' in col.upper():
            continue  # Skip AAVS samples
        if 'UFM' not in col.upper():
            continue
            
        # Parse column name
        # Example: UFM1_CYTO_DMSO_ANS_rep1, UFM1_NUCL_DOX_CTRL_rep2
        parts = col.upper().split('_')
        
        # Find fraction (NUCL/CYTO)
        frac = 'NUCL' if 'NUCL' in parts else 'CYTO' if 'CYTO' in parts else None
        # Find treatment (DMSO/DOX)
        treat = 'DMSO' if 'DMSO' in parts else 'DOX' if 'DOX' in parts else None
        # Find stress (ANS/CTRL)
        stress = 'ANS' if 'ANS' in parts else 'CTRL' if 'CTRL' in parts else None
        
        if frac and treat and stress:
            key = f"{frac}_{treat}_{stress}"
            if key not in conditions:
                conditions[key] = []
            conditions[key].append(col)
    
    # Calculate mean TPM for each condition
    result_df = df[['gene_id']].copy()
    for key, cols in conditions.items():
        if cols:
            result_df[key] = df[cols].mean(axis=1)
            print(f"[INFO] Condition {key}: {len(cols)} replicates")
    
    return result_df


def create_boxplot(data_dict, title, ylabel, outfile, gene_set_colors):
    """Create boxplot comparing groups with p-values."""
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Prepare data
    box_data = []
    box_labels = []
    box_colors = []
    positions = []
    
    pos = 0
    for gene_set, (dmso_vals, dox_vals, color) in data_dict.items():
        if len(dmso_vals) > 0:
            box_data.append(dmso_vals)
            box_labels.append(f'{gene_set}\nWT (DMSO)\n(n={len(dmso_vals)})')
            box_colors.append(color)
            positions.append(pos)
            pos += 1
        
        if len(dox_vals) > 0:
            box_data.append(dox_vals)
            box_labels.append(f'{gene_set}\nUFM1 (DOX)\n(n={len(dox_vals)})')
            box_colors.append(color)
            positions.append(pos)
            pos += 1
        
        pos += 0.5  # Gap between gene sets
    
    if len(box_data) == 0:
        plt.close(fig)
        return
    
    bp = ax.boxplot(box_data, positions=positions, patch_artist=True, widths=0.7)
    
    # Color boxes
    color_idx = 0
    for i, (patch, color) in enumerate(zip(bp['boxes'], box_colors)):
        patch.set_facecolor(color)
        patch.set_alpha(0.6 if i % 2 == 0 else 0.8)  # Lighter for DMSO, darker for DOX
    
    # Add p-values for DMSO vs DOX within each gene set
    y_max = max([np.percentile(d, 95) for d in box_data if len(d) > 0])
    y_range = y_max - min([np.percentile(d, 5) for d in box_data if len(d) > 0])
    
    gene_sets = list(data_dict.keys())
    stat_texts = []
    for idx, gene_set in enumerate(gene_sets):
        dmso_vals, dox_vals, _ = data_dict[gene_set]
        if len(dmso_vals) > 5 and len(dox_vals) > 5:
            _, p_val = stats.mannwhitneyu(dmso_vals, dox_vals, alternative='two-sided')
            stat_texts.append(f'{gene_set}: p={p_val:.2e}')
            
            # Draw significance bar
            x1 = idx * 2.5
            x2 = idx * 2.5 + 1
            y = y_max + (idx + 1) * 0.08 * y_range
            ax.plot([x1, x2], [y, y], 'k-', linewidth=1)
            sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
            ax.text((x1 + x2)/2, y + 0.01 * y_range, sig, ha='center', fontsize=9)
    
    ax.set_xticks(positions)
    ax.set_xticklabels(box_labels, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14)
    
    # Add p-value summary
    if stat_texts:
        ax.text(0.98, 0.02, '\n'.join(stat_texts), transform=ax.transAxes, fontsize=8,
               verticalalignment='bottom', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"[DONE] Saved {outfile}")


def create_combined_boxplot(all_data, outfile, gene_set_colors):
    """Create a 2x2 faceted boxplot figure."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    plt.subplots_adjust(hspace=0.4, wspace=0.2)
    
    # Define plot order
    plot_configs = [
        ('CYTO', 'CTRL', axes[0, 0]),
        ('CYTO', 'ANS', axes[0, 1]),
        ('NUCL', 'CTRL', axes[1, 0]),
        ('NUCL', 'ANS', axes[1, 1])
    ]
    
    for frac, stress, ax in plot_configs:
        key = f"{frac}_{stress}"
        if key not in all_data:
            ax.text(0.5, 0.5, f'No data for\n{frac} / {stress}', ha='center', va='center')
            continue
            
        data_dict = all_data[key]
        box_data = []
        box_labels = []
        box_colors = []
        positions = []
        
        pos = 0
        for gene_set, (dmso_vals, dox_vals, color) in data_dict.items():
            box_data.append(dmso_vals)
            box_labels.append(f'{gene_set}\nWT')
            box_colors.append(color)
            positions.append(pos)
            pos += 1
            
            box_data.append(dox_vals)
            box_labels.append(f'{gene_set}\nUFM1')
            box_colors.append(color)
            positions.append(pos)
            pos += 1
            
            pos += 0.5
            
        bp = ax.boxplot(box_data, positions=positions, patch_artist=True, widths=0.7)
        
        for i, (patch, color) in enumerate(zip(bp['boxes'], box_colors)):
            patch.set_facecolor(color)
            patch.set_alpha(0.6 if i % 2 == 0 else 0.8)
            
        y_max = max([np.percentile(d, 95) for d in box_data if len(d) > 0])
        y_min = min([np.percentile(d, 5) for d in box_data if len(d) > 0])
        y_range = y_max - y_min
        
        gene_sets = list(data_dict.keys())
        stat_texts = []
        for idx, gene_set in enumerate(gene_sets):
            dmso_vals, dox_vals, _ = data_dict[gene_set]
            if len(dmso_vals) > 5 and len(dox_vals) > 5:
                _, p_val = stats.mannwhitneyu(dmso_vals, dox_vals, alternative='two-sided')
                stat_texts.append(f'{gene_set}: p={p_val:.1e}')
                
                x1 = idx * 2.5
                x2 = idx * 2.5 + 1
                y = y_max + (idx + 1) * 0.08 * y_range
                ax.plot([x1, x2], [y, y], 'k-', linewidth=0.8)
                sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else 'ns'
                ax.text((x1 + x2)/2, y + 0.01 * y_range, sig, ha='center', fontsize=8)
        
        ax.set_xticks(positions)
        ax.set_xticklabels(box_labels, fontsize=7)
        ax.set_ylabel('log2(TPM + 1)')
        ax.set_title(f'{frac} / {stress}', fontsize=12, fontweight='bold')
        
    fig.suptitle('Expression Levels: WT (DMSO) vs UFM1 (DOX) across Conditions', fontsize=16, fontweight='bold')
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"[DONE] Saved combined plot to {outfile}")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    print("[INFO] === Step 15B: Expression Boxplots by Gene Set ===")
    
    # Load gene sets
    gene_sets = load_gene_sets(args.dependent_tsv, args.independent_tsv, args.motif_targets)
    
    # Load expression data
    expr_df = load_expression_data(args.expression_file)
    
    # Define colors for gene sets
    colors = {
        'Background': 'grey',
        'All UFM1-Dep': 'red',
        'All UFM1-Indep': 'blue',
        'Motif-Dep': 'darkgreen'
    }
    
    # Get background genes (not in any gene set)
    all_target_genes = gene_sets['all_dep'] | gene_sets['all_indep']
    
    # Generate plots for each combination
    fractions = ['NUCL', 'CYTO']
    stresses = ['ANS', 'CTRL']
    
    all_plot_data = {}
    
    for frac in fractions:
        for stress in stresses:
            print(f"[INFO] Processing {frac}_{stress}...")
            
            dmso_col = f"{frac}_DMSO_{stress}"
            dox_col = f"{frac}_DOX_{stress}"
            
            if dmso_col not in expr_df.columns or dox_col not in expr_df.columns:
                print(f"[WARN] Missing columns for {frac}_{stress}. Skipping.")
                continue
            
            # Prepare data for each gene set
            data_dict = {}
            
            # Background
            bg_mask = ~expr_df['gene_id'].isin(all_target_genes)
            bg_dmso = expr_df.loc[bg_mask, dmso_col].dropna().values
            bg_dox = expr_df.loc[bg_mask, dox_col].dropna().values
            data_dict['BG'] = (np.log2(bg_dmso + 1), np.log2(bg_dox + 1), colors['Background'])
            
            # All UFM1-dependent
            dep_mask = expr_df['gene_id'].isin(gene_sets['all_dep'])
            dep_dmso = expr_df.loc[dep_mask, dmso_col].dropna().values
            dep_dox = expr_df.loc[dep_mask, dox_col].dropna().values
            data_dict['Dep'] = (np.log2(dep_dmso + 1), np.log2(dep_dox + 1), colors['All UFM1-Dep'])
            
            # All UFM1-independent
            indep_mask = expr_df['gene_id'].isin(gene_sets['all_indep'])
            indep_dmso = expr_df.loc[indep_mask, dmso_col].dropna().values
            indep_dox = expr_df.loc[indep_mask, dox_col].dropna().values
            data_dict['Indep'] = (np.log2(indep_dmso + 1), np.log2(indep_dox + 1), colors['All UFM1-Indep'])
            
            # Motif-enriched dependent (not overlapping)
            motif_mask = expr_df['gene_id'].isin(gene_sets['motif_dep_specific'])
            motif_dmso = expr_df.loc[motif_mask, dmso_col].dropna().values
            motif_dox = expr_df.loc[motif_mask, dox_col].dropna().values
            data_dict['Motif'] = (np.log2(motif_dmso + 1), np.log2(motif_dox + 1), colors['Motif-Dep'])
            
            all_plot_data[f"{frac}_{stress}"] = data_dict
            
            # Create individual plot
            title = f'Expression: WT (DMSO) vs UFM1 (DOX) - {frac} / {stress}'
            ylabel = 'log2(TPM + 1)'
            outfile = os.path.join(args.outdir, f"Expression_Boxplot_{frac}_{stress}.pdf")
            create_boxplot(data_dict, title, ylabel, outfile, colors)
            
    # Create combined faceted plot
    combined_outfile = os.path.join(args.outdir, "Expression_Boxplots_Combined.png")
    create_combined_boxplot(all_plot_data, combined_outfile, colors)
    
    print("[INFO] === Step 15B Completed ====")


if __name__ == "__main__":
    main()
