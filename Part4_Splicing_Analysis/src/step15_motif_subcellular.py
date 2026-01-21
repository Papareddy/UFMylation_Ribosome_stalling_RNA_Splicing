#!/usr/bin/env python3
"""
Step 15: Motif Target Identification & Subcellular Distribution Analysis

This script:
1. Identifies UFM1-dependent retained introns containing SRSF3_SRp20 or PCBP2 motif matches
2. Processes gene expression data to calculate Nucleus-to-Cytosol (N/C) ratios
3. Generates N/C shift scatter plots comparing CTRL vs ANS conditions
"""

import argparse
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_args():
    parser = argparse.ArgumentParser(description="Step 15: Motif Target Identification & Subcellular Distribution")
    parser.add_argument("--ame_dir", required=True, help="Directory containing AME results (e.g., step10_motif_analysis/RI_DeepDive_combined)")
    parser.add_argument("--dependent_tsv", required=True, help="Path to UFM1_dependent.tsv from Step 1")
    parser.add_argument("--independent_tsv", required=True, help="Path to UFM1_independent.tsv from Step 1 (for filtering)")
    parser.add_argument("--expression_file", required=True, help="Path to gene expression TPM file")
    parser.add_argument("--outdir", required=True, help="Output directory for Step 15 results")
    parser.add_argument("--intron_fasta", help="Path to UFM1_dependent intron FASTA (for direct motif scanning)")
    parser.add_argument("--motif_db", help="Path to MEME motif database")
    return parser.parse_args()


def load_ame_results(ame_dir, region="intron"):
    """
    Load AME results and identify genes with significant SRSF3/PCBP2 motif matches.
    """
    ame_file = os.path.join(ame_dir, f"ame_{region}_UFM1_dependent_vs_constitutive", "ame.tsv")
    
    if not os.path.exists(ame_file):
        print(f"[WARN] AME file not found: {ame_file}")
        return pd.DataFrame()
    
    # Read AME results (skip header lines starting with #)
    df = pd.read_csv(ame_file, sep='\t', comment='#')
    
    # Filter for SRSF3 (SRp20) and PCBP2 motifs
    srsf3_pattern = re.compile(r'SRSF3|SRp20|Serine.*arginine.*rich.*CHWCMC', re.IGNORECASE)
    pcbp2_pattern = re.compile(r'PCBP2', re.IGNORECASE)
    
    if 'motif_ID' in df.columns:
        motif_col = 'motif_ID'
    elif 'motif_id' in df.columns:
        motif_col = 'motif_id'
    else:
        print(f"[WARN] Could not find motif ID column in AME results. Columns: {df.columns.tolist()}")
        return pd.DataFrame()
    
    srsf3_hits = df[df[motif_col].str.contains(srsf3_pattern, na=False)]
    pcbp2_hits = df[df[motif_col].str.contains(pcbp2_pattern, na=False)]
    
    print(f"[INFO] Found {len(srsf3_hits)} SRSF3/SRp20 motif hits")
    print(f"[INFO] Found {len(pcbp2_hits)} PCBP2 motif hits")
    
    return pd.concat([srsf3_hits, pcbp2_hits]).drop_duplicates()


def scan_fasta_for_motifs(fasta_file, dependent_df):
    """
    Scan intron sequences for SRSF3_SRp20 (CHWCMC) and PCBP2 (CCUYCCC, UWCCC) motifs.
    Returns a DataFrame with gene_id, dpsi, and motif_id for each match.
    """
    try:
        from Bio import SeqIO
    except ImportError:
        print("[WARN] BioPython not available. Skipping FASTA motif scanning.")
        return pd.DataFrame()
    
    # IUPAC codes: H=[ACT], W=[AT], Y=[CT], M=[AC], C=C, U=T (in DNA)
    # SRSF3_SRp20: CHWCMC -> C[ACT][AT]C[AC]C
    srsf3_regex = re.compile(r'C[ACT][AT]C[AC]C', re.IGNORECASE)
    # PCBP2 variants:
    # CCUYCCC -> CCT[CT]CCC
    # UWCCC -> [AT][AT]CCC
    pcbp2_regex1 = re.compile(r'CCT[CT]CCC', re.IGNORECASE)
    pcbp2_regex2 = re.compile(r'[AT][AT]CCC', re.IGNORECASE)
    
    results = []
    
    if not os.path.exists(fasta_file):
        print(f"[WARN] FASTA file not found: {fasta_file}")
        return pd.DataFrame()
    
    # Build mappings from FASTA ID to gene info
    # FASTA IDs are in format: ID=N::chr:start-end(strand)
    # where N is the 0-based row index in the dependent TSV
    
    # Map by row index (ID=N)
    index_map = {}
    # Map by coordinates (chr:start-end)
    coord_map = {}
    
    dpsi_col = 'IncLevelDifference.WT' if 'IncLevelDifference.WT' in dependent_df.columns else 'IncLevelDifference'
    chr_col = 'chr.WT' if 'chr.WT' in dependent_df.columns else 'chr'
    start_col = 'riExonStart_0base.WT' if 'riExonStart_0base.WT' in dependent_df.columns else 'riExonStart_0base'
    end_col = 'riExonEnd.WT' if 'riExonEnd.WT' in dependent_df.columns else 'riExonEnd'
    
    if 'GeneID' in dependent_df.columns:
        for idx, row in dependent_df.iterrows():
            gene_id = str(row['GeneID']).split('.')[0]  # Remove version
            dpsi = row.get(dpsi_col, np.nan)
            
            # Index-based key
            index_map[idx] = (gene_id, dpsi)
            
            # Coordinate-based key
            chrom = str(row.get(chr_col, '')).replace('chr', '')
            start = row.get(start_col, '')
            end = row.get(end_col, '')
            coord_key = f"{chrom}:{start}-{end}"
            coord_map[coord_key] = (gene_id, dpsi)
    
    print(f"[INFO] Built mappings for {len(index_map)} events by index, {len(coord_map)} by coordinates")
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        seq_id = record.id
        
        # Try to extract gene info from sequence ID
        gene_id = None
        dpsi = np.nan
        
        # Parse FASTA ID: ID=N::chr:start-end(strand)
        # Extract row index N
        id_match = re.match(r'ID=(\d+)', seq_id)
        if id_match:
            row_idx = int(id_match.group(1))
            if row_idx in index_map:
                gene_id, dpsi = index_map[row_idx]
        
        # If not found, try coordinate matching
        if gene_id is None:
            coord_match = re.search(r'::([^:]+):(\d+)-(\d+)', seq_id)
            if coord_match:
                chrom = coord_match.group(1).replace('chr', '')
                start = coord_match.group(2)
                end = coord_match.group(3)
                coord_key = f"{chrom}:{start}-{end}"
                if coord_key in coord_map:
                    gene_id, dpsi = coord_map[coord_key]
        
        if gene_id is None:
            continue
        
        # Scan for motifs
        if srsf3_regex.search(seq):
            results.append({'gene_id': gene_id, 'dpsi': dpsi, 'motif_id': 'SRSF3_SRp20'})
        
        if pcbp2_regex1.search(seq):
            results.append({'gene_id': gene_id, 'dpsi': dpsi, 'motif_id': 'PCBP2_CCUYCCC'})
            
        if pcbp2_regex2.search(seq):
            results.append({'gene_id': gene_id, 'dpsi': dpsi, 'motif_id': 'PCBP2_UWCCC'})
    
    df = pd.DataFrame(results)
    if len(df) > 0:
        df = df.drop_duplicates()
        print(f"[INFO] Found {len(df)} motif matches in {len(df['gene_id'].unique())} unique genes")
    return df


def process_expression_data(expr_file):
    """
    Process gene expression data:
    1. Clean gene IDs (remove version numbers)
    2. Filter columns (keep UFM, discard AAVS)
    3. Normalize column names
    4. Calculate mean TPM per condition
    """
    print(f"[INFO] Loading expression data: {expr_file}")
    df = pd.read_csv(expr_file, sep='\t')
    
    # Identify gene ID column
    gene_col = df.columns[0]
    print(f"[INFO] Gene ID column: {gene_col}")
    
    # Clean gene IDs (remove version numbers)
    df[gene_col] = df[gene_col].str.replace(r'\.\d+$', '', regex=True)
    
    # Filter columns: keep gene_id and UFM columns, discard AAVS
    keep_cols = [gene_col]
    for col in df.columns[1:]:
        if 'AAVS' in col.upper():
            continue
        if 'UFM' in col.upper():
            keep_cols.append(col)
    
    df = df[keep_cols]
    print(f"[INFO] Kept {len(keep_cols)-1} expression columns after filtering")
    
    # Normalize column names: Remove digit after UFM (UFM1_, UFM2_, UFM3_ -> UFM_)
    new_cols = {gene_col: 'gene_id'}
    for col in df.columns[1:]:
        new_col = re.sub(r'UFM\d+_', 'UFM_', col)
        new_cols[col] = new_col
    df = df.rename(columns=new_cols)
    
    # Group replicates and calculate mean TPM
    # Expected pattern: UFM_<condition>_<fraction>_rep<N>
    # e.g., UFM_DMSO_NUCL_rep1, UFM_ANS_CYTO_rep2
    
    conditions = {}
    for col in df.columns[1:]:
        # Extract base condition (remove _repN suffix)
        base = re.sub(r'_rep\d+$', '', col)
        if base not in conditions:
            conditions[base] = []
        conditions[base].append(col)
    
    # Calculate mean for each condition
    mean_df = df[['gene_id']].copy()
    for cond, cols in conditions.items():
        if len(cols) > 0:
            mean_df[f'Mean_{cond}'] = df[cols].mean(axis=1)
    
    print(f"[INFO] Calculated mean TPM for {len(conditions)} conditions")
    return mean_df


def calculate_nc_ratios(expr_df):
    """
    Calculate Nucleus-to-Cytosol (N/C) ratios.
    Formula: (Mean_NUCL + 1) / (Mean_CYTO + 1)
    
    Calculate for:
    - DMSO context: ANS ratio, CTRL ratio
    - DOX context: ANS ratio, CTRL ratio
    """
    result_df = expr_df[['gene_id']].copy()
    
    # Identify available conditions
    nucl_cols = [c for c in expr_df.columns if 'NUCL' in c.upper()]
    cyto_cols = [c for c in expr_df.columns if 'CYTO' in c.upper()]
    
    print(f"[INFO] Nuclear columns: {nucl_cols}")
    print(f"[INFO] Cytosolic columns: {cyto_cols}")
    
    # Define condition pairs
    # DMSO context
    dmso_ans_nucl = [c for c in nucl_cols if 'DMSO' in c.upper() and 'ANS' in c.upper()]
    dmso_ans_cyto = [c for c in cyto_cols if 'DMSO' in c.upper() and 'ANS' in c.upper()]
    dmso_ctrl_nucl = [c for c in nucl_cols if 'DMSO' in c.upper() and ('CTRL' in c.upper() or all(x not in c.upper() for x in ['ANS', 'DOX']))]
    dmso_ctrl_cyto = [c for c in cyto_cols if 'DMSO' in c.upper() and ('CTRL' in c.upper() or all(x not in c.upper() for x in ['ANS', 'DOX']))]
    
    # DOX context
    dox_ans_nucl = [c for c in nucl_cols if 'DOX' in c.upper() and 'ANS' in c.upper()]
    dox_ans_cyto = [c for c in cyto_cols if 'DOX' in c.upper() and 'ANS' in c.upper()]
    dox_ctrl_nucl = [c for c in nucl_cols if 'DOX' in c.upper() and ('CTRL' in c.upper() or all(x not in c.upper() for x in ['ANS']))]
    dox_ctrl_cyto = [c for c in cyto_cols if 'DOX' in c.upper() and ('CTRL' in c.upper() or all(x not in c.upper() for x in ['ANS']))]
    
    def calc_ratio(nucl_cols, cyto_cols, name):
        if nucl_cols and cyto_cols:
            nucl_val = expr_df[nucl_cols[0]] if len(nucl_cols) == 1 else expr_df[nucl_cols].mean(axis=1)
            cyto_val = expr_df[cyto_cols[0]] if len(cyto_cols) == 1 else expr_df[cyto_cols].mean(axis=1)
            result_df[name] = (nucl_val + 1) / (cyto_val + 1)
            print(f"[INFO] Calculated {name}")
    
    calc_ratio(dmso_ans_nucl if dmso_ans_nucl else [c for c in nucl_cols if 'ANS' in c.upper()][:1],
               dmso_ans_cyto if dmso_ans_cyto else [c for c in cyto_cols if 'ANS' in c.upper()][:1],
               'NC_DMSO_ANS')
    calc_ratio(dmso_ctrl_nucl if dmso_ctrl_nucl else [c for c in nucl_cols if 'CTRL' in c.upper() or 'DMSO' in c.upper()][:1],
               dmso_ctrl_cyto if dmso_ctrl_cyto else [c for c in cyto_cols if 'CTRL' in c.upper() or 'DMSO' in c.upper()][:1],
               'NC_DMSO_CTRL')
    calc_ratio(dox_ans_nucl if dox_ans_nucl else [], dox_ans_cyto if dox_ans_cyto else [], 'NC_DOX_ANS')
    calc_ratio(dox_ctrl_nucl if dox_ctrl_nucl else [], dox_ctrl_cyto if dox_ctrl_cyto else [], 'NC_DOX_CTRL')
    
    return result_df


def plot_nc_shift(merged_df, all_genes_df, motif, context, outdir):
    """
    Generate N/C shift scatter plot for a specific motif and context.
    X-axis: N/C Ratio (CTRL)
    Y-axis: N/C Ratio (ANS)
    """
    ctrl_col = f'NC_{context}_CTRL'
    ans_col = f'NC_{context}_ANS'
    
    if ctrl_col not in all_genes_df.columns or ans_col not in all_genes_df.columns:
        print(f"[WARN] Missing columns for {context}: {ctrl_col}, {ans_col}")
        return None
    
    fig, ax = plt.subplots(figsize=(6, 6))
    
    # Plot all genes as grey background
    all_valid = all_genes_df[[ctrl_col, ans_col]].dropna()
    ax.scatter(np.log2(all_valid[ctrl_col]), np.log2(all_valid[ans_col]), 
               c='grey', alpha=0.1, s=5, label='All Genes')
    
    # Plot motif-containing genes
    motif_genes = merged_df[merged_df['motif_id'] == motif]
    if len(motif_genes) > 0:
        motif_valid = motif_genes[[ctrl_col, ans_col, 'gene_id']].dropna()
        if len(motif_valid) > 0:
            color = 'red' if 'SRSF3' in motif else 'blue'
            ax.scatter(np.log2(motif_valid[ctrl_col]), np.log2(motif_valid[ans_col]),
                      c=color, alpha=0.7, s=30, label=f'{motif} ({len(motif_valid)} genes)')
    
    # Diagonal line
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1)
    
    ax.set_xlabel(f'log2(N/C Ratio, {context} CTRL)')
    ax.set_ylabel(f'log2(N/C Ratio, {context} ANS)')
    ax.set_title(f'{motif} - {context} Context')
    ax.legend(loc='upper left', fontsize=8)
    ax.set_aspect('equal', adjustable='box')
    
    return fig


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    print("[INFO] === Step 15: Motif Target Identification & Subcellular Distribution ===")
    
    # 1. Load dependent events
    print(f"[INFO] Loading dependent events: {args.dependent_tsv}")
    dependent_df = pd.read_csv(args.dependent_tsv, sep='\t')
    print(f"[INFO] Loaded {len(dependent_df)} UFM1-dependent events")
    
    # 1b. Load independent events and extract gene IDs for filtering
    print(f"[INFO] Loading independent events: {args.independent_tsv}")
    independent_df = pd.read_csv(args.independent_tsv, sep='\t')
    independent_genes = set(independent_df['GeneID'].str.replace(r'\.\d+$', '', regex=True).unique())
    print(f"[INFO] Loaded {len(independent_df)} UFM1-independent events ({len(independent_genes)} unique genes to exclude)")
    
    # 2. Identify motif-containing genes
    motif_targets = pd.DataFrame()
    
    # Try AME results first
    if args.ame_dir and os.path.exists(args.ame_dir):
        ame_results = load_ame_results(args.ame_dir)
        if len(ame_results) > 0:
            print(f"[INFO] Loaded {len(ame_results)} motif hits from AME")
    
    # Scan FASTA for direct motif matches
    if args.intron_fasta and os.path.exists(args.intron_fasta):
        fasta_targets = scan_fasta_for_motifs(args.intron_fasta, dependent_df)
        if len(fasta_targets) > 0:
            motif_targets = pd.concat([motif_targets, fasta_targets]).drop_duplicates()
            print(f"[INFO] Found {len(fasta_targets)} motif matches in FASTA")
    
    # If no motif scanning available, create targets from dependent events with placeholder
    if len(motif_targets) == 0:
        print("[WARN] No motif scanning results. Creating targets from dependent events...")
        # Find the correct column names
        dpsi_col = None
        for col in ['IncLevelDifference', 'IncLevelDifference.WT']:
            if col in dependent_df.columns:
                dpsi_col = col
                break
        
        if dpsi_col:
            motif_targets = dependent_df[['GeneID', dpsi_col]].copy()
            motif_targets.columns = ['gene_id', 'dpsi']
        else:
            motif_targets = dependent_df[['GeneID']].copy()
            motif_targets.columns = ['gene_id']
            motif_targets['dpsi'] = np.nan
        
        motif_targets['gene_id'] = motif_targets['gene_id'].str.replace(r'\.\d+$', '', regex=True)
        motif_targets['motif_id'] = 'ALL_DEPENDENT'
    
    # 2b. Separate genes into UFM1-dependent-specific vs overlap with UFM1-independent
    before_filter = len(motif_targets)
    before_genes = len(motif_targets['gene_id'].unique())
    
    # Track which genes are excluded (overlap with independent)
    excluded_genes_mask = motif_targets['gene_id'].isin(independent_genes)
    excluded_targets = motif_targets[excluded_genes_mask].copy()
    excluded_targets['group'] = 'Overlap_Independent'
    
    # Keep only UFM1-dependent-specific
    motif_targets = motif_targets[~excluded_genes_mask].copy()
    motif_targets['group'] = 'Dependent_Specific'
    
    after_filter = len(motif_targets)
    after_genes = len(motif_targets['gene_id'].unique()) if len(motif_targets) > 0 else 0
    excluded_count = len(excluded_targets['gene_id'].unique()) if len(excluded_targets) > 0 else 0
    print(f"[INFO] Filtered: {before_genes} genes -> {after_genes} UFM1-dependent-specific + {excluded_count} overlap with independent")
    
    # Save both files
    targets_file = os.path.join(args.outdir, "UFM1_dependent_SRSF3_PCBP2_targets.tsv")
    motif_targets.to_csv(targets_file, sep='\t', index=False)
    print(f"[DONE] Saved motif targets (dependent-specific only) to {targets_file}")
    
    excluded_file = os.path.join(args.outdir, "UFM1_overlap_independent_targets.tsv")
    excluded_targets.to_csv(excluded_file, sep='\t', index=False)
    print(f"[DONE] Saved overlap genes to {excluded_file}")
    
    
    # 3. Process expression data
    expr_df = process_expression_data(args.expression_file)
    
    # 4. Calculate N/C ratios
    nc_df = calculate_nc_ratios(expr_df)
    
    # Save N/C ratios
    nc_file = os.path.join(args.outdir, "NC_ratios.tsv")
    nc_df.to_csv(nc_file, sep='\t', index=False)
    print(f"[DONE] Saved N/C ratios to {nc_file}")
    
    # 5. Merge with motif targets
    merged_df = motif_targets.merge(nc_df, on='gene_id', how='left')
    merged_file = os.path.join(args.outdir, "motif_targets_with_NC_ratios.tsv")
    merged_df.to_csv(merged_file, sep='\t', index=False)
    print(f"[DONE] Saved merged data to {merged_file}")
    
    # 6. Generate plots
    print("[INFO] Generating N/C shift scatter plots...")
    
    motifs = motif_targets['motif_id'].unique()
    contexts = ['DMSO', 'DOX']
    
    # Merge excluded targets with NC ratios too
    excluded_merged = excluded_targets.merge(nc_df, on='gene_id', how='left') if len(excluded_targets) > 0 else pd.DataFrame()
    
    pdf_file = os.path.join(args.outdir, "NC_shift_plots.pdf")
    with PdfPages(pdf_file) as pdf:
        for context in contexts:
            for motif in motifs:
                fig = plot_nc_shift(merged_df, nc_df, motif, context, args.outdir)
                if fig:
                    pdf.savefig(fig, bbox_inches='tight')
                    plt.close(fig)
    
    print(f"[DONE] Saved plots to {pdf_file}")
    
    # Also save combined multi-panel figure showing both groups
    fig_combined, axes = plt.subplots(2, 2, figsize=(12, 12))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    
    for i, context in enumerate(contexts):
        for j, motif_type in enumerate(['SRSF3', 'PCBP2']):
            ax = axes[i, j]
            ctrl_col = f'NC_{context}_CTRL'
            ans_col = f'NC_{context}_ANS'
            
            if ctrl_col in nc_df.columns and ans_col in nc_df.columns:
                # Plot all genes as grey background
                all_valid = nc_df[[ctrl_col, ans_col]].dropna()
                ax.scatter(np.log2(all_valid[ctrl_col]), np.log2(all_valid[ans_col]),
                          c='grey', alpha=0.1, s=5, label='All genes')
                
                # Plot excluded/overlap genes (lighter color)
                if len(excluded_merged) > 0:
                    excl_genes = excluded_merged[excluded_merged['motif_id'].str.contains(motif_type, case=False, na=False)]
                    if len(excl_genes) > 0:
                        excl_valid = excl_genes[[ctrl_col, ans_col]].dropna()
                        if len(excl_valid) > 0:
                            excl_color = 'orange' if motif_type == 'SRSF3' else 'cyan'
                            ax.scatter(np.log2(excl_valid[ctrl_col]), np.log2(excl_valid[ans_col]),
                                      c=excl_color, alpha=0.6, s=25, marker='o',
                                      label=f'{motif_type} Overlap (n={len(excl_valid)})')
                
                # Plot UFM1-dependent-specific genes (darker color, on top)
                motif_genes = merged_df[merged_df['motif_id'].str.contains(motif_type, case=False, na=False)]
                if len(motif_genes) > 0:
                    motif_valid = motif_genes[[ctrl_col, ans_col]].dropna()
                    if len(motif_valid) > 0:
                        dep_color = 'darkred' if motif_type == 'SRSF3' else 'darkblue'
                        ax.scatter(np.log2(motif_valid[ctrl_col]), np.log2(motif_valid[ans_col]),
                                  c=dep_color, alpha=0.8, s=35, marker='s',
                                  label=f'{motif_type} Specific (n={len(motif_valid)})')
                
                lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
                ax.plot(lims, lims, 'k--', alpha=0.5)
                ax.set_xlabel(f'log2(N/C, {context} CTRL)')
                ax.set_ylabel(f'log2(N/C, {context} ANS)')
                ax.set_title(f'{motif_type} - {context}')
                ax.legend(loc='upper left', fontsize=7)
            else:
                ax.text(0.5, 0.5, f'No data for\n{context}', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{motif_type} - {context}')
    
    fig_combined.suptitle('N/C Ratio Shift: UFM1-Dependent-Specific vs Overlap with Independent', fontsize=12, fontweight='bold')
    combined_file = os.path.join(args.outdir, "NC_shift_combined.png")
    fig_combined.savefig(combined_file, dpi=150, bbox_inches='tight')
    plt.close(fig_combined)
    print(f"[DONE] Saved combined plot to {combined_file}")
    
    # 7. CDF Plot and Boxplot with Shift Metric
    print("[INFO] Generating CDF and Boxplot for Shift Metric...")
    from scipy import stats
    
    for context in contexts:
        ctrl_col = f'NC_{context}_CTRL'
        ans_col = f'NC_{context}_ANS'
        
        if ctrl_col not in nc_df.columns or ans_col not in nc_df.columns:
            continue
        
        # Calculate Shift Metric: log2(NC_ANS / NC_CTRL)
        nc_df_valid = nc_df[[ctrl_col, ans_col, 'gene_id']].dropna()
        nc_df_valid['Shift_Metric'] = np.log2(nc_df_valid[ans_col] / nc_df_valid[ctrl_col])
        
        # Prepare groups
        # Background: all genes not in dependent or independent
        all_dep_genes = set(motif_targets['gene_id'].unique())
        all_excl_genes = set(excluded_targets['gene_id'].unique()) if len(excluded_targets) > 0 else set()
        
        bg_mask = ~nc_df_valid['gene_id'].isin(all_dep_genes | all_excl_genes)
        dep_mask = nc_df_valid['gene_id'].isin(all_dep_genes)
        excl_mask = nc_df_valid['gene_id'].isin(all_excl_genes)
        
        bg_shift = nc_df_valid.loc[bg_mask, 'Shift_Metric'].dropna()
        dep_shift = nc_df_valid.loc[dep_mask, 'Shift_Metric'].dropna()
        excl_shift = nc_df_valid.loc[excl_mask, 'Shift_Metric'].dropna()
        
        # A. CDF Plot
        fig_cdf, ax_cdf = plt.subplots(figsize=(8, 6))
        
        if len(bg_shift) > 0:
            sorted_bg = np.sort(bg_shift)
            cdf_bg = np.arange(1, len(sorted_bg) + 1) / len(sorted_bg)
            ax_cdf.plot(sorted_bg, cdf_bg, color='grey', linewidth=2, label=f'Background (n={len(bg_shift)})')
        
        if len(dep_shift) > 0:
            sorted_dep = np.sort(dep_shift)
            cdf_dep = np.arange(1, len(sorted_dep) + 1) / len(sorted_dep)
            ax_cdf.plot(sorted_dep, cdf_dep, color='red', linewidth=2, label=f'UFM1-Dependent-Specific (n={len(dep_shift)})')
        
        if len(excl_shift) > 0:
            sorted_excl = np.sort(excl_shift)
            cdf_excl = np.arange(1, len(sorted_excl) + 1) / len(sorted_excl)
            ax_cdf.plot(sorted_excl, cdf_excl, color='blue', linewidth=2, label=f'Overlap with Independent (n={len(excl_shift)})')
        
        ax_cdf.axvline(x=0, color='black', linestyle='--', alpha=0.5)
        ax_cdf.set_xlabel('Shift Metric (log2[N/C ANS / N/C CTRL])', fontsize=12)
        ax_cdf.set_ylabel('Cumulative Probability', fontsize=12)
        ax_cdf.set_title(f'N/C Ratio Shift Distribution - {context}', fontsize=14)
        ax_cdf.legend(loc='lower right', fontsize=10)
        ax_cdf.grid(True, alpha=0.3)
        
        cdf_file = os.path.join(args.outdir, f"NC_shift_CDF_{context}.pdf")
        fig_cdf.savefig(cdf_file, dpi=150, bbox_inches='tight')
        plt.close(fig_cdf)
        print(f"[DONE] Saved CDF plot to {cdf_file}")
        
        # B. Boxplot with Statistics
        fig_box, ax_box = plt.subplots(figsize=(8, 6))
        
        # Prepare data for boxplot
        box_data = []
        box_labels = []
        box_colors = []
        
        if len(bg_shift) > 0:
            box_data.append(bg_shift.values)
            box_labels.append(f'Background\n(n={len(bg_shift)})')
            box_colors.append('grey')
        
        if len(dep_shift) > 0:
            box_data.append(dep_shift.values)
            box_labels.append(f'Dependent\nSpecific\n(n={len(dep_shift)})')
            box_colors.append('red')
        
        if len(excl_shift) > 0:
            box_data.append(excl_shift.values)
            box_labels.append(f'Overlap\nIndep.\n(n={len(excl_shift)})')
            box_colors.append('blue')
        
        bp = ax_box.boxplot(box_data, labels=box_labels, patch_artist=True)
        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        
        # Add statistical tests
        stat_results = []
        y_max = max([np.percentile(d, 95) for d in box_data if len(d) > 0])
        y_offset = 0.15 * (y_max - np.min([np.percentile(d, 5) for d in box_data if len(d) > 0]))
        
        if len(bg_shift) > 0 and len(dep_shift) > 0:
            stat_bg_dep, p_bg_dep = stats.mannwhitneyu(bg_shift, dep_shift, alternative='two-sided')
            stat_results.append(f'BG vs Dep: p={p_bg_dep:.2e}')
            # Add significance bar
            ax_box.plot([1, 2], [y_max + y_offset, y_max + y_offset], 'k-', linewidth=1)
            sig_text = '***' if p_bg_dep < 0.001 else '**' if p_bg_dep < 0.01 else '*' if p_bg_dep < 0.05 else 'ns'
            ax_box.text(1.5, y_max + y_offset + 0.02, sig_text, ha='center', fontsize=10)
        
        if len(bg_shift) > 0 and len(excl_shift) > 0:
            stat_bg_excl, p_bg_excl = stats.mannwhitneyu(bg_shift, excl_shift, alternative='two-sided')
            stat_results.append(f'BG vs Overlap: p={p_bg_excl:.2e}')
            # Add significance bar
            ax_box.plot([1, 3], [y_max + 2*y_offset, y_max + 2*y_offset], 'k-', linewidth=1)
            sig_text = '***' if p_bg_excl < 0.001 else '**' if p_bg_excl < 0.01 else '*' if p_bg_excl < 0.05 else 'ns'
            ax_box.text(2, y_max + 2*y_offset + 0.02, sig_text, ha='center', fontsize=10)
        
        if len(dep_shift) > 0 and len(excl_shift) > 0:
            stat_dep_excl, p_dep_excl = stats.mannwhitneyu(dep_shift, excl_shift, alternative='two-sided')
            stat_results.append(f'Dep vs Overlap: p={p_dep_excl:.2e}')
        
        ax_box.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        ax_box.set_ylabel('Shift Metric (log2[N/C ANS / N/C CTRL])', fontsize=12)
        ax_box.set_title(f'N/C Ratio Shift Comparison - {context}', fontsize=14)
        
        # Add p-values as text
        stat_text = '\n'.join(stat_results)
        ax_box.text(0.98, 0.02, stat_text, transform=ax_box.transAxes, fontsize=9,
                   verticalalignment='bottom', horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        box_file = os.path.join(args.outdir, f"NC_shift_Boxplot_{context}.pdf")
        fig_box.savefig(box_file, dpi=150, bbox_inches='tight')
        plt.close(fig_box)
        print(f"[DONE] Saved Boxplot to {box_file}")
    
    # 8. Additional Plots: All UFM1-dependent genes and Motif-specific plots
    print("[INFO] Generating additional CDF/Boxplot analyses...")
    
    # Get all UFM1-dependent genes (from the original dependent_df, not just motif-containing)
    all_dep_gene_ids = set(dependent_df['GeneID'].str.replace(r'\.\d+$', '', regex=True).unique())
    all_indep_gene_ids = set(independent_genes)  # Already cleaned
    
    for context in contexts:
        ctrl_col = f'NC_{context}_CTRL'
        ans_col = f'NC_{context}_ANS'
        
        if ctrl_col not in nc_df.columns or ans_col not in nc_df.columns:
            continue
        
        # Calculate Shift Metric
        nc_df_valid = nc_df[[ctrl_col, ans_col, 'gene_id']].dropna().copy()
        nc_df_valid['Shift_Metric'] = np.log2(nc_df_valid[ans_col] / nc_df_valid[ctrl_col])
        
        # --- PLOT A: All UFM1-Dependent vs Background ---
        print(f"[INFO] Creating All UFM1-Dependent plot for {context}...")
        
        # Masks for all dependent genes
        all_dep_mask = nc_df_valid['gene_id'].isin(all_dep_gene_ids)
        all_indep_mask = nc_df_valid['gene_id'].isin(all_indep_gene_ids)
        bg_all_mask = ~(all_dep_mask | all_indep_mask)
        
        bg_all_shift = nc_df_valid.loc[bg_all_mask, 'Shift_Metric'].dropna()
        all_dep_shift = nc_df_valid.loc[all_dep_mask, 'Shift_Metric'].dropna()
        all_indep_shift = nc_df_valid.loc[all_indep_mask, 'Shift_Metric'].dropna()
        
        # CDF for all dependent
        fig_cdf_all, ax_cdf_all = plt.subplots(figsize=(8, 6))
        
        if len(bg_all_shift) > 0:
            sorted_bg = np.sort(bg_all_shift)
            cdf_bg = np.arange(1, len(sorted_bg) + 1) / len(sorted_bg)
            ax_cdf_all.plot(sorted_bg, cdf_bg, color='grey', linewidth=2, label=f'Background (n={len(bg_all_shift)})')
        
        if len(all_dep_shift) > 0:
            sorted_dep = np.sort(all_dep_shift)
            cdf_dep = np.arange(1, len(sorted_dep) + 1) / len(sorted_dep)
            ax_cdf_all.plot(sorted_dep, cdf_dep, color='red', linewidth=2, label=f'All UFM1-Dependent (n={len(all_dep_shift)})')
        
        if len(all_indep_shift) > 0:
            sorted_indep = np.sort(all_indep_shift)
            cdf_indep = np.arange(1, len(sorted_indep) + 1) / len(sorted_indep)
            ax_cdf_all.plot(sorted_indep, cdf_indep, color='blue', linewidth=2, label=f'All UFM1-Independent (n={len(all_indep_shift)})')
        
        ax_cdf_all.axvline(x=0, color='black', linestyle='--', alpha=0.5)
        ax_cdf_all.set_xlabel('Shift Metric (log2[N/C ANS / N/C CTRL])', fontsize=12)
        ax_cdf_all.set_ylabel('Cumulative Probability', fontsize=12)
        ax_cdf_all.set_title(f'All UFM1-Dependent: N/C Shift Distribution - {context}', fontsize=14)
        ax_cdf_all.legend(loc='lower right', fontsize=10)
        ax_cdf_all.grid(True, alpha=0.3)
        
        cdf_all_file = os.path.join(args.outdir, f"NC_shift_CDF_AllDependent_{context}.pdf")
        fig_cdf_all.savefig(cdf_all_file, dpi=150, bbox_inches='tight')
        plt.close(fig_cdf_all)
        print(f"[DONE] Saved CDF plot (All Dependent) to {cdf_all_file}")
        
        # Boxplot for all dependent
        fig_box_all, ax_box_all = plt.subplots(figsize=(8, 6))
        box_data_all = [bg_all_shift.values, all_dep_shift.values, all_indep_shift.values]
        box_labels_all = [f'Background\n(n={len(bg_all_shift)})', f'UFM1-Dep\n(n={len(all_dep_shift)})', f'UFM1-Indep\n(n={len(all_indep_shift)})']
        box_colors_all = ['grey', 'red', 'blue']
        
        bp_all = ax_box_all.boxplot(box_data_all, tick_labels=box_labels_all, patch_artist=True)
        for patch, color in zip(bp_all['boxes'], box_colors_all):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        
        # Stats
        stat_results_all = []
        if len(bg_all_shift) > 0 and len(all_dep_shift) > 0:
            _, p_val = stats.mannwhitneyu(bg_all_shift, all_dep_shift, alternative='two-sided')
            stat_results_all.append(f'BG vs Dep: p={p_val:.2e}')
        if len(bg_all_shift) > 0 and len(all_indep_shift) > 0:
            _, p_val = stats.mannwhitneyu(bg_all_shift, all_indep_shift, alternative='two-sided')
            stat_results_all.append(f'BG vs Indep: p={p_val:.2e}')
        if len(all_dep_shift) > 0 and len(all_indep_shift) > 0:
            _, p_val = stats.mannwhitneyu(all_dep_shift, all_indep_shift, alternative='two-sided')
            stat_results_all.append(f'Dep vs Indep: p={p_val:.2e}')
        
        ax_box_all.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        ax_box_all.set_ylabel('Shift Metric (log2[N/C ANS / N/C CTRL])', fontsize=12)
        ax_box_all.set_title(f'All UFM1-Dependent: N/C Shift Comparison - {context}', fontsize=14)
        ax_box_all.text(0.98, 0.02, '\n'.join(stat_results_all), transform=ax_box_all.transAxes, fontsize=9,
                       verticalalignment='bottom', horizontalalignment='right',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        box_all_file = os.path.join(args.outdir, f"NC_shift_Boxplot_AllDependent_{context}.pdf")
        fig_box_all.savefig(box_all_file, dpi=150, bbox_inches='tight')
        plt.close(fig_box_all)
        print(f"[DONE] Saved Boxplot (All Dependent) to {box_all_file}")
        
        # --- PLOT B: Separate SRSF3 and PCBP2 ---
        for motif_name, motif_pattern, color in [('SRSF3', 'SRSF3', 'darkred'), ('PCBP2', 'PCBP2', 'darkblue')]:
            print(f"[INFO] Creating {motif_name} motif plot for {context}...")
            
            # Get genes with this motif (from motif_targets and excluded_targets)
            motif_dep_genes = set(motif_targets[motif_targets['motif_id'].str.contains(motif_pattern, case=False, na=False)]['gene_id'].unique())
            motif_excl_genes = set(excluded_targets[excluded_targets['motif_id'].str.contains(motif_pattern, case=False, na=False)]['gene_id'].unique()) if len(excluded_targets) > 0 else set()
            
            motif_dep_mask = nc_df_valid['gene_id'].isin(motif_dep_genes)
            motif_excl_mask = nc_df_valid['gene_id'].isin(motif_excl_genes)
            motif_bg_mask = ~nc_df_valid['gene_id'].isin(motif_dep_genes | motif_excl_genes | all_dep_gene_ids | all_indep_gene_ids)
            
            motif_bg_shift = nc_df_valid.loc[motif_bg_mask, 'Shift_Metric'].dropna()
            motif_dep_shift = nc_df_valid.loc[motif_dep_mask, 'Shift_Metric'].dropna()
            motif_excl_shift = nc_df_valid.loc[motif_excl_mask, 'Shift_Metric'].dropna()
            
            if len(motif_dep_shift) == 0 and len(motif_excl_shift) == 0:
                print(f"[WARN] No genes found for {motif_name} motif. Skipping.")
                continue
            
            # CDF for this motif
            fig_cdf_motif, ax_cdf_motif = plt.subplots(figsize=(8, 6))
            
            if len(motif_bg_shift) > 0:
                sorted_bg = np.sort(motif_bg_shift)
                cdf_bg = np.arange(1, len(sorted_bg) + 1) / len(sorted_bg)
                ax_cdf_motif.plot(sorted_bg, cdf_bg, color='grey', linewidth=2, label=f'Background (n={len(motif_bg_shift)})')
            
            if len(motif_dep_shift) > 0:
                sorted_dep = np.sort(motif_dep_shift)
                cdf_dep = np.arange(1, len(sorted_dep) + 1) / len(sorted_dep)
                ax_cdf_motif.plot(sorted_dep, cdf_dep, color=color, linewidth=2, label=f'{motif_name} Dep-Specific (n={len(motif_dep_shift)})')
            
            if len(motif_excl_shift) > 0:
                sorted_excl = np.sort(motif_excl_shift)
                cdf_excl = np.arange(1, len(sorted_excl) + 1) / len(sorted_excl)
                ax_cdf_motif.plot(sorted_excl, cdf_excl, color='orange', linewidth=2, linestyle='--', label=f'{motif_name} Overlap (n={len(motif_excl_shift)})')
            
            ax_cdf_motif.axvline(x=0, color='black', linestyle='--', alpha=0.5)
            ax_cdf_motif.set_xlabel('Shift Metric (log2[N/C ANS / N/C CTRL])', fontsize=12)
            ax_cdf_motif.set_ylabel('Cumulative Probability', fontsize=12)
            ax_cdf_motif.set_title(f'{motif_name} Motif: N/C Shift Distribution - {context}', fontsize=14)
            ax_cdf_motif.legend(loc='lower right', fontsize=10)
            ax_cdf_motif.grid(True, alpha=0.3)
            
            cdf_motif_file = os.path.join(args.outdir, f"NC_shift_CDF_{motif_name}_{context}.pdf")
            fig_cdf_motif.savefig(cdf_motif_file, dpi=150, bbox_inches='tight')
            plt.close(fig_cdf_motif)
            print(f"[DONE] Saved CDF plot ({motif_name}) to {cdf_motif_file}")
            
            # Boxplot for this motif
            fig_box_motif, ax_box_motif = plt.subplots(figsize=(8, 6))
            
            box_data_motif = []
            box_labels_motif = []
            box_colors_motif = []
            
            if len(motif_bg_shift) > 0:
                box_data_motif.append(motif_bg_shift.values)
                box_labels_motif.append(f'Background\n(n={len(motif_bg_shift)})')
                box_colors_motif.append('grey')
            
            if len(motif_dep_shift) > 0:
                box_data_motif.append(motif_dep_shift.values)
                box_labels_motif.append(f'{motif_name}\nDep-Specific\n(n={len(motif_dep_shift)})')
                box_colors_motif.append(color)
            
            if len(motif_excl_shift) > 0:
                box_data_motif.append(motif_excl_shift.values)
                box_labels_motif.append(f'{motif_name}\nOverlap\n(n={len(motif_excl_shift)})')
                box_colors_motif.append('orange')
            
            bp_motif = ax_box_motif.boxplot(box_data_motif, tick_labels=box_labels_motif, patch_artist=True)
            for patch, col in zip(bp_motif['boxes'], box_colors_motif):
                patch.set_facecolor(col)
                patch.set_alpha(0.6)
            
            # Stats
            stat_results_motif = []
            if len(motif_bg_shift) > 0 and len(motif_dep_shift) > 0:
                _, p_val = stats.mannwhitneyu(motif_bg_shift, motif_dep_shift, alternative='two-sided')
                stat_results_motif.append(f'BG vs {motif_name}: p={p_val:.2e}')
            
            ax_box_motif.axhline(y=0, color='black', linestyle='--', alpha=0.5)
            ax_box_motif.set_ylabel('Shift Metric (log2[N/C ANS / N/C CTRL])', fontsize=12)
            ax_box_motif.set_title(f'{motif_name} Motif: N/C Shift Comparison - {context}', fontsize=14)
            ax_box_motif.text(0.98, 0.02, '\n'.join(stat_results_motif), transform=ax_box_motif.transAxes, fontsize=9,
                             verticalalignment='bottom', horizontalalignment='right',
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            box_motif_file = os.path.join(args.outdir, f"NC_shift_Boxplot_{motif_name}_{context}.pdf")
            fig_box_motif.savefig(box_motif_file, dpi=150, bbox_inches='tight')
            plt.close(fig_box_motif)
            print(f"[DONE] Saved Boxplot ({motif_name}) to {box_motif_file}")
    
    print("[INFO] === Step 15 Completed ===")


if __name__ == "__main__":
    main()
