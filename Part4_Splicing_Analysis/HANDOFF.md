# Project Handoff: Splicing Analysis Pipeline

## Project Status (As of Jan 21, 2026)
This repository contains a comprehensive splicing analysis pipeline for Human, Mouse, and Arabidopsis. It performs motif enrichment, functional impact, and subcellular distribution analyses.

## Key Features Implemented
-   **Multi-Species Support**: Human, Mouse, Arabidopsis.
-   **15-Step Pipeline**: From data prep through subcellular distribution analysis.
-   **Persistent Caching**: TxDb and Biomart data cached in `data/cache`.
-   **Configurable Thresholds**: `--fdr`, `--dpsi`, `--min-reads`, `--event_types`.

## Recent Changes (Jan 21, 2026)

### Step 15: Motif Target Identification & Subcellular Distribution (Human Only)
-   **New Script**: `src/step15_motif_subcellular.py`
    -   Scans UFM1-dependent introns for SRSF3/PCBP2 motifs.
    -   Calculates Nucleus/Cytosol (N/C) ratios from gene expression data.
    -   Generates N/C shift scatter plots (DMSO vs DOX contexts).
    -   Produces CDF and Boxplot comparisons (Background vs Dependent-Specific vs Overlap).
    -   Outputs: `NC_shift_combined.png`, `NC_shift_CDF_*.pdf`, `NC_shift_Boxplot_*.pdf`.

### Step 15B: Expression Boxplots by Gene Set
-   **New Script**: `src/step15b_expression_boxplots.py`
    -   Compares WT (DMSO) vs UFM1 (DOX) expression levels.
    -   Gene Sets: Background, All UFM1-Dependent, All UFM1-Independent, Motif-Dependent.
    -   Generates faceted 2x2 boxplot: `Expression_Boxplots_Combined.png`.
    -   Statistical tests (Wilcoxon) annotated on plots.

### Functional Enrichment Analysis
-   **New Script**: `src/run_functional_enrichment.R`
    -   Uses `gprofiler2` for GO/KEGG enrichment.
    -   Outputs: `*_enrichment_results.csv`, `*_top_terms.csv`, `*_enrichment_plot.pdf`.

### Pipeline Enhancements
-   **Log Directory**: Logs now saved in `<outdir>/logs/` (previously just `logs/`).
-   **dPSI Threshold**: Restored to original event-specific defaults (SE=0.2, others=0.1).

## Example Runs

### Human RI with Strict dPSI (0.25)
```bash
python3 run_pipeline.py --species human --event_types RI --dpsi 0.25 --outdir hs_RI_0.25
```

### All Events Except SE (FDR=0.01, dPSI=0.1)
```bash
python3 run_pipeline.py --species human --event_types A3SS A5SS MXE RI --fdr 0.01 --steps 1 --outdir hs_allExceptSE_fdr0.01_dpsi0.1
```

### Machine Learning Pipeline: Deep Feature Hunter
-   **New Directory**: `machine_learning/`
-   **Main Script**: `machine_learning/run_deep_feature_hunter.R`
    -   Classifies UFM1-dependent vs independent introns using Random Forest.
    -   **Feature Evolution**:
        -   **v1**: Sequence composition (GC, PPT, CpG) and basic motifs (SRSF3, PCBP2).
        -   **v2**: 45 features including all SRSF motifs, PCBP/HNRNP/RBM families, and Branch Point scoring.
        -   **v3**: 53 features adding **Contextual/Kinetic** features (Upstream Exon GC, Codon Optimality, Decoy Competition, RNA Accessibility).
    -   **Multi-Direction Analysis (FDR 0.05)**:
        -   **dPSI Positive (CTRL > ANS)**: 443 introns. Anisomycin *decreases* retention (Kinetic Window Hypothesis). AUC 0.553.
        -   **dPSI Negative (ANS > CTRL)**: 144 introns. Anisomycin *increases* retention (Blocking Hypothesis). AUC 0.607.
    -   **Top Features**: `GC_intron`, `intron_length`, `PCBP2_density`, `SRSF3_density`, `upstream_exon_GC`.

## Scientific Insights (Jan 21, 2026)
1.  **Kinetic Window Hypothesis**: The majority of UFM1-dependent introns are "fast" introns that are normally retained because the ribosome moves too quickly. Anisomycin stalling "rescues" their splicing.
2.  **RBP Coupling**: SRSF3 and PCBP2 density are consistently among the most discriminative features, suggesting these RBPs may be sensitive to ribosome-mediated splicing regulation.
3.  **Predictability**: Introns inhibited by stalling (dPSI Negative) are more sequence-predictable (AUC 0.6) than those rescued by stalling (AUC 0.55).

## Key Outputs

| Step | Output Directory | Key Files |
|------|------------------|-----------|
| 1 | `step01_data_prep/` | `UFM1_dependent.tsv`, `UFM1_independent.tsv` |
| 10 | `step10_motif_analysis/` | AME results, motif comparison plots |
| 14 | `step14_genomic_associations/` | miRNA, NMD, EJC/PTC, Feature Lengths |
| 15 | `step15_subcellular_distribution/` | N/C shift plots, Expression boxplots |
| ML | `machine_learning/outputs_v3/` | `Variable_Importance.png`, `ROC_Curve.pdf` |

## Maintenance Protocol
1.  **Start**: Read this `HANDOFF.md`.
2.  **Work**: Update `task.md` and code.
3.  **End**: Update this file with new changes.
