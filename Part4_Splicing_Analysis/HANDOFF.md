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

## Key Outputs

| Step | Output Directory | Key Files |
|------|------------------|-----------|
| 1 | `step01_data_prep/` | `UFM1_dependent.tsv`, `UFM1_independent.tsv` |
| 10 | `step10_motif_analysis/` | AME results, motif comparison plots |
| 14 | `step14_genomic_associations/` | miRNA, NMD, EJC/PTC, Feature Lengths |
| 15 | `step15_subcellular_distribution/` | N/C shift plots, Expression boxplots |

## dPSI Thresholds (Default)
| Event Type | Threshold |
|------------|-----------|
| SE | ≥ 0.2 |
| RI, A3SS, A5SS, MXE | ≥ 0.1 |

## Next Steps for New Agent
1.  Review `hs_RI_0.25/` results (stricter dPSI threshold analysis).
2.  Review `functional_analysis/` for GO enrichment insights.
3.  Consider extending Step 15/15B to Mouse data if expression data available.

## Maintenance Protocol
1.  **Start**: Read this `HANDOFF.md`.
2.  **Work**: Update `task.md` and code.
3.  **End**: Update this file with new changes.
