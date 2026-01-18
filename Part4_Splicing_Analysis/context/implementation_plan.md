# Biophysical Feature Enrichment (SignalP, TMHMM, NCOILS)

## Goal
Calculate and visualize the enrichment of specific biophysical features (Signal Peptides, Transmembrane Helices, Coiled-Coils) in splicing events.

## Proposed Changes

### R Script
#### [MODIFY] [src/get_splice_impact_features.R](file:///Users/ranjith.papareddy/.gemini/antigravity/scratch/Figure4_splicing_analysis/src/get_splice_impact_features.R)
- **Data Fetching**: Update `biomaRt` query to include `signalp`, `tmhmm`, `ncoils`.
- **Processing**:
    -   Create binary flags for each gene: `HasSignalP`, `HasTMHMM`, `HasNCOILS` (TRUE if attribute is non-empty).
- **Statistics**:
    -   Perform Fisher's Exact Test for each feature:
        1.  **Lost vs Genome**
        2.  **Preserved vs Genome**
        3.  **Lost vs Preserved** (Direct comparison)
- **Output**: Save `biophysical_enrichment.tsv`.

### Python Script
#### [MODIFY] [src/analyze_domain_enrichment.py](file:///Users/ranjith.papareddy/.gemini/antigravity/scratch/Figure4_splicing_analysis/src/analyze_domain_enrichment.py)
- **Input**: Read `biophysical_enrichment.tsv`.
- **Visualization**:
    -   Generate a Grouped Bar Chart (`Biophysical_Enrichment.png`).
    -   Y-axis: Odds Ratio (Log2 scale?).
    -   Groups: SignalP, TMHMM, NCOILS.
    -   Bars: Lost vs Genome, Preserved vs Genome.
    -   Annotate with significance stars.

## Directionality Support (New Request)

### Goal
Enable `--direction` flag to split analyses by splicing direction (dPSI > 0 vs dPSI < 0).

### R Script
#### [MODIFY] [src/get_splice_impact_features.R](file:///Users/ranjith.papareddy/.gemini/antigravity/scratch/Figure4_splicing_analysis/src/get_splice_impact_features.R)
-   Accept `--direction` flag.
-   If enabled:
    -   Split `Lost` into `Lost_Inc` (dPSI > 0) and `Lost_Exc` (dPSI < 0).
    -   Split `Preserved` into `Preserved_Inc` (dPSI > 0) and `Preserved_Exc` (dPSI < 0).
    -   Calculate enrichment for all 4 subgroups (plus distinct comparisons if needed, e.g. Lost_Inc vs Genome).
-   Update `domain_enrichment.tsv` and `biophysical_enrichment.tsv` to use these new set names (or add a direction column).

### Python Script
#### [MODIFY] [src/analyze_domain_enrichment.py](file:///Users/ranjith.papareddy/.gemini/antigravity/scratch/Figure4_splicing_analysis/src/analyze_domain_enrichment.py)
-   Handle 4 sets (`Lost_Inc`, `Lost_Exc`, `Preserved_Inc`, `Preserved_Exc`) instead of 2.
-   **Volcano Plot**: Can we fit 4 on one plot? Or maybe 2 separate plots (Inclusion vs Low)?
    -   Proposed: Keep side-by-side, but maybe split into 4 facets or 2 rows if direction is on.
    -   Or just color code? Faceting (2x2 or 1x4) is safer.
-   **Biophysical Plot**: Grouped bar chart with 4 groups instead of 2.

### Pipeline
#### [MODIFY] [run_pipeline.py](file:///Users/ranjith.papareddy/.gemini/antigravity/scratch/Figure4_splicing_analysis/run_pipeline.py)
-   Pass `--direction` to `run_splice_impact`.

## Verification
-   Run with `--direction` and check output plots.
