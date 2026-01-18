# Project Handoff: Splicing Analysis Pipeline

## Project Status (As of Jan 18, 2026)
This repository contains a splicing analysis pipeline that allows for:
1.  **Domain Enrichment**: Identifying Pfam domains enriched in alternatively spliced genes.
2.  **Protein Attribute Enrichment** (formerly Biophysical): Analyzing Signal Peptides, Transmembrane Domains, and Coiled-Coils.
3.  **Visualization**: Generating Volcano plots and grouped bar charts.
4.  **Microsome Enrichment**: (Arabidopsis) Analyzing enrichment of UFM1-dependent targets in microsome fractions using DESeq2.

## Key Features Implemented
-   **Multi-Species Support**: Validated for **Human**, **Mouse**, and **Arabidopsis**.
    -   Automatically selects Ensembl datasets (`hsapiens_gene_ensembl`, `mmusculus_gene_ensembl`, `athaliana_eg_gene`).
    -   Handles directory structures without 'nucleus' subfolder (e.g. Mouse data).
-   **Persistent Caching**:
    -   All Biomart downloads and TxDb objects are cached in `data/cache` (configurable via `--cache-dir`).
    -   Arabidopsis DESeq2 results are cached in `ar_results/.../DESeq2_results.RData`.
    -   Reduces runtime and dependency on Ensembl availability.
-   **Streamlined Pipeline**:
    -   Refactored into 12 sequential, clearly named steps.
    -   Consistent terminology ("Protein Attributes" instead of "Biophysical").
-   **Directionality (`--direction`)**:
    -   Splits events into **psI** (Promote Splicing Inclusion, dPSI > 0) and **psE** (Promote Splicing Exclusion, dPSI < 0).
    -   Generates faceted plots for these subsets.
-   **Background Selection (`--background`)**:
    -   `genome` (Default): Uses all genes in the annotation.
    -   `rmats`: Uses only genes tested in rMATS.

## How to Run

### 1. Standard Run (Human)
```bash
mamba run -n splicing-functional python3 run_pipeline.py --species human --background genome
```

### 2. Mouse Analysis (Full Features)
```bash
mamba run -n splicing-functional python3 run_pipeline.py --species mouse --background genome --direction
```

### 3. Arabidopsis Analysis (Microsome Enrichment)
```bash
mamba run -n splicing-functional python3 run_pipeline.py --species arabidopsis --event_types RI SE --fraction nucleus
```
*Note: Step 12 is specific to Arabidopsis and requires raw count files in `data/arabidopsis/GSE82041_RAW/` (or configured base dir).*

## Repository Structure
-   `run_pipeline.py`: Main entry point (Steps 1-12).
-   `src/get_splice_impact_features.R`: Core R logic for enrichment (Impact + Protein Attributes).
-   `src/Arabidopsis_microsome_enrichment.R`: DESeq2 analysis for Arabidopsis microsome data.
-   `src/splicing_functional_impat.py`: Upstream annotation script. Outputs `inclusion_transcript_ids` and `exclusion_transcript_ids`.
-   `src/protein_primary_sequence_impact.py`: Refactored protein impact analysis. Uses `pysam` for efficient FASTA handling.
-   `src/analyze_domain_enrichment.py`: Python plotting logic for Volcano and Bar charts.
-   `data/cache/`: Directory containing persistent `.sqlite` (TxDb) and `.rds` (Biomart) files.

## Recent Changes (Jan 18, 2026)
-   **Arabidopsis Microsome Enrichment (Step 12)**:
    -   **New Script**: `src/Arabidopsis_microsome_enrichment.R` implements DESeq2 analysis.
    -   **Pipeline Integration**: Added Step 12 to `run_pipeline.py`. Runs automatically for `--species arabidopsis`.
    -   **Analysis**: Compares `Microsome vs Cytosol` and `RibosomeMBP vs RibosomeFP`.
    -   **Visualization**: Generates boxplots comparing log2FoldChange of Preserved vs Lost vs Genome genes.
-   **Refactored Protein Impact Analysis** (Jan 17):
    -   **Biological Logic**: Requires explicit matching of `inclusion_transcript_ids` and `exclusion_transcript_ids`.
    -   **Performance**: Swapped dictionary-based FASTA loading for `pysam.FastaFile`.
    -   **Mouse Pipeline Fixes**: Auto-detects data directories and handles case-insensitive file matching.

## Next Steps for New Agent
-   Check `context/task.md` for the full checklist.
-   Refer to `context/walkthrough.md` for visual examples of Key Results.

## Future Maintenance Protocol
To ensure this context survives for the *next* session, please follow this workflow:
1.  **Start**: Read this file (`HANDOFF.md`) to get context.
2.  **Work**: update `task.md` and code as usual.
3.  **End**: Before the user ends the session, update this `HANDOFF.md` file with your new changes (new features, fixed bugs).
4.  **Sync**: Copy the updated artifacts (`HANDOFF.md`, `task.md`, etc.) back to `context/` and git push.
