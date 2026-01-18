# Project Handoff: Splicing Analysis Pipeline

## Project Status (As of Jan 16, 2026)
This repository contains a splicing analysis pipeline that allows for:
1.  **Domain Enrichment**: Identifying Pfam domains enriched in alternatively spliced genes.
2.  **Protein Attribute Enrichment** (formerly Biophysical): Analyzing Signal Peptides, Transmembrane Domains, and Coiled-Coils.
3.  **Visualization**: Generating Volcano plots and grouped bar charts.

## Key Features Implemented
-   **Multi-Species Support**: Validated for **Human** and **Mouse**.
    -   Automatically selects Ensembl datasets (`hsapiens_gene_ensembl`, `mmusculus_gene_ensembl`).
-   **Persistent Caching**:
    -   All Biomart downloads and TxDb objects are cached in `data/cache` (configurable via `--cache-dir`).
    -   Reduces runtime and dependency on Ensembl availability.
-   **Streamlined Pipeline**:
    -   Refactored into 9 sequential, clearly named steps.
    -   Consistent terminology ("Protein Attributes" instead of "Biophysical").
-   **Directionality (`--direction`)**:
    -   Splits events into **Inclusion (Inc)** (dPSI > 0) and **Exclusion (Exc)** (dPSI < 0).
    -   Generates faceted plots for these subsets.
-   **Background Selection (`--background`)**:
    -   `genome` (Default): Uses all genes in the annotation.
    -   `rmats`: Uses only genes tested in rMATS.
-   **Robustness**:
    -   Handles missing FASTA files gracefully (skips Step 0.5/Step 2).
    -   Split FDR control: `--fdr` (Event filtering) vs `--fdr_domain` (Enrichment significance).
    -   Robust BioMart fetching (smaller chunks, retries).

## How to Run

### 1. Standard Run (Human)
```bash
mamba run -n splicing-functional python3 run_pipeline.py --species human --no-fasta --background genome
```

### 2. Mouse Analysis
```bash
mamba run -n splicing-functional python3 run_pipeline.py --species mouse --no-fasta --background genome
```

### 3. Usage with Caching
By default, the pipeline uses `data/cache`. To specify a custom location:
```bash
python3 run_pipeline.py --species human --cache-dir /path/to/custom/cache
```

## Repository Structure
-   `run_pipeline.py`: Main entry point (Steps 1-9).
-   `src/get_splice_impact_features.R`: Core R logic for enrichment (Impact + Protein Attributes).
-   `src/analyze_domain_enrichment.py`: Python plotting logic for Volcano and Bar charts.
-   `src/biophysical_properties.py`: Computes biophysical traits (legacy/internal).
-   `data/cache/`: Directory containing persistent `.sqlite` (TxDb) and `.rds` (Biomart) files.

## Recent Changes
-   **Refactoring**: `run_pipeline.py` steps renamed and reordered for clarity.
-   **Renaming**: "Biophysical" -> "Protein Attributes" in plots and logs.
-   **Caching**: Implemented persistent caching for TxDb and Biomart queries to prevent redundant downloads.
-   **Fixes**: Corrected `download_genomes.sh` argument parsing and R script error handling.

## Next Steps for New Agent
-   Check `context/task.md` for the full checklist of completed items.
-   Refer to `context/walkthrough.md` for visual examples of expected outputs.

## Future Maintenance Protocol
To ensure this context survives for the *next* session, please follow this workflow:
1.  **Start**: Read this file (`HANDOFF.md`) to get context.
2.  **Work**: update `task.md` and code as usual.
3.  **End**: Before the user ends the session, update this `HANDOFF.md` file with your new changes (new features, fixed bugs).
4.  **Sync**: Copy the updated artifacts (`HANDOFF.md`, `task.md`, etc.) back to `context/` and git push.
