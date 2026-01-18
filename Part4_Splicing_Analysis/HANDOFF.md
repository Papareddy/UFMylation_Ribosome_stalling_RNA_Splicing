# Project Handoff: Splicing Analysis Pipeline

## Project Status (As of Jan 17, 2026)
This repository contains a splicing analysis pipeline that allows for:
1.  **Domain Enrichment**: Identifying Pfam domains enriched in alternatively spliced genes.
2.  **Protein Attribute Enrichment** (formerly Biophysical): Analyzing Signal Peptides, Transmembrane Domains, and Coiled-Coils.
3.  **Visualization**: Generating Volcano plots and grouped bar charts.

## Key Features Implemented
-   **Multi-Species Support**: Validated for **Human** and **Mouse**.
    -   Automatically selects Ensembl datasets (`hsapiens_gene_ensembl`, `mmusculus_gene_ensembl`).
    -   Handles directory structures without 'nucleus' subfolder (e.g. Mouse data).
-   **Persistent Caching**:
    -   All Biomart downloads and TxDb objects are cached in `data/cache` (configurable via `--cache-dir`).
    -   Reduces runtime and dependency on Ensembl availability.
-   **Streamlined Pipeline**:
    -   Refactored into 9 sequential, clearly named steps.
    -   Consistent terminology ("Protein Attributes" instead of "Biophysical").
-   **Directionality (`--direction`)**:
    -   Splits events into **psI** (Promote Splicing Inclusion, dPSI > 0) and **psE** (Promote Splicing Exclusion, dPSI < 0).
    -   Generates faceted plots for these subsets.
-   **Background Selection (`--background`)**:
    -   `genome` (Default): Uses all genes in the annotation.
    -   `rmats`: Uses only genes tested in rMATS.
-   **Robustness**:
    -   Handles missing FASTA files gracefully (skips Step 0.5/Step 2).
    -   Split FDR control: `--fdr` (Event filtering) vs `--fdr_domain` (Enrichment significance).
    -   Robust BioMart fetching (smaller chunks, retries).
    -   Robust Protein FASTA detection (handles `*pc_translations*` and `*.pep.*fa`).

## How to Run

### 1. Standard Run (Human)
```bash
mamba run -n splicing-functional python3 run_pipeline.py --species human --background genome
```

### 2. Mouse Analysis (Full Features)
```bash
mamba run -n splicing-functional python3 run_pipeline.py --species mouse --background genome --direction
```
Requires `*.pep.all.fa` or `*pc_translations.fa` in `data/mouse/`.

## Repository Structure
-   `run_pipeline.py`: Main entry point (Steps 1-9).
-   `src/get_splice_impact_features.R`: Core R logic for enrichment (Impact + Protein Attributes).
-   `src/splicing_functional_impat.py`: Upstream annotation script. Outputs `inclusion_transcript_ids` and `exclusion_transcript_ids`.
-   `src/protein_primary_sequence_impact.py`: Refactored protein impact analysis. Uses `pysam` for efficient FASTA handling and precise isoform alignment.
-   `src/analyze_domain_enrichment.py`: Python plotting logic for Volcano and Bar charts.
-   `src/biophysical_properties.py`: Computes biophysical traits (legacy/internal).
-   `data/cache/`: Directory containing persistent `.sqlite` (TxDb) and `.rds` (Biomart) files.

## Recent Changes (Jan 17, 2026)
-   **Refactored Protein Impact Analysis**:
    -   **Biological Logic**: `make_fig1F` now requires explicit matching of `inclusion_transcript_ids` and `exclusion_transcript_ids` (provided by `splicing_functional_impat.py`). Alignment is performed only on valid opposing pairs.
    -   **Performance**: Swapped dictionary-based FASTA loading for `pysam.FastaFile` (indexed access). Supports uncompressed and bgzip-compressed FASTA.
    -   **Priority**: Updated `Start_Stop_disruption` to have higher priority than frameshifts.
    -   **Robustness**: Improved FASTA header parsing to handle space-delimited Ensembl headers (`transcript:ENSMUST...`). Fixed `.fai` index file incorrectly being selected as protein FASTA.
    -   **Stability**: Refactored `pysam` usage to prevent segmentation faults (single-open/passing objects). Replaced deprecated `Bio.pairwise2` with `Bio.Align.PairwiseAligner`.
    -   **Visualization**: Updated Alignment Score Histogram to use **Frequency** (probability) instead of counts.
-   **Label Renaming**:
    -   Renamed directional labels from `Inc`/`Exc` to `psI`/`psE` across R scripts, pipeline runner, and plot outputs.
-   **Mouse Pipeline Fixes**:
    -   Modified `run_pipeline.py` to auto-detect data directories (fixing 'nucleus' subfolder requirement).
    -   Updated glob patterns to be case-insensitive (matching `ufm` and `UFM`).
    -   Passed `--species` argument to R script to ensure correct BioMart dataset selection (`mmusculus_gene_ensembl`).
-   **Log & Output Cleanliness**:
    -   Renamed output TSVs in Step 6 to be descriptive (e.g., `Alignment_Scores.tsv` instead of `fig1F...`).
    -   Renamed process log names in Step 6.
    -   Suppressed trivial Warnings in R steps.

## Next Steps for New Agent
-   Check `context/task.md` for the full checklist.
-   Refer to `context/walkthrough.md` for visual examples of Key Results.

## Future Maintenance Protocol
To ensure this context survives for the *next* session, please follow this workflow:
1.  **Start**: Read this file (`HANDOFF.md`) to get context.
2.  **Work**: update `task.md` and code as usual.
3.  **End**: Before the user ends the session, update this `HANDOFF.md` file with your new changes (new features, fixed bugs).
4.  **Sync**: Copy the updated artifacts (`HANDOFF.md`, `task.md`, etc.) back to `context/` and git push.
