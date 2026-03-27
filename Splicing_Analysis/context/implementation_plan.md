# Implementation Plan - Arabidopsis Microsome Enrichment Fix

## User Review Required
> [!IMPORTANT]
> I will modify `src/Arabidopsis_microsome_enrichment.R` to **always** run analysis for both "RI" and "Others" (SE/MXE/etc) groups, regardless of command line flags. This ensures both plots are generated consistently.
> I will also update `run_pipeline.py` to call this script once, simplifying the pipeline.

## Proposed Changes

### [Arabidopsis Script]
#### [MODIFY] [src/Arabidopsis_microsome_enrichment.R](file:///Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/src/Arabidopsis_microsome_enrichment.R)
-   Remove `event_group` argument dependency for controlling flow.
-   Implement a loop/list strategy to process:
    1.  **RI**: `EventType == "RI"`
    2.  **Others**: `EventType %in% c("SE", "MXE", "A3SS", "A5SS")`
-   For each group:
    -   Filter `genes_preserved` and `genes_lost`.
    -   Generate specific PDF output (`microsome_enrichment_RI.pdf`, `microsome_enrichment_Others.pdf`).

### [Pipeline Runner]
#### [MODIFY] [run_pipeline.py](file:///Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/run_pipeline.py)
-   Steps 19-20 (logic for Step 12):
    -   Replace the double call (RI vs Others) with a single call to `src/Arabidopsis_microsome_enrichment.R`.
    -   Remove `--event_group` flag passing (or keep it as dummy if needed, but script will ignore it for selection).

## Verification Plan

### Automated Tests
-   Run Step 12 via `run_pipeline.py`:
    ```bash
    mamba run -n splicing-functional python3 run_pipeline.py --species arabidopsis --fraction nucleus --steps 12
    ```
-   Check output directory `results/arabidopsis/nucleus/step12_microsome_enrichment/` for existence of:
    -   `microsome_enrichment_RI.pdf`
    -   `microsome_enrichment_Others.pdf`

### Manual Verification
-   Inspect the logs to ensure "Event Type(s): RI" and "Event Type(s): Others" sections appear and report correct gene counts (expected ~27/100 for RI, 1/5 for Others based on current data).
