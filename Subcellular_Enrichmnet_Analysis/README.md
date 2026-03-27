# Figure 2: Subcellular Compartment Enrichment Analysis

This project analyzes subcellular compartment proteomics and transcriptomics data to identify enrichment patterns of proteins and transcripts across different cellular fractions. The pipeline is designed to process quantitative data, perform statistical analysis for enrichment, and visualize the results.

## Requirements

### Software
- R (version 4.0 or higher)
- R packages: `ggplot2`, `dplyr`, `tidyr`, `readr`
- bash shell

### Data
- Input data should be a CSV or TSV file containing columns for identifiers (e.g., GeneID, ProteinID), log2 fold change, and p-values.
- A sample data file is located in `data/Microsome_ANS_vsDMSO_Col0_ufm1_Limma_with_log2counts.csv`.

## Directory Structure

```
.
├── data/
│   └── Microsome_ANS_vsDMSO_Col0_ufm1_Limma_with_log2counts.csv
├── results/
│   ├── enrichment_results.csv
│   └── enrichment_plot.png
├── src/
│   ├── 01_preprocess.R
│   ├── 02_enrichment.R
│   └── 03_visualize.R
└── run_pipeline.sh
```

- `data/`: Contains raw input data.
- `results/`: Stores output files from the pipeline, such as tables and plots.
- `src/`: Contains the R scripts for each step of the analysis.

## Pipeline Execution

The pipeline is executed via the `run_pipeline.sh` script.

```bash
bash run_pipeline.sh --pval 0.05 --fc_threshold 1.5 --top_n 20
```

### Control Flags

- `--pval`: Sets the p-value threshold for determining significance. Default is `0.05`.
- `--fc_threshold`: Sets the log2 fold change threshold for enrichment. Default is `1.0`.
- `--top_n`: Specifies the number of top enriched items to highlight in visualizations. Default is `20`.

## Pipeline Steps

1.  **Data Pre-processing (`src/01_preprocess.R`):**
    - Loads the input data from the `data/` directory.
    - Filters the data based on user-defined p-value and fold-change thresholds.
    - Prepares the data for enrichment analysis.

2.  **Compartment Enrichment (`src/02_enrichment.R`):**
    - Performs a Fisher's Exact Test or Hypergeometric Test to determine the statistical significance of enrichment for each subcellular compartment.
    - Generates a table of enriched compartments with corresponding statistics.

3.  **Visualization (`src/03_visualize.R`):**
    - Creates visualizations (e.g., volcano plots, bar charts) of the enrichment results.
    - Highlights the top `N` enriched compartments or items.
    - Saves the plots to the `results/` directory.
