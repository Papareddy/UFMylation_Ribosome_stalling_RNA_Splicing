# Splicing Functional Impact Pipeline

This pipeline analyzes alternative splicing events from rMATS output to determine their functional impact, particularly in the context of UFMylation stress. It identifies domains enriched in differentially spliced genes, analyzes protein attributes (SignalP, TMHMM, NCOILS), and visualizes the results.

## Pipeline Steps

The pipeline is orchestrated by `run_pipeline.py` and executes the following analysis steps sequentially:

1.  **Data Preparation (`step01_prep`)**:
    -   Classifies rMATS splicing events into 'Lost' (UFMylation-dependent) and 'Preserved' (UFMylation-independent) categories.
    -   **Filtering Logic**: Applies event-specific |dPSI| thresholds: **≥ 0.2 for SE**, and **≥ 0.1 for other types** (RI, A3SS, A5SS, MXE). FDR < 0.05.
2.  **SpliceImpactR Integration (`step02_splice_impact_features`)**:
    -   Fetches functional domains (Pfam/InterPro) and protein attributes (SignalP, TMHMM, NCOILS) from Ensembl Biomart.
    -   Uses persistent caching (`data/cache`) to avoid redundant downloads.
3.  **Domain Enrichment Plotting (`step03_domain_enrichment`)**:
    -   Calculates Fisher's exact test statistics for domain enrichment in differentially spliced genes.
    -   Generates Volcano plots visualizing enrichment significance vs. odds ratio.
4.  **Protein Attribute Enrichment Plotting (`step04_protein_attributes`)**:
    -   Visualizes enrichment of biophysical features (Signal Peptides, Transmembrane domains, Coiled-coils) using grouped bar charts.
5.  **Functional Impact Classification (`step05_functional_impact`)**:
    -   Predicts the consequence of splicing events, such as Frameshifts, NMD induction, or Start/Stop codon disruption.
6.  **Protein Sequence Impact (`step06_protein_sequence_impact`)**:
    -   Synthesizes functional annotations into high-level categories (e.g., `Frameshift_likely_NMD`, `Inframe_CDS_change`).
    -   Generates "Impact Fraction" mirror plots comparing Lost vs. Preserved events.
    -   Generates **Alignment Score Histogram** (Frequency distribution) verifying isoform similarity.
7.  **Frame Shift Density Analysis (`step07_frameshift_density`)**:
    -   Maps splicing events to relative positions within the coding sequence (0.0 to 1.0).
    -   Analyzes and plots the density of frameshifting events across the gene body to identify positional biases.
8.  **AA Feature Extraction (`step08_aa_features`)**:
    -   Extracts nucleotide sequences from the reference genome for the identified splicing events.
    -   Translates sequences to amino acids, identifying the best Open Reading Frame (ORF).
9.  **Protein Attribute Calculation (`step09_biophysical_properties`)**:
    -   Computes detailed physicochemical properties (MW, pI, hydrophobicity) for the extracted amino acid sequences.
10. **Deep Dive Motif Analysis (`step10_motif_analysis`)**:
    -   **Consolidated Feature Analysis**: Visualizes GC content, Sequence Length, and MaxEntScan splice site scores (5'/3') for proper splicing regulation (Lost vs Preserved vs Constitutive).
    -   **Motif Enrichment**: Performs targeted motif enrichment analysis using **AME** (with relaxed E-value reporting) and DREME. Supports **CisBP** and **RBPDB** databases.
    -   **Motif Comparison**: Generates scatter plots comparing motif enrichment significance (Lost vs Preserved).
11. **Mechanism Investigation (`step11_mechanism_investigation`)**:
    -   **Stalling Analysis**: Tests for "Hard-to-translate" peptide features (Poly-Basic stretches, Proline density, Rare codons).
    -   **Adjacency Analysis**: Checks for spatial co-occurrence of SE and RI events to test regulatory coupling.
    -   **Sequence Properties**: Verifies GC content and exon length stability.


## Directory Structure

```text
.
├── download_genomes.sh       # Script to download reference genomes
├── run_pipeline.py           # Main execution script
├── splicing-functional.yml   # Conda/Mamba environment file
├── data/                     # Input data (GTF, rMATS output)
│   ├── cache/                # Persistent cache (auto-generated)
│   ├── human/
│   └── mouse/
├── results/                  # Pipeline outputs (Plots, Tables)
└── src/                      # Source code for analysis steps
```

## Requirements

### Software

All required software and libraries are specified in the `splicing-functional.yml` environment file.

It is highly recommended to use [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) to create the environment, as it is significantly faster than Conda.

Key dependencies include:
- Mamba (or Conda)
- Python (with pandas, numpy, biopython, matplotlib)
- R (with tidyverse)
- `bedtools`
- `samtools`

### Environment Setup
```bash
mamba env create -f splicing-functional.yml
mamba activate splicing-functional

```
If you do not have Mamba, you can use Conda, although it will be much slower: `conda env create -f splicing-functional.yml`.

### Data

The following data files are required to run the pipeline:
- **rMATS output directories**: The `data` directory contains the required rMATS output files, organized by species.
- **Gene Annotation**: The required GTF files (compressed) are included in the `data` directory.
- **Reference Genome**: A FASTA file (`.fa`) for the organism of interest. These are not stored in Git and must be downloaded separately.
- **(Optional) Protein Translations**: A FASTA file containing protein sequences. This is used for optional alignment-based analysis in the protein impact step.

#### Genome Files (Required Download)

The reference genome FASTA files (`.fa`) are too large to be stored in this repository and must be downloaded separately. After creating the `splicing-functional` environment, run the following script from within the `Figure4_splicing_analysis` directory:

```bash
bash download_genomes.sh
```
This will download and index the required human and mouse genomes from Ensembl.

### Annotation Versions

To ensure reproducibility, the specific versions of the gene annotation files (GTF) used in our analysis are listed below. If you plan to re-analyze the data with a different annotation version, be aware that the results may vary.

-   **Human:** GENCODE Release 45 (GRCh38) - `gencode.v45.pc_translations.fa` and `gencode.v45.annotation.gtf`
-   **Mouse:** Ensembl Release 112 (GRCm39) - `Mus_musculus.GRCm39.112.gtf`
-   **Arabidopsis:** Ensembl Plants Release 56 (TAIR10) - `Arabidopsis_thaliana.TAIR10.56.gtf` (from [Ensembl Plants FTP](ftp://ftp.ensemblgenomes.org/pub/release-56/plants/gtf/arabidopsis_thaliana/))

## How to Run

1.  Create and activate the Mamba environment as described in the [Software](#software) section.

2.  Download and index the required reference genomes (after activatig mamba env).
    ```bash
    bash download_genomes.sh
    ```

3.  Execute the main pipeline script, specifying the species to analyze:
    ```bash
    # For human
    python run_pipeline.py --species human

    # For mouse
    python run_pipeline.py --species mouse
    ```
    You can adjust other pipeline arguments as needed.

### Parameters for `run_pipeline.py`

The main pipeline script `run_pipeline.py` accepts several arguments to customize its behavior:

*   `--species` (required):
    *   **Description**: Specifies the species for analysis.
    *   **Choices**: `human`, `mouse`, `arabidopsis`.
    *   **Example**: `--species human`

*   `--fraction` (optional):
    *   **Description**: Specifies the cellular fraction to analyze.
    *   **Choices**: `nucleus`, `cytosol`.
    *   **Default**: `nucleus` (for human and arabidopsis if not specified).
    *   **Example**: `--fraction cytosol`

*   `--outdir` (optional):
    *   **Description**: The main output directory where all results will be stored.
    *   **Default**: `results`
    *   **Example**: `--outdir my_analysis_results`

*   `--no-fasta` (optional):
    *   **Description**: If set, skips protein alignment even if a protein FASTA file is found in the data directory. This will skip "Step 2: Protein Primary Sequence Impact" and related plotting.
    *   **Type**: Flag (no value needed).
    *   **Example**: `--no-fasta`

*   `--fdr` (optional):
    *   **Description**: False Discovery Rate (FDR) threshold for rMATS filtering.
    *   **Default**: `0.05`
    *   **Example**: `--fdr 0.01`

*   `--fdr_domain` (optional):
    *   **Description**: FDR threshold for domain enrichment significance.
    *   **Default**: `0.01`
    *   **Example**: `--fdr_domain 0.05`

*   `--dpsi` (optional):
    -   **Description**: *Legacy argument*. The pipeline now applies event-specific thresholds defined internally in Step 1 (SE ≥ 0.2, Others ≥ 0.1).
    -   **Default**: `0.15` (Effective: SE=0.2, others=0.1)
    -   **Example**: `--dpsi` (value ignored for step1 filtering)

*   `--min-reads` (optional):
    *   **Description**: Minimum reads threshold for splicing events.
    *   **Default**: `20`
    *   **Example**: `--min-reads 10`
    
*   `--background` (optional):
    *   **Description**: Background for enrichment analysis. 'genome' (all genes) or 'rmats' (only genes tested in rMATS).
    *   **Default**: `genome`
    *   **Example**: `--background rmats`

*   `--normalize` (optional, for `add_frame_shify_density.py`):
    *   **Description**: Normalization method for frame shift density calculation.
    *   **Default**: `log2ratio`
    *   **Example**: `--normalize count`

*   `--nperm` (optional, for `add_frame_shify_density.py`):
    *   **Description**: Number of permutations for statistical testing in frame shift density.
    *   **Default**: `1000`
    *   **Example**: `--nperm 5000`

*   `--nbins` (optional, for `add_frame_shify_density.py`):
    *   **Description**: Number of bins for frame shift density plots.
    *   **Default**: `5`
    *   **Example**: `--nbins 10`

*   `--pool` (optional, for `add_frame_shify_density.py`):
    *   **Description**: If set, pools data for frame shift density calculation.
    *   **Type**: Flag (no value needed).
    *   **Example**: `--pool`

*   `--event_types` (optional, for `add_frame_shify_density.py`):
    *   **Description**: List of event types to include in frame shift density analysis.
    *   **Default**: `SE` (Skipped Exon). Can be multiple types.
    *   **Example**: `--event_types SE A3SS A5SS`

*   `--direction` (optional, for `AAfeatures.sh`):
    *   **Description**: If set, enables direction-based splitting in the `AAfeatures.sh` script.
    *   **Type**: Flag (no value needed).
    *   **Example**: `--direction`
    
*   `--cache-dir` (optional):
    *   **Description**: Directory to store persistent caching (TxDb, Biomart).
    *   **Default**: `data/cache`
    *   **Example**: `--cache-dir /tmp/mycache`

#### Example Usage with Parameters:
```bash
python run_pipeline.py --species human --fraction cytosol --fdr 0.01 --dpsi 0.1 --event_types SE MXE --direction
```

## Outputs

Results are saved in `results/{species}/{fraction}/`. The output is organized into summary plots, data tables, and specific step directories.

### 1. Summary Visualizations (Step-Specific Directories)
*   **`step03_domain_enrichment/Volcano_Domain_Enrichment_Comparison.png`**: A volcano plot showing Pfam/InterPro domains enriched in UFMylation-dependent ('Lost') vs. Independent ('Preserved') splicing events.
*   **`step04_protein_attributes/Protein_Attributes_Enrichment.png`**: A grouped bar chart displaying the enrichment of Signal Peptides, Transmembrane Domains, and Coiled-coils in the event subsets.
*   **`step06_protein_sequence_impact/impact_fractions.png`**: A "Mirror" bar chart comparing the distribution of functional impacts (e.g., Frameshift, NMD, UTR only) between Lost and Preserved events.

### 2. Data Tables (Subdirectories)
*   **`step03_domain_enrichment/domain_enrichment.tsv`**: Full statistical results (Fisher's exact test) for domain enrichment, including P-values, Odds Ratios, and FDR.
*   **`step04_protein_attributes/biophysical_enrichment.tsv`**: Statistical results for Protein Attribute enrichment.
*   **`step01_prep/lost.tsv`**: The subset of rMATS events classified as 'Lost' (dependent on UFMylation).
*   **`step01_prep/preserved.tsv`**: The subset of rMATS events classified as 'Preserved' (independent).
*   **`step09_biophysical_properties/protein_attributes_properties.tsv`**: Detailed physicochemical properties (Molecular Weight, pI, Hydrophobicity, Instability Index) calculated for every protein isoform derived from the splicing events.

### 3. Frameshift Analysis (`step07_frameshift_density/`)
*   **`fig4__event_SE__class_*.png`**: Density plots visualizing *where* frameshifts occur along the coding sequence (normalized 0.0 to 1.0).
    *   Example: `fig4__event_SE__class_NMD_likely.png` shows the positional distribution of frameshifts predicted to trigger NMD.
*   **`*.bin_fisher.tsv`**: Statistical tests for positional enrichment of frameshifts (e.g., N-terminal vs C-terminal bias).

### 4. Intermediate Files
*   **`step08_aa_features/`**: Contains extracted FASTA files:
    *   `*.nt.fa`: Nucleotide sequences of the spliced regions.
    *   `*.aa.fa`: Translated Amino Acid sequences for the best Open Reading Frame.
*   **`logs/`**: Execution logs for each step of the pipeline.
*   **`step10_motif_analysis/`**:
    *   **`master_features.tsv`**: Consolidated table of GC content, length, and MaxEnt scores.
    *   **`combined_features.pdf`**: Multi-panel plot visualizing feature distributions.
    *   **`motif_comparison.png`**: Scatter plot of AME enrichment (if run with RBPDB/CisBP).


