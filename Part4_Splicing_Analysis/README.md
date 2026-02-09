# UFM1 Splicing Functional Impact Pipeline

This repository contains a comprehensive analytical pipeline for defining the **UFM1-dependent splicing program** and its evolutionary conservation across the eukaryotic kingdom (**Human**, **Mouse**, **Arabidopsis**).

## Overview

The project has evolved from a mammalian analysis into a unified **Eukaryotic Hub**, revealing a fundamental split in the cellular response to ribosome stalling:

1. **UFM1-Dependent ("ER Architects")**: Splicing specifically remodels the ER environment and lipid membranes (e.g., `ESYT1`, `TAFAZZIN`).
2. **UFM1-Independent ("Core OS")**: Housekeeping maintenance (e.g., RNA Splicing factors) is universally conserved but independent of UFM1 stress.
3. **Signal Peptide Eraser**: UFM1-dependent splicing acts as a functional switch with 100% penetrance, systematically stripping Signal Peptides from targeted transcripts.

---

## Directory Structure

```
Part4_Splicing_Analysis/
├── data/                          # Reference genomes and annotations
│   ├── arabidopsis/               # Arabidopsis data (TAIR10)
│   ├── human/                     # Human data (GRCh38/GENCODE v45)
│   ├── mouse/                     # Mouse data (GRCm39/Ensembl 112)
│   ├── cache/                     # Cached TxDb and BioMart results
│   └── motifs/                    # RBP motif databases (CisBP, RBPDB)
├── results/                       # Pipeline outputs by species
│   ├── arabidopsis/nucleus/       # Arabidopsis analysis results
│   ├── human/nucleus/             # Human nucleus fraction results
│   ├── mouse/total/               # Mouse total fraction results
│   └── RI_Consolidated_Aggregation.tsv  # Cross-species RI consolidation
├── src/                           # Analysis scripts (R/Python)
│   ├── prepare_rmats_data.R       # Step 1: rMATS data preparation
│   ├── get_splice_impact_features.R  # Step 2: Domain/protein features
│   ├── protein_primary_sequence_impact.py  # Step 4-5: Sequence impact
│   ├── add_frame_shift_density.py # Step 6: Frame shift analysis
│   ├── analyze_go_enrichment.R    # Step 7: GO enrichment
│   ├── plot_rna_map.py            # Step 10: RNA binding maps
│   ├── aggregate_ri_events.py     # RI data aggregation utility
│   └── run_GO_Final_pipeline.sh   # Cross-species GO pipeline
├── logs/                          # Execution logs for each step
├── run_pipeline.py                # Main pipeline orchestrator
├── README.md                      # This file
└── splicing-functional.yml        # Conda environment specification
```

---

## Reference Annotations

### Arabidopsis thaliana
- **Genome Assembly**: TAIR10
- **Gene Annotation**: Ensembl Plants Release 56
- **GTF File**: `no_plastid_no_rRNA.Arabidopsis_thaliana.TAIR10.56.gtf`
- **FASTA**: `Arabidopsis_thaliana_TAIR10.dna.primary_assembly.fa`
- **Protein Translations**: `Arabidopsis_thaliana.TAIR10.pc_translations.fa`

### Homo sapiens (Human)
- **Genome Assembly**: GRCh38
- **Gene Annotation**: GENCODE v45
- **GTF File**: `pcg_gencode.v45.annotation.gtf.gz` (protein-coding genes only)
- **FASTA**: `Homo_sapiens.GRCh38.dna.primary_assembly.fa`
- **Protein Translations**: `gencode.v45.pc_translations.fa`

### Mus musculus (Mouse)
- **Genome Assembly**: GRCm39
- **Gene Annotation**: Ensembl Release 112
- **GTF File**: `Mus_musculus.GRCm39.112.gtf.gz`
- **FASTA**: `Mus_musculus.GRCm39.dna.primary_assembly.fa`
- **Protein Translations**: `Mus_musculus.GRCm39.pc_translations.fa`

---

## Installation and Setup

### 1. Environment Setup
The pipeline uses a multi-language (Python/R) environment managed via Conda/Mamba.

```bash
# Create the environment
mamba env create -f splicing-functional.yml

# Activate the environment
mamba activate splicing-functional
```

### 2. Data Preparation
Ensure reference genomes and annotations are placed in the appropriate `data/<species>/` directories as shown in the directory structure above.

---

## Usage

### Single-Species Analysis
Run the complete pipeline for an individual species:

```bash
# Human (nucleus fraction)
python run_pipeline.py --species human --fraction nucleus

# Mouse (total fraction)
python run_pipeline.py --species mouse --fraction total

# Arabidopsis (nucleus fraction)
python run_pipeline.py --species arabidopsis --fraction nucleus
```

### Cross-Species Conservation Analysis
To regenerate the integrated results across all three species:

```bash
# Automated GO enrichment and landscape plotting
bash src/run_GO_Final_pipeline.sh
```

### Retained Intron (RI) Data Aggregation
A consolidated cross-species RI dataset is available:

```bash
# Generate consolidated RI events
python src/aggregate_ri_events.py

# Output: results/RI_Consolidated_Aggregation.tsv
# Columns: Species, GeneID, GeneName, PSI_Percentage, Dependency
```

---

## Pipeline Overview

The pipeline consists of multiple steps organized into functional phases:

| Step | Analysis Phase | Description | Key Outputs |
|:-----|:---------------|:------------|:------------|
| **1** | Data Preparation | Filter rMATS events by FDR and dPSI; categorize UFM1-dependent/independent | `UFM1_dependent.tsv`, `UFM1_independent.tsv` |
| **2-3** | Protein Impact | Extract protein domains (InterPro) and attributes | `domain_enrichment.tsv` |
| **4-5** | Sequence Impact | Analyze frame shifts, start/stop codon disruption, NMD prediction | `per_event_compact_forplotting.tsv` |
| **6** | Frame Shift Density | Calculate positional frame shift density | `frame_shift_density.tsv` |
| **7** | GO Enrichment | Perform Gene Ontology enrichment analysis | `GO_enrichment_results.tsv` |
| **8-9** | Sequence Features | GC content and MaxEnt splice site scores | `combined_features.pdf` |
| **10** | RNA Binding Maps | Positional motif enrichment (CisBP/RBPDB) | `RNA_Map_*.png` |
| **11** | Translation Stalling | Ribosome stalling metrics analysis | `stalling_stats.tsv` |
| **12** | Genomic Associations | Annotate genomic context and create BED files | `events_annotated.bed` |
| **13** | Cross-Species Hub | Integrate GO results across species | `Eukaryotic_Conserved_Landscape.png` |
| **14** | Signal Peptide Loss | Quantify signal peptide loss in RI transcripts | `SignalP_Loss_Analysis.tsv` |

---

## Key Scientific Findings

### 1. Conserved ER Sectors
The 3-species landscape identifies **ER Organization** and **Lipid Metabolism** as universal eukaryotic targets of UFM1-dependent splicing. This suggests UFM1-mediated splicing is an ancient mechanism for protecting the secretory factory.

### 2. Signal Peptide Erasing (100% Penetrance)
Analysis reveals that 100% of Retained Intron (RI) isoforms derived from Signal Peptide-positive UFM1-dependent genes in Human and Mouse lose their signal peptide. This confirms splicing is used to systematically reroute proteins away from the ER during stress.

### 3. Splicing Factors Autoregulation
Baseline splicing machinery (e.g., `SF3B1`, `SRSF4`) is a conserved target in the UFM1-independent group, likely representing a universal homeostatic mechanism to brake RNA processing during ribosome stalling.

---

## Output Files

### Consolidated RI Dataset
- **File**: `results/RI_Consolidated_Aggregation.tsv`
- **Description**: Cross-species Retained Intron events with PSI percentages
- **Columns**:
  - `Species`: human, mouse, or arabidopsis
  - `GeneID`: Gene identifier (e.g., ENSG..., ENSMUSG..., AT...)
  - `GeneName`: Gene symbol
  - `PSI_Percentage`: Percent Spliced In (positive = retention in control/WT)
  - `Dependency`: UFM1-Dependent or UFM1-Independent

### Species-Specific Results
Each species has results organized under `results/<species>/<fraction>/`:
- `step01_data_prep/`: Filtered rMATS events
- `step02_domain_enrichment/`: Protein domain analysis
- `step04_protein_impact/`: Sequence impact predictions
- `step07_go_enrichment/`: GO enrichment results
- `step10_rna_maps/`: RNA binding protein positional enrichment
- Additional step-specific outputs

---

## Manual Analysis Tools

Standalone scripts are provided in `src/` for custom analyses:

- **`src/aggregate_ri_events.py`**: Generate consolidated RI dataset across species
- **`src/analyze_signalp_loss.R`**: Background-controlled Odds Ratio calculations for SignalP
- **`src/analyze_go_enrichment.R`**: Custom GO enrichment analysis
- **`src/plot_rna_map.py`**: RNA binding protein positional enrichment visualization

---

## Performance Optimization

- **Caching**: All TxDb and BioMart results are cached in `data/cache/` to improve performance
- **Parallel Processing**: Many R scripts support multi-core execution via `BiocParallel`
- **Logging**: Detailed execution logs are saved to `logs/` for debugging

---

## Citation

If you use this pipeline in your research, please cite:

[Citation information pending publication]

---

## Contact and Support

For questions, issues, or contributions, please contact the project maintainers or open an issue in the repository.

---

## License

[Add appropriate license information]
