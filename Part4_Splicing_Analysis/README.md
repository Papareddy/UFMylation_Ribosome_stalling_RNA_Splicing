# 🧬 UFM1 Splicing Functional Impact Pipeline

This repository contains an advanced analytical pipeline for defining the **UFM1-dependent splicing program** and its evolutionary conservation across the eukaryotic kingdom (**Human**, **Mouse**, **Arabidopsis**).

## 🌍 Key Focus: The Eukaryotic Conservation Hub
The project has evolved from a mammalian analysis into a unified **Eukaryotic Hub**, revealing a fundamental split in the cellular response to ribosome stalling:
1.  **UFM1-Dependent ("ER Architects")**: Splicing specifically remodels the ER environment and lipid membranes (e.g., `ESYT1`, `TAFAZZIN`).
2.  **UFM1-Independent ("Core OS")**: Housekeeping maintenance (e.g., **RNA Splicing factors**) is universally conserved but independent of UFM1 stress.
3.  **Signal Peptide Eraser**: UFM1-dependent splicing acts as a functional switch with **100% penetrance**—systematically stripping Signal Peptides from targeted transcripts.

---

## 🚀 Execution Guide

### 1. Environment Setup
The pipeline uses a multi-language (Python/R) environment.
```bash
mamba env create -f splicing-functional.yml
mamba activate splicing-functional
```

### 2. Main Mammalian Pipeline
Run the 15-step pipeline for individual species:
```bash
python run_pipeline.py --species human
python run_pipeline.py --species mouse
```

### 3. Eukaryotic Conservation Hub (Batch Execution)
To regenerate the cross-kingdom integrated results (Human + Mouse + Arabidopsis):
```bash
# This script automates GO enrichment, metadata mapping, and integrated landscape plotting
bash src/run_GO_Final_pipeline.sh
```

---

## 📂 Pipeline Overview

| Step | Analysis Phase | Description | Key Artifacts |
| :--- | :--- | :--- | :--- |
| **1-6** | **Impact Analysis** | Classification (Lost/Preserved), Frameshifts, NMD, Start/Stop codon disruption. | `impact_fractions.png`, `Alignment_Scores_Histogram.png` |
| **8-10** | **Sequence Feature** | GC content, MaxEnt scores, and AME Motif Enrichment (CisBP/RBPDB). | `combined_features.pdf`, `motif_comparison.png` |
| **11-12** | **Mechanistic** | Translation stalling metrics & Positional RNA Maps of RNA-Binding Proteins. | `stalling_stats.tsv`, `RNA_Map_*.png` |
| **14** | **Eukaryotic Hub** | Cross-species GO enrichment & Conservation Landscapes. | `Eukaryotic_Conserved_Landscape.png` |
| **15** | **Signal Peptide** | Quantitative analysis of Signal Peptide loss in RI transcripts. | `SignalP_Loss_Analysis.tsv` |

---

## 📈 Major Scientific Results

### 1. Conserved ER Sectors
Our 3-species landscape identifies **ER Organization** and **Lipid Metabolism** as the universal eukaryotic targets of UFM1-dependent splicing. This suggests UFM1-mediated splicing is an ancient mechanism for protecting the secretory factory.

### 2. Signal Peptide Erasing (100% Penetrance)
We found that **100%** of Retained Intron (RI) isoforms derived from SP-positive UFM1-dependent genes in Human and Mouse lose their signal peptide. This confirms splicing is used to systematically reroute proteins away from the ER during stress.

### 3. Splicing Factors Autoregulation
Baseline splicing machinery (e.g., `SF3B1`, `SRSF4`) is a conserved target in the UFM1-independent group, likely representing a universal homeostatic "emergency brake" on RNA processing during ribosome stalling.

---

## 🛠️ Manual Tools & Data Visualization
Standalone scripts are provided in `src/` for custom analysis:
*   **[`src/user_plot_eukaryotic_landscape.R`](file:///Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/src/user_plot_eukaryotic_landscape.R)**: Customizable bubble plots for 3-species conservation.
*   **[`src/analyze_signalp_loss.R`](file:///Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/src/analyze_signalp_loss.R)**: Background-controlled Odds Ratio calculations for SignalP.

---
**Maintenance**: All TxDb and BioMart results are cached in `data/cache/` to ensure performance.
