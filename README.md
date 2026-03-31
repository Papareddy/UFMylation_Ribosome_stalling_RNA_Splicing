# UFM1-Mediated Splicing Regulation in Response to Ribosome Stalling

This repository provides a comprehensive multi-omic analytical framework for studying the UFM1-dependent retrograde signaling pathway that links ER-associated translation stress to nuclear mRNA splicing.

## Scientific Overview

Our research uncovers a non-canonical mechanism where ER-bound ribosomes, under stalling conditions, serve as a dynamic docking platform. This platform captures SR shuttling proteins, establishing a spatial regulatory axis between the Endoplasmic Reticulum and the Nucleus. This UFMylation-driven process couples translational stress to nuclear mRNA processing, selectively reshaping the expression of membrane-associated genes and potentially driving membrane remodeling.

This repository integrates phylogenetic profiling, quantitative proteomics, and cross-species transcriptomics across **Human**, **Mouse**, and **Arabidopsis**.

---

## Repository Structure

The project is organized into five core analysis modules, each corresponding to specific phases of the study:

### 1. Phylogenetic Profiling (`phylo_profiling/`)
*   **Purpose**: Identification of proteins that have co-evolved with the UFM1 system.
*   **Key Files**: MCL clusters, UFM1-correlated orthogroups, and functional GO annotations.

### 2. Ribosome Proteomics (`Ribosome_Proteomics_Analysis/`)
*   **Purpose**: Quantitative analysis of ribosome-associated proteomes under stalling stress.
*   **Outputs**: Differential expression (volcano plots) and GO enrichment for stalled ribosome clusters (Figure 3).

### 3. Subcellular Enrichment (`Subcellular_Enrichmnet_Analysis/`)
*   **Purpose**: Subcellular relocalisation of proteins under ribosome stalling in UFM1 dependent manner.
*   **Methods**: Statistical enrichment testing and visualization of fraction-specific targets.

### 4. Splicing Functional Impact (`Splicing_Analysis/`)
*   **Purpose**: A multi-species pipeline to quantify the functional consequences of UFM1-dependent alternative splicing.
*   **Features**: Frame-shift density analysis, protein domain enrichment, signal peptide loss quantification, and cross-species (Eukaryotic-wide) consolidation.

### 5. Supplemental Data (`Tables_S1-S9/`)
*   **Purpose**: Consolidated datasets supporting the manuscript's findings, including orthology mappings and cross-species splicing events.

---

## Pipeline Overview

This analysis integrates phylogenetic profiling, quantitative proteomics, and transcriptomics across human, mouse, and Arabidopsis models to quantify protein and transcript sequestering at the ER vs. nuclear localization.

---

## Contact & Citation

**Author**: Ranjith K. Papareddy

**DOI**: https://doi.org/10.64898/2026.03.30.715226

For technical inquiries regarding the pipelines, please refer to the README files within each specific sub-module directory.
