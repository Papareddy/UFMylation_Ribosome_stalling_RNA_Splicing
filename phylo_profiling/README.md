# Phylogenetic Profiling: UFM1 Co-evolution Analysis

This repository contains supplemental datasets, custom scripts, and analysis results for the phylogenetic profiling study identifying proteins that co-evolved with **UFM1**. 

Ref: *Zhan, Papareddy et al.*

---

## Directory Structure

### 1. Datasets
Contains the house-predicted proteomes and associated metadata used for analysis.
*   `eukaryotes.fasta`: Combined FASTA file of all predicted proteomes.
*   `eukaryotes_metadata.tab`: Tab-delimited metadata including source and accession numbers.

### 2. Orthogroups
Contains the orthology clusters (MCL clusters) forming the basis for profiling.
*   `dump.refined_mcl_clusters.filt_expanded`: Filtered and expanded orthogroup file.
*   `UFM1.orthogroup_assignment.tab`: Mapping of proteins to their respective clusters.

### 3. Phylogenetic Profiling
Execution results from the profiling analysis.
*   `UFM1.orthogroup_assignment.tab`: Detailed profiling results (p-values, significance class, PBC correlation, etc.).
*   `UFM1.clade_file.txt`: Clade definitions used to categorize species.

### 4. Annotations
Functional Gene Ontology (GO) annotations based on **Wei2Go**.
*   `dump.refined_mcl_clusters.filt_expanded.go_annotation`: Background GO annotations for all clusters.
*   `UFM1.correlated.go_annotation`: GO annotations for UFM1-correlated clusters.
*   `UFM1.correlated.go_annotation.go_enrichment.tab`: Results of the GO enrichment analysis (fold enrichment, p-values).

### 5. Scripts
Custom Python implementations for the profiling pipeline.
*   `profiling.py`: Script for phylogenetic profiling calculation.
*   `go_enrichment.py`: Script for GO-enrichment analysis execution.

---

## Requirements & Usage

### Language
*   **Python 3**

### Dependencies
*   `numpy`
*   `pandas`
*   `scipy`
*   `multiprocess`

> [!NOTE]
> Usage instructions and example commands are provided within the header of each script in the `5.Scripts` directory.

---

## Contact Information

*   **Author:** Nick Irwin
*   **Email:** nicholas.irwin@gmi.oeaw.ac.at
*   **Institution:** Gregor Mendel Institute (GMI), Vienna BioCenter, Vienna, AT
