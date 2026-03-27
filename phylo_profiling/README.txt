# README: Supplemental scripts and Data for Zhan, Papareddy et al.

## Project Overview
This supplemental dataset provides the proteomes, orthology clusters, results, and custom scripts for the phylogenetic profiling analyses used to identify proteins co-evolving with UFM1.

---

## Directory Structure

### 1. Datasets
This directory contains the predicted proteomes used in the analysis and associated metadata.
*   `eukaryotes.fasta`: A combined FASTA file containing all predicted proteomes used in the study.
*   `eukaryotes_metadata.tab`: A tab-delimited table containing the source and accession numbers for each predicted proteome.

### 2. Orthogroups
This directory contains the orthogroups (clusters) on which the profiling was based.
*   `dump.refined_mcl_clusters.filt_expanded`: Orthogroup file containing each cluster and its associated proteins.
*   `UFM1.orthogroup_assignment.tab`: A list mapping each protein in the dataset to its associated cluster.

### 3. Phylogenetic_profiling
Data resulting from the phylogenetic profiling analysis.
*   `UFM1.orthogroup_assignment.tab`: Profiling results including clusters, significance class, p-values, median/average homolog numbers for both classes, consistency (% of species in the class with the protein), and the PBC correlation.
*   `UFM1.clade_file.txt`: The clade definitions used to categorize species in the profiling analysis.

### 4. Annotations
Gene Ontology (GO) annotations for each cluster based on Wei2Go.
*   `dump.refined_mcl_clusters.filt_expanded.go_annotation`: Every cluster and its associated GO terms (background set).
*   `UFM1.correlated.go_annotation`: UFM1-correlated clusters and their associated GO terms (test set).
*   `UFM1.correlated.go_annotation.go_enrichment.tab`: Results of the GO-enrichment analysis including GO IDs, p-values, fold enrichment, observed counts (test vs. background), and expected counts.

### 5. Scripts
Custom Python code used to run the profiling and GO-enrichment analyses. Usage and example instructions are provided within the header of each script.
*   `profiling.py`: Script for conducting the phylogenetic profiling analysis.
*   `go_enrichment.py`: Script for conducting the GO-enrichment analysis.

---

## Requirements
*   **Language:** Python 3
*   **Dependencies:** numpy, pandas, scipy, multiprocess

## Contact Information
*   **Author:** Nick Irwin
*   **Email:** nicholas.irwin@gmi.oeaw.ac.at
*   **Institution:** Gregor Mendel Institute (GMI), Vienna BioCenter, Vienna, AT
