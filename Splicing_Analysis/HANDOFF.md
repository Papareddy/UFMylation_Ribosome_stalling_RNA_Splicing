# Project Handoff: Splicing Analysis Pipeline

## Project Status (As of Jan 23, 2026)
This repository contains a comprehensive splicing analysis pipeline for Human, Mouse, and Arabidopsis. It performs motif enrichment, functional impact, and subcellular distribution analyses. The pipeline has recently undergone a major reorganization for better performance and consistency.

## Key Features
- **Multi-Species Support**: Human, Mouse, Arabidopsis.
- **11-Step Continuous Pipeline**: Streamlined structure with optimized execution flow.
- **Persistent Caching**: TxDb and Biomart data cached in `data/cache`.
- **Event-Specific dPSI**: SE (0.2), all others (0.1) - hardcoded in `prepare_rmats_data.R`.
- **Directional Analysis**: Optional directional DeepDive via `--direct` flag.
- **Optional miRNA**: miRNA isoform analysis via `--run-mirna` (skipped by default).

## Pipeline Structure (Renumbered Jan 22-23, 2026)

All steps are now continuously numbered 1-11. Step 5 (Functional Impact) has been consolidated into Step 4.

| # | Step Name | Description | Output Directory |
|---|-----------|-------------|------------------|
| 1 | Data Preparation | Pre-filters rMATS results by FDR/dPSI | `step01_data_prep/` |
| 2 | Domain Enrichment | PFAM/InterPro domain enrichment | `step02_domain_enrichment/` |
| 3 | Protein Attributes | IUPRED, Codon Optimality, etc. | `step03_protein_attributes/` |
| 4 | Protein Sequence Impact | Splicing impact and functional annotation | `step04_protein_sequence_impact/` |
| 5 | Frameshift Density | Metagene analysis of frameshifts | `step05_frameshift_density/` |
| 6 | AA Features | Amino acid composition analysis | `step06_aa_features/` |
| 7 | Biophysical Properties | Charge, Hydropathy, Disorder | `step07_biophysical_properties/` |
| 8 | Motif Analysis | AME/ESE motif enrichment | `step08_motif_analysis/` |
| 9 | Mechanism Investigation | Kinetic window vs Blocking hypothesis | `step09_mechanism_investigation/` |
| 10 | RNA Maps | Positional distribution of motifs | `step10_rna_maps/` |
| 11 | Genomic Associations | miRNA, NMD, EJC/PTC, GC content | `step11_genomic_associations/` |

> [!NOTE]
> Step 15 (Subcellular Distribution) remains as an additional specialized step for Human analysis.

## Recent Changes & Optimizations

### Command-Line Flags
- `--direct`: Runs directional motif analysis (DeepDive_dPSI_positive/negative) in Step 8.
- `--run-mirna`: Enables Step 11 miRNA analysis (computationally expensive).
- `--anchored_window`: Sets window size for genomic associations (default 1000bp).

### Optimizations
- **Anchored Stop Codon Analysis**: `src/anchored_stopcodon_intron_density.R` optimized for speed (~3x). Now uses 25bp bins and outputs raw counts.
- **Directory Consolidation**: `step05` (Functional Impact) is now a sub-directory within `step04` (`step04_protein_sequence_impact/annotated/`).

## Machine Learning Pipeline: Deep Feature Hunter
Located in `miscellaneous/machine_learning/`. Classifies UFM1-dependent vs independent introns using Random Forest.
- **v3 Features (53)**: Includes RBP densities, sequence composition, and **Contextual/Kinetic** features (Upstream Exon GC, Codon Optimality, Decoy Competition).
- **Key Findings**: Kinetic Window Hypothesis supported; SRSF3/PCBP2 density are top discriminators.

## Summary of Completed Runs

### Full Mammalian Pipeline (RI events, FDR 0.05)
**Execution Date**: Jan 23, 2026
- **Human**: `mammalian_RI_dpsi01_fdr05/human/nucleus/` (587 dependent, 345 independent)
- **Mouse**: `mammalian_RI_dpsi01_fdr05/mouse/total/` (207 dependent, 150 independent)

## Maintenance Protocol
1. **Start**: Read this `HANDOFF.md`.
2. **Execution**: Use `conda run -n splicing-functional python run_pipeline.py`.
3. **Location**: `run_pipeline.py` is in the root; all other scripts are in `src/`.
4. **Rules**: Follow `GEMINI.md` for execution order and directory rules.
5. **End**: Update this file with new changes.
