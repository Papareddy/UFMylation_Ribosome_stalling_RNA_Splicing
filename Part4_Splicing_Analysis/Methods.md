# Methods

## Defining the UFM1-Dependent Splicing Program and Its Evolutionary Conservation

### Reference Genomes and Annotations

For human, we utilized the GRCh38 primary assembly genome and GENCODE v45 annotations, restricted to protein-coding genes. For mouse, we used the GRCm39 assembly with Ensembl Release 112 annotations. For *Arabidopsis thaliana*, we employed the TAIR10 reference genome and Ensembl Plants Release 56 annotations, filtered to exclude plastid and ribosomal RNA genes. Genome sequence and annotation files were obtained from Ensembl^1^ and GENCODE^2^.

### Differential Splicing Analysis

Differential splicing between control and anisomycin-treated conditions was quantified using rMATS v4.1.2^3^ with default parameters. Events were retained if they met the following criteria: false discovery rate (FDR) < 0.05, absolute delta percent-spliced-in (|ΔPSI|) ≥ 0.1 for retained intron (RI) events or ≥ 0.2 for skipped exon (SE) events, and a minimum of 10 junction reads supporting the event. Splicing events were classified as UFM1-dependent if they exhibited significant changes in wild-type cells but not in UFM1 knockout cells, and UFM1-independent if changes occurred in both genetic backgrounds.

### Splice Site Strength Scoring

Splice site strength was evaluated using the Maximum Entropy Model implemented in MaxEntScan^4^, which scores 9-nt 5' splice sites and 23-nt 3' splice sites based on their adherence to the position weight matrix derived from human constitutive splice sites. For each splicing event, sequences flanking the annotated splice junctions were extracted from the reference genome using the Biopython library^5^ and scored using the Python port of MaxEntScan. Higher scores indicate stronger splice site consensus.

### Motif Enrichment Analysis

De novo and known motif enrichment was performed using the Analysis of Motif Enrichment (AME) tool from the MEME Suite v5.x^6^. For retained intron events, intronic sequences were extracted and compared against constitutive intron background sequences using Fisher's exact test with E-value threshold of 0.05. Position frequency matrices for RNA-binding protein motifs were obtained from the CisBP-RNA database^7^ and custom curated RBPDB motif sets. GC content was calculated for each sequence to control for compositional biases.

### Protein Domain and Functional Impact Analysis

To assess the functional consequences of alternative splicing, we mapped splicing events to protein domain annotations using biomaRt^8^ to query Ensembl proteome data and InterPro domain coordinates. For each event, we determined whether the alternative region overlapped with annotated functional domains and classified events by predicted functional impact: frame-preserving (multiple of 3 nucleotides), frameshift-inducing, or nonsense-mediated decay (NMD)-triggering based on the position of premature termination codons relative to the last exon-exon junction. Reading frame analysis was performed using custom Python scripts leveraging the pandas^9^ and NumPy^10^ libraries for data manipulation.

### Signal Peptide Analysis

Signal peptide presence and loss was predicted using SignalP 6.0^11^, which employs deep neural network architectures for sequence-based signal peptide prediction. For each gene with retained intron events, we extracted the translated protein sequences of both the canonical transcript and the predicted RI-containing isoform. Signal peptide loss was defined as the presence of a high-confidence signal peptide (probability > 0.5) in the canonical isoform but not in the alternative isoform. Odds ratios comparing signal peptide loss between UFM1-dependent and UFM1-independent events were calculated using Fisher's exact test implemented in SciPy^12^.

### Gene Ontology Enrichment Analysis

Functional enrichment analysis was performed using clusterProfiler v4.0^13^ in R. For each category of splicing events (UFM1-dependent, UFM1-independent) and direction (increased or decreased inclusion), we tested for over-representation of Gene Ontology (GO) Biological Process, Cellular Component, and Molecular Function terms. Background gene sets were defined as all genes expressed above the detection threshold in each species. P-values were corrected for multiple testing using the Benjamini-Hochberg method with FDR < 0.05 considered significant. Cross-species conservation of enriched GO terms was assessed by identifying terms significantly enriched (FDR < 0.05) in at least two of the three species analyzed.

### Genomic Feature Annotation

Genomic coordinates of splicing events were manipulated using BEDTools^14^ and its Python wrapper pybedtools^15^ for interval arithmetic operations. Events were annotated with overlapping genomic features including gene biotype, transcript support level, and chromosomal location. Sequence extraction from reference genomes was performed using SAMtools^16^ faidx functionality. Data visualization was performed using the tidyverse suite^17^ in R, with complex set visualizations generated using ComplexUpset^18^. All statistical analyses were performed in R (≥v4.3) and Python (v3.11), with environment management via Mamba/Conda.

### Cross-Species Conservation Analysis

To identify evolutionarily conserved splicing responses, we performed ortholog mapping between human, mouse, and *Arabidopsis* using Ensembl Compara homology data accessed via biomaRt. For each species, genes undergoing UFM1-dependent or UFM1-independent splicing were mapped to their orthologs, and we identified genes exhibiting conserved splicing patterns across species. GO term enrichment was independently performed for each species, and terms reaching significance in multiple species were designated as conserved functional modules.

### Data Visualization and Statistical Analysis

All statistical tests were two-sided unless otherwise specified. Multiple testing correction was applied using the Benjamini-Hochberg procedure. Box plots display median (center line), interquartile range (box limits), and 1.5× interquartile range (whiskers). Violin plots show kernel density estimates of data distribution. Heatmaps were generated using hierarchical clustering with complete linkage. All computational analyses were implemented using custom Python and R scripts, with reproducible execution managed through a unified pipeline framework.

---

## References

1. Cunningham, F. *et al.* Ensembl 2022. *Nucleic Acids Res.* **50**, D988–D995 (2022).
2. Frankish, A. *et al.* GENCODE 2021. *Nucleic Acids Res.* **49**, D916–D923 (2021).
3. Shen, S. *et al.* rMATS: robust and flexible detection of differential alternative splicing from replicate RNA-Seq data. *Proc. Natl. Acad. Sci. USA* **111**, E5593–E5601 (2014).
4. Yeo, G. & Burge, C. B. Maximum entropy modeling of short sequence motifs with applications to RNA splicing signals. *J. Comput. Biol.* **11**, 377–394 (2004).
5. Cock, P. J. A. *et al.* Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* **25**, 1422–1423 (2009).
6. McLeay, R. C. & Bailey, T. L. Motif enrichment analysis: a unified framework and an evaluation on ChIP data. *Bioinformatics* **26**, 46–54 (2010).
7. Ray, D. *et al.* A compendium of RNA-binding motifs for decoding gene regulation. *Nature* **499**, 172–177 (2013).
8. Durinck, S. *et al.* BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis. *Bioinformatics* **21**, 3439–3440 (2005).
9. McKinney, W. Data structures for statistical computing in Python. *Proc. 9th Python Sci. Conf.* 56–61 (2010).
10. Harris, C. R. *et al.* Array programming with NumPy. *Nature* **585**, 357–362 (2020).
11. Teufel, F. *et al.* SignalP 6.0 predicts all five types of signal peptides using protein language models. *Nat. Biotechnol.* **40**, 1023–1025 (2022).
12. Virtanen, P. *et al.* SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat. Methods* **17**, 261–272 (2020).
13. Wu, T. *et al.* clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *Innovation* **2**, 100141 (2021).
14. Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* **26**, 841–842 (2010).
15. Dale, R. K., Pedersen, B. S. & Quinlan, A. R. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. *Bioinformatics* **27**, 3423–3424 (2011).
16. Danecek, P. *et al.* Twelve years of SAMtools and BCFtools. *GigaScience* **10**, giab008 (2021).
17. Wickham, H. *et al.* Welcome to the Tidyverse. *J. Open Source Softw.* **4**, 1686 (2019).
18. Krassowski, M. ComplexUpset: a Python package for complex upset plots. Zenodo (2020).
