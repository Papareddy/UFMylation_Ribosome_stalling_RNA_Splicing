# Walkthrough: SpliceImpactR Integration and Domain Enrichment Analysis

I have successfully integrated `SpliceImpactR` logic into the splicing analysis pipeline and performed a comprehensive domain enrichment analysis.

## 1. Feature Integration
I modified the pipeline to include two new steps:
- **Step 0.3**: `src/get_splice_impact_features.R`
    - Parses rMATS output (`lost.tsv`, `preserved.tsv`).
    - Annotates splicing events with:
        - **CDS Overlap**: Determines if the event affects the protein-coding sequence.
        - **Domain Mapping**: Fetches Pfam/InterPro domains from Ensembl via `biomaRt`.
        - **Biophysical Features**: Fetches SignalP, TMHMM, and NCOILS annotations.
    - Performs **Enrichment Analysis**: Calculates statistical enrichment of domains and biophysical features relative to the genomic background.
- **Step 0.4**: `src/analyze_domain_enrichment.py`
    - Visualizes the results.

## 2. Issues Resolved
During integration, I addressed several challenges:
- **Dependency Hell**: Manually resolved missing R dependencies (`pwalign`, `txdbmaker`, `ComplexUpset`).
- **Coordinate Mismatch**: Implemented automatic detection and conversion between UCSC (`chr1`) and Ensembl (`1`) chromosome naming styles to ensure correct overlaps.
- **Mixed Event Types**: Updated logic to robustly handle different rMATS event types (SE, A3SS, etc.) seamlessly.

## 3. Results Generated

### Files
The following files are available in `results/human/nucleus/`:
- **Annotated Data**:
    - `lost_with_domains.tsv`: Contains `affects_cds` (TRUE/FALSE) and `impact_domains` (comma-separated list).
    - `preserved_with_domains.tsv`: Same for the preserved set.
- **Statistics**:
    - `domain_enrichment.tsv`: Full statistical table (P-value, FDR, Log2FC) for all domains.
    - `biophysical_enrichment.tsv`: Enrichment statistics for SignalP, TMHMM, and NCOILS.

### Visualizations

#### Biophysical Feature Enrichment
A grouped bar chart showing the enrichment/depletion (Odds Ratio) of Signal Peptides, Transmembrane Helices, and Coiled-Coils in Lost/Preserved sets vs **rMATS Tested Genes** (Background).
![Biophysical Enrichment](/Users/ranjith.papareddy/.gemini/antigravity/brain/08812b67-8c98-46b1-a143-264d44a43a50/Biophysical_Enrichment.png)



#### Domain Enrichment Landscape (Volcano Plot)
A side-by-side comparison of domain enrichment in Lost vs Preserved sets relative to the genome.
- **Red points**: Statistically significant enrichment/depletion (FDR < 0.01).
- **Parameters**: Excluded RI events; FDR Threshold = 0.01.
![Volcano Comparison](/Users/ranjith.papareddy/.gemini/antigravity/brain/08812b67-8c98-46b1-a143-264d44a43a50/Volcano_Domain_Enrichment_Comparison.png)

### Mouse Analysis Results
I also ran the pipeline for **Mouse** (using `mmusculus_gene_ensembl`).

#### Mouse Biophysical Enrichment
![Mouse Biophysical Enrichment](/Users/ranjith.papareddy/.gemini/antigravity/brain/08812b67-8c98-46b1-a143-264d44a43a50/Mouse_Biophysical_Enrichment.png)

#### Mouse Volcano Plot
![Mouse Volcano Comparison](/Users/ranjith.papareddy/.gemini/antigravity/brain/08812b67-8c98-46b1-a143-264d44a43a50/Mouse_Volcano_Domain_Enrichment_Comparison.png)

### Directional Analysis (Human)
Analysis split by splicing direction (Inclusion vs Exclusion).
- **Inc**: Increased Inclusion (dPSI > 0)
- **Exc**: Increased Exclusion (dPSI < 0)

#### Directional Biophysical Enrichment
![Directional Biophysical](/Users/ranjith.papareddy/.gemini/antigravity/brain/08812b67-8c98-46b1-a143-264d44a43a50/Directional_Biophysical_Enrichment.png)

#### Directional Volcano Plots
![Directional Volcano](/Users/ranjith.papareddy/.gemini/antigravity/brain/08812b67-8c98-46b1-a143-264d44a43a50/Directional_Volcano_Domain_Enrichment_Comparison.png)

## 4. Verification
I verified the pipeline execution:
```bash
mamba run -n splicing-functional python3 run_pipeline.py --species human --fraction nucleus --no-fasta
```
- Steps 0.3 and 0.4 complete successfully.
- Output files are populated with valid data (verified >2000 overlaps found).
