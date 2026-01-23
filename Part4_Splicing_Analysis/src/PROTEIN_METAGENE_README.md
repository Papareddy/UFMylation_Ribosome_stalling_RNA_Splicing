# Protein Property Metagene Analysis

## Overview
This script calculates protein biochemical properties using Bio.SeqUtils.ProtParam and creates metagene plots showing how properties vary along protein length.

## Features

### Protein Properties Calculated (via ProtParam)
1. **GRAVY** - Hydropathy (Grand Average of Hydropathicity)
2. **Aromaticity** - Fraction of aromatic amino acids
3. **Instability Index** - Protein stability prediction
4. **Isoelectric Point** - pI
5. **Charge at pH 7** - Net charge
6. **Helix Fraction** - Predicted α-helix content
7. **Sheet Fraction** - Predicted β-sheet content
8. **Aliphatic Index** - Thermostability indicator

### Analysis Approach
- **Full-length proteins** extracted (longest CDS isoform per gene)
- **Sliding window** analysis (default: 30 AA window, 10 AA step)
- **Normalized to 0-100%** protein length (metagene)
- **Three groups**: Control, UFM1-dependent, UFM1-independent

## Usage

```bash
python src/protein_properties_metagene.py \
  --gtf data/arabidopsis/no_plastid_no_rRNA.Arabidopsis_thaliana.TAIR10.56.gtf \
  --genome data/arabidopsis/Arabidopsis_thaliana_TAIR10.dna.primary_assembly.fa \
  --dep_genes dep_gene_list.txt \
  --indep_genes indep_gene_list.txt \
  --outdir results/protein_metagene \
  --window 30 \
  --step 10 \
  --nbins 100
```

## Input Files

### Gene Lists (create these first)
```bash
# Extract gene lists from your events
awk -F'\t' '{print $2}' results/arabidopsis/nucleus/step01_data_prep/UFM1_dependent.tsv | tail -n +2 > dep_genes.txt
awk -F'\t' '{print $2}' results/arabidopsis/nucleus/step01_data_prep/UFM1_independent.tsv | tail -n +2 > indep_genes.txt
```

## Outputs

### Plots
- `protein_properties_metagene.pdf` - 8-panel figure showing all properties
- `protein_properties_metagene.png` - High-res PNG version

### Data Tables
- `metagene_Control.tsv` - Control protein properties by position
- `metagene_Dependent.tsv` - UFM1-dependent properties
- `metagene_Independent.tsv` - UFM1-independent properties

Each table contains:
- Position_pct (0-100%)
- Property values at each position
- Property_SEM (standard error of mean)

## Parameters

- `--window` - Sliding window size in amino acids (default: 30)
  - Larger = smoother curves, less resolution
  - Smaller = more detail, noisier

- `--step` - Step size for sliding window (default: 10)
  - Smaller = more datapoints, slower

- `--nbins` - Number of bins for normalized metagene (default: 100)
  - Standard for metagene plots

## Dependencies

```bash
conda install -c conda-forge biopython matplotlib seaborn pandas numpy
```

## Example Analysis Pipeline

```bash
# 1. Extract gene lists
awk -F'\t' '{print $2}' mammalian_RI_dpsi01_fdr05/human/nucleus/step01_data_prep/UFM1_dependent.tsv | tail -n +2 | sort -u > dep_genes.txt
awk -F'\t' '{print $2}' mammalian_RI_dpsi01_fdr05/human/nucleus/step01_data_prep/UFM1_independent.tsv | tail -n +2 | sort -u > indep_genes.txt

# 2. Run analysis
python src/protein_properties_metagene.py \
  --gtf data/human/pcg_gencode.v45.annotation.gtf.gz \
  --genome data/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --dep_genes dep_genes.txt \
  --indep_genes indep_genes.txt \
  --outdir results/protein_metagene_human

# 3. View results
open results/protein_metagene_human/protein_properties_metagene.pdf
```

## Interpretation

### What the plots show:
- **X-axis**: Relative position in protein (0% = N-terminus, 100% = C-terminus)
- **Y-axis**: Property value
- **Lines**: Mean value across all proteins in group
- **Shaded areas**: Standard error (uncertainty)

### Expected patterns:
- **Hydropathy (GRAVY)**: May vary if proteins have different membrane association
- **Charge**: Often higher at termini (localization signals)
- **Instability**: May reveal protein stability differences
- **Secondary structure**: Structural preferences along protein length

## Performance

- ~1-2 minutes for 500 genes
- ~5-10 minutes for 2000 genes
- Scales linearly with number of genes

## Troubleshooting

**"No module named 'Bio'"**
```bash
conda install -c conda-forge biopython
```

**"Too few proteins extracted"**
- Check GTF gene_id format matches gene lists
- Verify genome chromosome names match GTF

**"All properties are NaN"**
- Check that proteins contain standard amino acids
- Verify window size isn't larger than shortest protein
