#!/bin/bash

# ==============================================================================
# Master Pipeline: GO Enrichment & Cross-Species Conservation (Final Test)
# ==============================================================================

BASE_OUT="GO_Final_test"
mkdir -p $BASE_OUT

# Define Inputs
HUMAN_RDS="mammalian_RI_dpsi01_fdr05/human/nucleus/step01_data_prep/UFM1_events_rich.rds"
MOUSE_RDS="mammalian_RI_dpsi01_fdr05/mouse/total/step01_data_prep/UFM1_events_rich.rds"
ARAB_RDS="mammalian_RI_dpsi01_fdr05/arabidopsis/nucleus/step01_data_prep/UFM1_events_rich.rds"

# ------------------------------------------------------------------------------
# STEP 1: Species-Specific Enrichment (Steps 13 & 14)
# ------------------------------------------------------------------------------

# Human
echo ">>> Running Human (Nucleus)..."
mamba run -n splicing-functional Rscript src/analyze_subcellular_localization.R --species human --rds $HUMAN_RDS --outdir $BASE_OUT/human/nucleus/step13_subcellular_localization
mamba run -n splicing-functional Rscript src/analyze_go_enrichment.R --species human --rds $HUMAN_RDS --outdir $BASE_OUT/human/nucleus/step14_go_enrichment

# Mouse
echo ">>> Running Mouse (Total)..."
mamba run -n splicing-functional Rscript src/analyze_subcellular_localization.R --species mouse --rds $MOUSE_RDS --outdir $BASE_OUT/mouse/total/step13_subcellular_localization
mamba run -n splicing-functional Rscript src/analyze_go_enrichment.R --species mouse --rds $MOUSE_RDS --outdir $BASE_OUT/mouse/total/step14_go_enrichment

# Arabidopsis
echo ">>> Running Arabidopsis (Nucleus)..."
mamba run -n splicing-functional Rscript src/analyze_subcellular_localization.R --species arabidopsis --rds $ARAB_RDS --outdir $BASE_OUT/arabidopsis/nucleus/step13_subcellular_localization
mamba run -n splicing-functional Rscript src/analyze_go_enrichment.R --species arabidopsis --rds $ARAB_RDS --outdir $BASE_OUT/arabidopsis/nucleus/step14_go_enrichment

# ------------------------------------------------------------------------------
# STEP 2: Eukaryotic Conservation (Human, Mouse, Arabidopsis)
# ------------------------------------------------------------------------------

echo ">>> Running Eukaryotic Conservation Analysis..."

# 1. Venn Diagram (Overlap of Mammalian Dependent GO Terms)
mamba run -n splicing-functional Rscript src/plot_species_go_venn.R

# 2. Eukaryotic Conserved Dependent Landscape (3-species Bubble Plot)
mamba run -n splicing-functional Rscript src/plot_conserved_go_landscape.R

# 3. Eukaryotic Conserved Independent Landscape (3-species Bubble Plot)
mamba run -n splicing-functional Rscript src/plot_conserved_independent_go_landscape.R

# 4. Extract Genes associated with Eukaryotic Conserved Terms (Unified Hub)
mamba run -n splicing-functional Rscript src/extract_eukaryotic_go_genes.R

echo ">>> Pipeline Complete. Results in: $BASE_OUT"
