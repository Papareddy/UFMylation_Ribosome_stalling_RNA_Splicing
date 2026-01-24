#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT: analyze_signalp_loss.R
# PURPOSE: Quantify Signal Peptide presence in UFM1-Dep/Ind genes and check for loss in RI isoforms.
# ==============================================================================

suppressPackageStartupMessages({
  library(biomaRt)
  library(tidyverse)
  library(argparse)
})

# --- Arguments ---
parser <- ArgumentParser(description='Signal Peptide Loss Analysis')
parser$add_argument('--species', required=TRUE, help='Species (human, mouse, arabidopsis)')
parser$add_argument('--rds', required=TRUE, help='Path to UFM1_events_rich.rds')
parser$add_argument('--outdir', required=TRUE, help='Output directory')
args <- parser$parse_args()

SPECIES <- tolower(args$species)
RDS_PATH <- args$rds
OUTDIR <- args$outdir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# --- Configuration ---
SPECIES_META <- list(
  human = list(
    mart_name = "ensembl", mart_data = "hsapiens_gene_ensembl", mart_host = "https://www.ensembl.org",
    id_type = "ensembl_gene_id", tx_id = "ensembl_transcript_id", strip_version = TRUE
  ),
  mouse = list(
    mart_name = "ensembl", mart_data = "mmusculus_gene_ensembl", mart_host = "https://www.ensembl.org",
    id_type = "ensembl_gene_id", tx_id = "ensembl_transcript_id", strip_version = TRUE
  ),
  arabidopsis = list(
    mart_name = "plants_mart", mart_data = "athaliana_eg_gene", mart_host = "https://plants.ensembl.org",
    id_type = "ensembl_gene_id", tx_id = "ensembl_transcript_id", strip_version = FALSE
  )
)

if (!(SPECIES %in% names(SPECIES_META))) stop("Unsupported species: ", SPECIES)
meta <- SPECIES_META[[SPECIES]]

# --- Load Genes ---
get_gene_lists <- function(rds_path, strip_version) {
  if(!file.exists(rds_path)) return(NULL)
  events <- readRDS(rds_path)
  
  # Group RI events (Dependent vs Independent)
  ri <- events[events$EventType == "RI",]
  genes_dep   <- unique(ri$GeneID[ri$Group == "UFM1_dependent"])
  genes_indep <- unique(ri$GeneID[ri$Group == "UFM1_independent"])
  
  # Also Background (All analyzed genes)
  genes_bg <- unique(ri$GeneID)
  
  if(strip_version) {
    genes_dep   <- sub("\\..*", "", genes_dep)
    genes_indep <- sub("\\..*", "", genes_indep)
    genes_bg    <- sub("\\..*", "", genes_bg)
  }
  
  return(list(dep = genes_dep, ind = genes_indep, bg = genes_bg))
}

message(">>> Loading Gene Lists...")
lists <- get_gene_lists(RDS_PATH, meta$strip_version)
all_genes <- unique(c(lists$dep, lists$ind))

message("Found Genes: Dep=", length(lists$dep), " Indep=", length(lists$ind))

# --- BioMart Query ---
message(">>> Connecting to BioMart...")
mart <- useMart(meta$mart_name, dataset = meta$mart_data, host = meta$mart_host)

# 1. Fetch SignalP for ALL genes in the genome to create a background baseline
message("Fetching genome-wide SignalP status...")
# We use lists$bg as the "universe" of analyzed genes if available, 
# otherwise we fetch for all event genes + background sample.
genome_bg_data <- getBM(attributes = c(meta$id_type, "signalp_start"), 
                        mart = mart)
colnames(genome_bg_data) <- c("GeneID", "SignalP_Start")
genome_bg_data$Has_SignalP <- !is.na(genome_bg_data$SignalP_Start) & genome_bg_data$SignalP_Start != ""

# Clean up genome data to gene level
gene_sp_status_bg <- genome_bg_data %>%
  group_by(GeneID) %>%
  summarise(Has_SP_Gene = any(Has_SignalP))

# 2. Fetch specific transcript details (Biotype) only for our event genes
message("Fetching transcript biotypes for event genes...")
event_bm_data <- getBM(attributes = c(meta$id_type, meta$tx_id, "transcript_biotype", "signalp_start"), 
                       filters = meta$id_type, 
                       values = all_genes, 
                       mart = mart)
colnames(event_bm_data) <- c("GeneID", "TxID", "Biotype", "SignalP_Start")
event_bm_data$Has_SignalP <- !is.na(event_bm_data$SignalP_Start) & event_bm_data$SignalP_Start != ""

# --- Analysis 1: Gene-Level SignalP Enrichment (vs Genome) ---

calc_enrichment_vs_genome <- function(target_ids, bg_df, label) {
  # Target stats
  target_df <- bg_df %>% filter(GeneID %in% target_ids)
  a <- sum(target_df$Has_SP_Gene)
  b <- nrow(target_df) - a
  
  # Non-target (Genome background) stats
  nontarget_df <- bg_df %>% filter(!(GeneID %in% target_ids))
  c <- sum(nontarget_df$Has_SP_Gene)
  d <- nrow(nontarget_df) - c
  
  # Fisher Test
  mat <- matrix(c(a, c, b, d), nrow=2)
  ft <- fisher.test(mat)
  
  return(data.frame(
    Group = label, 
    SP_Pos = a, 
    Total = nrow(target_df), 
    Rate = a / nrow(target_df),
    Genome_SP_Pos = c,
    Genome_Total = nrow(nontarget_df),
    Genome_Rate = c / nrow(nontarget_df),
    OddsRatio = ft$estimate,
    Log2OR = log2(ft$estimate),
    Pval = ft$p.value
  ))
}

res_gene <- rbind(
  calc_enrichment_vs_genome(lists$dep, gene_sp_status_bg, "Dependent"),
  calc_enrichment_vs_genome(lists$ind, gene_sp_status_bg, "Independent")
)

message(">>> Gene-Level SignalP Enrichment (vs Genome):")
print(res_gene)
write.table(res_gene, file.path(OUTDIR, "SignalP_Gene_Enrichment_vs_Genome.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# --- Analysis 2: Transcript-Level Loss (RI vs Canonical) ---
# For genes that HAVE a Signal Peptide (in some transcript), do they have an RI transcript that LOSES it?
# Definition of RI transcript: Biotype = "retained_intron"
# Definition of Canonical: Biotype = "protein_coding" (and Has SignalP)

# Filter for SP+ Genes Only (based on genome-wide list)
sp_pos_genes <- gene_sp_status_bg$GeneID[gene_sp_status_bg$Has_SP_Gene]
sp_pos_event_data <- event_bm_data %>% filter(GeneID %in% sp_pos_genes)

# Analyze Loss
analyze_loss <- function(group_genes, group_label) {
  # Genes in this group that are SP+
  target_genes <- intersect(group_genes, sp_pos_genes)
  
  loss_events <- 0
  total_ri_tx <- 0
  genes_with_loss <- 0
  
  for(g in target_genes) {
    sub <- sp_pos_event_data %>% filter(GeneID == g)
    
    # Check if has 'retained_intron' biotype transcript
    ri_txs <- sub[sub$Biotype == "retained_intron",]
    
    if(nrow(ri_txs) > 0) {
      # We count how many RI transcripts exist
      total_ri_tx <- total_ri_tx + nrow(ri_txs)
      
      # How many lack SignalP?
      loss_count <- sum(!ri_txs$Has_SignalP)
      loss_events <- loss_events + loss_count
      
      if(loss_count > 0) genes_with_loss <- genes_with_loss + 1
    }
  }
  
  return(data.frame(Group=group_label,  
                    SP_Pos_Genes=length(target_genes),
                    Genes_With_RI_Tx=genes_with_loss, 
                    RI_Tx_Signal_Lost=loss_events,
                    Total_RI_Tx=total_ri_tx,
                    Loss_Rate = ifelse(total_ri_tx>0, loss_events/total_ri_tx, 0)))
}

res_loss <- rbind(
  analyze_loss(lists$dep, "Dependent"),
  analyze_loss(lists$ind, "Independent")
)

message(">>> Transcript-Level Signal Loss (RI isoforms):")
print(res_loss)
write.table(res_loss, file.path(OUTDIR, "SignalP_Loss_Analysis.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# --- Save Raw Data ---
write.table(event_bm_data, file.path(OUTDIR, "BioMart_Transcripts_SignalP_EventGenes.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
message("Done.")
