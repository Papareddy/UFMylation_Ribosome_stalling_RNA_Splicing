#!/usr/bin/env Rscript

# ==============================================================================
# SUBCELLULAR LOCALIZATION ENRICHMENT ANALYSIS (BioMart / GO)
# ==============================================================================
# Analyzes enrichment of specific subcellular compartments in UFM1-dependent 
# and independent gene sets vs Genome Background.
# ==============================================================================

suppressPackageStartupMessages({
  library(biomaRt)
  library(tidyverse)
  library(ggpubr)
  library(GenomicFeatures)
  library(argparse)
})

# --- Arguments ---
parser <- ArgumentParser(description='Subcellular Localization Enrichment')
parser$add_argument('--species', required=TRUE, help='Species (human, mouse, arabidopsis)')
parser$add_argument('--rds', required=TRUE, help='Path to UFM1_events_rich.rds')
parser$add_argument('--outdir', required=TRUE, help='Output directory')
args <- parser$parse_args()

# --- Configuration ---
SPECIES <- tolower(args$species)
RDS_PATH <- args$rds
OUTDIR <- args$outdir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Species Meta
SPECIES_META <- list(
  human = list(
    mart_name = "ensembl", mart_data = "hsapiens_gene_ensembl", mart_host = "https://www.ensembl.org",
    id_type = "ensembl_gene_id", go_attr = "go_id", strip_version = TRUE, bg_tot = 20000
  ),
  mouse = list(
    mart_name = "ensembl", mart_data = "mmusculus_gene_ensembl", mart_host = "https://www.ensembl.org",
    id_type = "ensembl_gene_id", go_attr = "go_id", strip_version = TRUE, bg_tot = 20000
  ),
  arabidopsis = list(
    mart_name = "plants_mart", mart_data = "athaliana_eg_gene", mart_host = "https://plants.ensembl.org",
    id_type = "ensembl_gene_id", go_attr = "go_id", strip_version = FALSE, bg_tot = 27500
  )
)

if (!(SPECIES %in% names(SPECIES_META))) stop("Unsupported species: ", SPECIES)
meta <- SPECIES_META[[SPECIES]]

# Localization GO Terms
LOC_TERMS <- list(
  "Nucleus"      = "GO:0005634",
  "Cytoplasm"    = "GO:0005737",
  "Mitochondria" = "GO:0005739",
  "ER / Network" = "GO:0005783",
  "Secreted"     = c("GO:0005576", "GO:0005615"),
  "Membrane"     = "GO:0016020",
  "SRP"          = "GO:0005786"
)
if(SPECIES == "arabidopsis") LOC_TERMS[["Chloroplast"]] <- "GO:0009507"

COLORS <- c("Dependent" = "#E57373", "Independent" = "#64B5F6")

# --- Functions ---

get_gene_lists <- function(rds_path, strip_version) {
  if(!file.exists(rds_path)) return(NULL)
  events <- readRDS(rds_path)
  
  # Group RI events
  ri <- events[events$EventType == "RI",]
  genes_dep   <- unique(ri$GeneID[ri$Group == "UFM1_dependent"])
  genes_indep <- unique(ri$GeneID[ri$Group == "UFM1_independent"])
  
  if(strip_version) {
    genes_dep   <- sub("\\..*", "", genes_dep)
    genes_indep <- sub("\\..*", "", genes_indep)
  }
  
  return(list(dep = genes_dep, ind = genes_indep))
}

# --- Execution ---

message(">>> Processing species: ", toupper(SPECIES))

# 1. Load Gene Lists
lists <- get_gene_lists(RDS_PATH, meta$strip_version)
if(is.null(lists)) stop("RDS file missing or empty for ", SPECIES)

# 2. BioMart Connection
message("Connecting to BioMart...")
mart <- useMart(meta$mart_name, dataset = meta$mart_data, host = meta$mart_host)

# 3. Fetch GO for our genes
all_event_genes <- unique(c(lists$dep, lists$ind))
go_data <- getBM(attributes = c(meta$id_type, meta$go_attr), 
                 filters = meta$id_type, 
                 values = all_event_genes, 
                 mart = mart)
colnames(go_data) <- c(meta$id_type, "go") # Rename for consistent access

# 4. Enrichment per Location
results <- list()
for(loc in names(LOC_TERMS)) {
  terms <- LOC_TERMS[[loc]]
  labeled_genes <- go_data[[meta$id_type]][go_data$go %in% terms]
  
  # Background Stats
  bg_hits_data <- getBM(attributes = c(meta$id_type), filters = 'go', values = terms, mart = mart)
  n_bg_hits <- length(unique(bg_hits_data[[1]]))
  n_bg_tot <- meta$bg_tot
  
  calc_stats <- function(target, label) {
    if(length(target) == 0) return(NULL)
    hits <- sum(target %in% labeled_genes)
    total <- length(target)
    
    # Fisher Test against Whole Genome Background
    ft <- fisher.test(matrix(c(hits, total - hits, n_bg_hits, n_bg_tot - n_bg_hits), nrow = 2))
    
    data.frame(
      Species = toupper(SPECIES), Category = loc, Group = label,
      OR = ft$estimate, Lower = ft$conf.int[1], Upper = ft$conf.int[2],
      P = ft$p.value, Count = hits, Total = total
    )
  }
  results[[paste0(loc, "_dep")]] <- calc_stats(lists$dep, "Dependent")
  results[[paste0(loc, "_ind")]] <- calc_stats(lists$ind, "Independent")
}

final_df <- do.call(rbind, results)
final_df$Sig <- ifelse(final_df$P < 0.05, "*", "")
final_df$LogP <- -log10(final_df$P)

# --- Output ---

# Save Raw Data
write.table(final_df, file.path(OUTDIR, "subcellular_localization_results.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Plotting
p <- ggplot(final_df, aes(x = Category, y = OR, color = Group)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = Lower, ymax = Upper), 
                  position = position_dodge(width = 0.6), size = 0.7) +
  geom_text(aes(y = Upper, label = Sig), 
            position = position_dodge(width = 0.6), vjust = -0.5, size = 6, show.legend = FALSE) +
  coord_flip() +
  scale_y_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10)) +
  scale_color_manual(values = COLORS) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  labs(title = paste(toupper(SPECIES), "Localization Enrichment"),
       subtitle = paste0("n: Dep=", length(lists$dep), ", Indep=", length(lists$ind)),
       y = "Odds Ratio (OR)", x = "")

pdf(file.path(OUTDIR, "subcellular_localization_plot.pdf"), width = 7, height = 5)
print(p)
dev.off()

png(file.path(OUTDIR, "subcellular_localization_plot.png"), width = 800, height = 600, res = 120)
print(p)
dev.off()

message("Done! Results saved to ", OUTDIR)
