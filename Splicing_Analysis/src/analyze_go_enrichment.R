#!/usr/bin/env Rscript

# ==============================================================================
# GO ENRICHMENT ANALYSIS (clusterProfiler) - Comparative Mode
# ==============================================================================
# Performs comparative enrichment between Dependent and Independent groups.
# Follows user-provided snippet for bitr + compareCluster.
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(GenomicRanges)
  library(tidyverse)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript analyze_go_enrichment.R --species <sp> --rds <path> --outdir <path>")
}

# Parse Args
parse_args <- function(args) {
  params <- list()
  for (i in seq(1, length(args), by=2)) {
    key <- gsub("--", "", args[i])
    params[[key]] <- args[i+1]
  }
  return(params)
}
params <- parse_args(args)

SPECIES <- tolower(params$species)
RDS_PATH <- params$rds
OUTDIR <- params$outdir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Species Meta for clusterProfiler
if (SPECIES == "human") {
  library(org.Hs.eg.db)
  ORG_DB <- "org.Hs.eg.db"
  KEY_TYPE <- "ENSEMBL"
} else if (SPECIES == "mouse") {
  library(org.Mm.eg.db)
  ORG_DB <- "org.Mm.eg.db"
  KEY_TYPE <- "ENSEMBL"
} else if (SPECIES == "arabidopsis") {
  library(org.At.tair.db)
  ORG_DB <- "org.At.tair.db"
  KEY_TYPE <- "TAIR"
} else {
  stop("Unsupported species: ", SPECIES)
}

# --- Functions ---

get_gene_lists <- function(rds_path, species) {
  if(!file.exists(rds_path)) return(NULL)
  events <- readRDS(rds_path)
  df_events <- if (inherits(events, "GRanges")) as.data.frame(events) else events
  # ri <- df_events[df_events$EventType == "RI", ] # Deleted hardcoded RI
  genes_dep   <- unique(df_events$GeneID[df_events$Group == "UFM1_dependent"])
  genes_indep <- unique(df_events$GeneID[df_events$Group == "UFM1_independent"])
  
  # Clean IDs
  if(species != "arabidopsis") {
    genes_dep <- sub("\\..*", "", genes_dep); genes_indep <- sub("\\..*", "", genes_indep)
  }
  return(list(Dependent = genes_dep, Independent = genes_indep))
}

# --- Execution ---

message(">>> Step 14: Comparative GO Enrichment for ", toupper(SPECIES))
lists <- get_gene_lists(RDS_PATH, SPECIES)

tryCatch({
  message("Mapping IDs to Entrez...")
  entrez_dep <- bitr(lists$Dependent, fromType=KEY_TYPE, toType="ENTREZID", OrgDb=ORG_DB)$ENTREZID
  entrez_ind <- bitr(lists$Independent, fromType=KEY_TYPE, toType="ENTREZID", OrgDb=ORG_DB)$ENTREZID
  
  message("Running compareCluster (BP)...")
  ck <- compareCluster(geneCluster = list(Dependent=entrez_dep, Independent=entrez_ind),
                       fun = "enrichGO", OrgDb = ORG_DB, ont = "BP", 
                       pvalueCutoff = 0.5, qvalueCutoff = 0.5) # Discovery Mode
  
  if (is.null(ck) || nrow(as.data.frame(ck)) == 0) {
    message("No significant pathways found in comparison.")
  } else {
    # 1. Save Full Results
    write.table(as.data.frame(ck), file.path(OUTDIR, "Comparative_GO_Results.tsv"), 
                sep="\t", row.names=FALSE, quote=FALSE)
    
    # 2. Plot
    p_enrich <- dotplot(ck, showCategory=10, title=paste(toupper(SPECIES), "Pathways")) + 
      theme(axis.text.x = element_text(angle=45, hjust=1))
    
    ggsave(file.path(OUTDIR, "Comparative_GO_DotPlot.png"), plot = p_enrich, width = 10, height = 10, dpi = 150)
    
    # 3. Save Drivers (as requested)
    ego_df <- as.data.frame(ck) %>% filter(Cluster == "Dependent")
    write.csv(ego_df, file.path(OUTDIR, paste0("UFM1_", toupper(SPECIES), "_Pathway_Drivers.csv")), row.names=FALSE)
  }
  
}, error = function(e) {
  message("Enrichment failed: ", e$message)
})

message("Done! Results in ", OUTDIR)
