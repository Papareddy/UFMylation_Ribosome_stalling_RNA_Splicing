#!/usr/bin/env Rscript

library(gprofiler2)
library(dplyr)
library(readr)
library(ggplot2)

# Arguments
args <- commandArgs(trailingOnly = TRUE)
dep_file <- args[1] # hs_RI_0.25/human/nucleus/step01_data_prep/UFM1_dependent.tsv
indep_file <- args[2] # hs_RI_0.25/human/nucleus/step01_data_prep/UFM1_independent.tsv
outdir <- args[3] # hs_RI_0.25/human/nucleus/functional_analysis

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load genes
dep_df <- read_tsv(dep_file, show_col_types = FALSE)
indep_df <- read_tsv(indep_file, show_col_types = FALSE)

dep_genes <- unique(na.omit(dep_df$geneSymbol))
indep_genes <- unique(na.omit(indep_df$geneSymbol))

cat("[INFO] Dependent genes (n):", length(dep_genes), "\n")
cat("[INFO] Independent genes (n):", length(indep_genes), "\n")

# Save gene lists
writeLines(dep_genes, file.path(outdir, "dependent_genes.txt"))
writeLines(indep_genes, file.path(outdir, "independent_genes.txt"))

# Enrichment
run_enr <- function(genes, name) {
  cat("[INFO] Running enrichment for:", name, "...\n")
  res <- gost(query = genes, organism = "hsapiens", ordered_query = FALSE, 
              multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
              measure_underrepresentation = FALSE, evcodes = TRUE, 
              user_threshold = 0.05, correction_method = "gSCS", 
              domain_scope = "annotated", custom_bg = NULL, 
              numeric_ns = "", sources = c("GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC"))
  
  if (!is.null(res$result)) {
    write_csv(res$result, file.path(outdir, paste0(name, "_enrichment_results.csv")))
    
    # Save a simple plot if possible
    p <- gostplot(res, capped = TRUE, interactive = FALSE)
    ggsave(file.path(outdir, paste0(name, "_enrichment_plot.pdf")), p, width = 10, height = 7)
    
    # Extract top terms
    top_terms <- res$result %>% 
      filter(precision > 0.1) %>%
      arrange(p_value) %>% 
      head(20) %>%
      select(source, term_id, term_name, p_value, query_size, intersection_size)
    
    write_csv(top_terms, file.path(outdir, paste0(name, "_top_terms.csv")))
    return(top_terms)
  } else {
    cat("[WARN] No significant enrichment found for:", name, "\n")
    return(NULL)
  }
}

dep_top <- run_enr(dep_genes, "UFM1_dependent")
indep_top <- run_enr(indep_genes, "UFM1_independent")

cat("[INFO] Functional analysis completed in:", outdir, "\n")
