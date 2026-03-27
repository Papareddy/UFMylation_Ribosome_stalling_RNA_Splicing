#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(optparse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.At.tair.db)

# --- Argument Parsing ---
option_list <- list(
  make_option(c("-o", "--outdir"), type="character", default="results/cross_species_visualization", 
              help="Output directory [default= %default]", metavar="character"),
  make_option(c("-r", "--results_dir"), type="character", default="results", 
              help="Base results directory to look for species folders [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

out_dir <- opt$outdir
res_base <- opt$results_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Expected subpaths for GO results in the pipeline structure:
# {res_base}/{species}/nucleus/step14_go_enrichment/Comparative_GO_Results.tsv
# Note: Mouse might be in 'total' instead of 'nucleus' depending on pipeline run
species_list <- c("human", "mouse", "arabidopsis")

# Helper to load and tag
load_res <- function(sp) {
  # Try nucleus first, then total (common variation in this pipeline)
  paths <- c(
    file.path(res_base, sp, "nucleus", "step14_go_enrichment", "Comparative_GO_Results.tsv"),
    file.path(res_base, sp, "total", "step14_go_enrichment", "Comparative_GO_Results.tsv"),
    file.path(res_base, sp, "Comparative_GO_Results.tsv") # Fallback
  )
  
  for (path in paths) {
    if (file.exists(path)) {
      cat(sprintf("[INFO] Loading %s data from %s\n", sp, path))
      df <- read_tsv(path, show_col_types = FALSE) %>% mutate(Species = sp)
      return(df)
    }
  }
  cat(sprintf("[WARN] No GO results found for %s in expected locations.\n", sp))
  return(NULL)
}

# Load all
data_list <- lapply(species_list, load_res)
data_list <- data_list[!sapply(data_list, is.null)]

if (length(data_list) == 0) {
    stop("[ERROR] No GO results found for any species. Cannot proceed.")
}

all_data <- bind_rows(data_list)

# --- Gene ID to Symbol Converters ---
get_symbols <- function(ids, species) {
    if (length(ids) == 0) return("")
    id_list <- unique(unlist(strsplit(ids, "/")))
    
    symbols <- tryCatch({
        if (species == "human") {
            mapIds(org.Hs.eg.db, keys = id_list, column = "SYMBOL", keytype = "ENTREZID")
        } else if (species == "mouse") {
            mapIds(org.Mm.eg.db, keys = id_list, column = "SYMBOL", keytype = "ENTREZID")
        } else if (species == "arabidopsis") {
            mapIds(org.At.tair.db, keys = id_list, column = "SYMBOL", keytype = "TAIR")
        }
    }, error = function(e) return(id_list))
    
    symbols <- symbols[!is.na(symbols)]
    return(paste(symbols, collapse = ", "))
}

# --- Process per Cluster ---
clusters <- c("Dependent", "Independent")
gene_by_term_dir <- file.path(out_dir, "Genes_By_Eukaryotic_Term")
dir.create(gene_by_term_dir, showWarnings = FALSE, recursive = TRUE)

for (clust in clusters) {
    cat(sprintf("\n[INFO] Processing Cluster: %s\n", clust))
    
    # 1. Get Top Human Terms for this cluster (pvalue < 0.01)
    top_human_pool <- all_data %>%
        dplyr::filter(Species == "human", Cluster == clust, pvalue < 0.01) %>%
        arrange(pvalue) %>%
        head(100)
    
    if (nrow(top_human_pool) == 0) {
        cat(sprintf("[WARN] No pvalue < 0.01 terms for %s in Human. Skipping.\n", clust))
        next
    }
    
    top_human_ids <- top_human_pool$ID
    
    # 2. Rule: Human (pvalue < 0.01), Mouse (pvalue < 0.01), Arabidopsis (Count > 1)
    plot_data <- all_data %>%
        dplyr::filter(Cluster == clust, ID %in% top_human_ids) %>%
        dplyr::filter(
            (Species == "human" & pvalue < 0.01) |
            (Species == "mouse" & pvalue < 0.01) |
            (Species == "arabidopsis" & Count > 1)
        ) %>%
        mutate(neg_log10p = -log10(pvalue))
    
    if (nrow(plot_data) == 0) next

    # Clean up factor levels
    active_ids <- unique(plot_data$ID)
    top_human_pool_subset <- top_human_pool %>% dplyr::filter(ID %in% active_ids)
    
    plot_data$Description <- factor(plot_data$Description, levels = rev(top_human_pool_subset$Description))
    plot_data$Species <- factor(plot_data$Species, levels = c("human", "mouse", "arabidopsis"))
    
    # 3. Create Eukaryotic Landscape Plot (Viridis Aesthetic)
    p <- ggplot(plot_data, aes(x = Species, y = Description, size = Count, color = neg_log10p)) +
        geom_point() +
        scale_color_viridis_c(option = "viridis", name = "-log10(p-value)") +
        scale_size_continuous(range = c(2, 10), name = "Gene Count") +
        theme_bw() +
        labs(title = sprintf("Eukaryotic Conserved Landscape: %s RI", clust),
             subtitle = "Human/Mouse (p < 0.01), Arabidopsis (Count > 1)",
             x = "", y = "GO Term") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face="bold"),
              axis.text.y = element_text(size = 8))
    
    ggsave(file.path(out_dir, sprintf("Eukaryotic_Landscape_%s.png", clust)), p, width = 12, height = 20)
    ggsave(file.path(out_dir, sprintf("Eukaryotic_Landscape_%s.pdf", clust)), p, width = 12, height = 20)
    
    # 4. Extract Gene Details
    cat("[INFO] Mapping gene symbols and exporting term-level files...")
    gene_details <- all_data %>%
        dplyr::filter(Cluster == clust, ID %in% active_ids) %>%
        dplyr::filter(
            (Species == "human" & pvalue < 0.01) |
            (Species == "mouse" & pvalue < 0.01) |
            (Species == "arabidopsis" & Count > 1)
        )
    
    gene_details$GeneSymbols <- mapply(function(ids, sp) get_symbols(ids, sp), 
                                     gene_details$geneID, gene_details$Species)
    
    write_tsv(gene_details, file.path(out_dir, sprintf("Eukaryotic_Genes_%s.tsv", clust)))
    
    # Save individual term files
    cluster_term_dir <- file.path(gene_by_term_dir, clust)
    dir.create(cluster_term_dir, showWarnings = FALSE)
    
    for (go_id in unique(gene_details$ID)) {
        term_data <- gene_details %>% dplyr::filter(ID == go_id)
        clean_id <- gsub(":", "_", go_id)
        write_tsv(term_data, file.path(cluster_term_dir, paste0(clean_id, ".tsv")))
    }
}

cat(sprintf("\n[DONE] Eukaryotic Landscape Analysis complete. Results in %s\n", out_dir))
