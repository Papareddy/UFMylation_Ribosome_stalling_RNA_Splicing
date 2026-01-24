#!/usr/bin/env Rscript

# ==============================================================================
# FUNCTIONAL CONSERVATION MAPPING (Step 16) - Enhanced with Venns
# ==============================================================================
# Intersects GO pathways between Human and Mouse.
# For shared pathways, calculates gene-level overlap.
# Generates per-theme Venn diagrams for high-priority categories.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(VennDiagram)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript map_functional_conservation.R --human_go <path> --mouse_go <path> --outdir <path>")
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

HUMAN_GO <- params$human_go
MOUSE_GO <- params$mouse_go
OUTDIR <- params$outdir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Create subfolder for Venns
VENN_DIR <- file.path(OUTDIR, "per_theme_venns")
dir.create(VENN_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load Data ---
message("[INFO] Loading GO results...")
h_go <- read.table(HUMAN_GO, sep="\t", header=TRUE, comment.char = "", quote = "")
m_go <- read.table(MOUSE_GO, sep="\t", header=TRUE, comment.char = "", quote = "")

# --- 2. Helper: Map Entrez to Symbols ---
map_ids <- function(id_string, species) {
  ids <- unlist(strsplit(id_string, "/"))
  db <- if(species == "human") org.Hs.eg.db else org.Mm.eg.db
  suppressMessages({
    res <- AnnotationDbi::select(db, keys = ids, columns = "SYMBOL", keytype = "ENTREZID")
  })
  return(unique(toupper(res$SYMBOL)))
}

# --- 3. Functional Breakdown Logic ---
# High-priority themes
PRIORITY_THEMES <- c("endoplasmic reticulum organization", "DNA repair", 
                    "glycerophospholipid metabolic process", "response to radiation",
                    "cellular response to DNA damage stimulus", "organelle organization")

# Helper to shorten theme names for display
shorten_theme <- function(x) {
  x <- gsub("endoplasmic reticulum", "ER", x)
  x <- gsub("glycerophospholipid metabolic process", "Lipid Met", x)
  x <- gsub("cellular response to DNA damage stimulus", "DNA Damage", x)
  x <- gsub("response to radiation", "Radiation", x)
  x <- gsub("organelle organization", "Organelle Org", x)
  return(x)
}

# Custom UpSet-style Plot Function for Aggregated Data
draw_aggregated_upset <- function(h_genes, m_genes, conserved_map, title, color_main, filename) {
  # h_genes / m_genes: vectors of unique gene symbols
  # conserved_map: named vector where names are genes and values are Themes
  
  common_genes <- intersect(h_genes, m_genes)
  only_h <- setdiff(h_genes, m_genes)
  only_m <- setdiff(m_genes, h_genes)
  
  # Format Labels for Conserved
  # e.g. "ESYT1 (ER)"
  labeled_genes <- c()
  for(g in common_genes) {
    if(g %in% names(conserved_map)) {
      lbl <- paste0(g, " (", conserved_map[g], ")")
      labeled_genes <- c(labeled_genes, lbl)
    } else {
      labeled_genes <- c(labeled_genes, g)
    }
  }
  
  # Limit label count if too many (though for this set usually small)
  label_text <- paste(labeled_genes, collapse="\n")
  
  # Plot Data
  plot_data <- data.frame(
    Category = c("Human Only", "Shared", "Mouse Only"),
    Count = c(length(only_h), length(common_genes), length(only_m)),
    Label = c("", label_text, "")
  )
  plot_data$Category <- factor(plot_data$Category, levels = c("Human Only", "Shared", "Mouse Only"))
  
  # 1. Bar Plot
  p_bar <- ggplot(plot_data, aes(x = Category, y = Count)) +
    geom_bar(stat = "identity", fill = color_main, alpha = 0.8, width = 0.6) +
    geom_text(aes(label = Count), vjust = -0.5, fontface = "bold") +
    geom_text(aes(label = Label), y = plot_data$Count[2] * 0.5, size = 3, color = "white", fontface="bold") +
    theme_minimal() +
    labs(title = title, subtitle = "Combined Priority Themes", y = "Gene Count", x = "") +
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
  
  # 2. Matrix Plot
  matrix_data <- data.frame(
    Category = rep(c("Human Only", "Shared", "Mouse Only"), each=2),
    Set = rep(c("Human", "Mouse"), 3),
    Present = c(1, 0,  1, 1,  0, 1)
  )
  matrix_data$Category <- factor(matrix_data$Category, levels = c("Human Only", "Shared", "Mouse Only"))
  
  p_matrix <- ggplot(matrix_data, aes(x = Category, y = Set)) +
    geom_point(aes(color = as.factor(Present)), size = 5) +
    scale_color_manual(values = c("0" = "transparent", "1" = "black")) +
    annotate("segment", x = 2, xend = 2, y = 1, yend = 2, size = 1.2, color = "black") +
    theme_minimal() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_text(face = "bold", size = 10)) +
    labs(y = "")
  
  # Combine
  library(patchwork)
  final_plot <- p_bar / p_matrix + plot_layout(heights = c(3, 1))
  ggsave(filename, plot = final_plot, width = 6, height = 7)
}

generate_aggregated_plot <- function(cluster_name) {
  # Filter for priority themes
  h_sub <- h_go %>% filter(Cluster == cluster_name, Description %in% PRIORITY_THEMES)
  m_sub <- m_go %>% filter(Cluster == cluster_name, Description %in% PRIORITY_THEMES)
  
  if(nrow(h_sub) == 0 || nrow(m_sub) == 0) return()
  
  # Collect all genes
  all_h_genes <- c()
  all_m_genes <- c()
  gene_theme_map <- c() # For annotation
  
  # Process Human
  for(i in 1:nrow(h_sub)) {
    term <- h_sub$Description[i]
    term_short <- shorten_theme(term)
    syms <- map_ids(h_sub$geneID[i], "human")
    all_h_genes <- c(all_h_genes, syms)
    
    # Map gene to theme (if not already mapped)
    # If gene in multiple themes, maybe list? For now overwrite or combine
    # Let's combine: "ESYT1 (ER, Org)"
    for(s in syms) {
      if(s %in% names(gene_theme_map)) {
         # append if new
         curr <- gene_theme_map[s]
         if(!grepl(term_short, curr)) {
           gene_theme_map[s] <- paste(curr, term_short, sep="/")
         }
      } else {
        gene_theme_map[s] <- term_short
      }
    }
  }
  
  # Process Mouse (just for the set)
  for(i in 1:nrow(m_sub)) {
    syms <- map_ids(m_sub$geneID[i], "mouse")
    all_m_genes <- c(all_m_genes, syms)
  }
  
  unique_h <- unique(all_h_genes)
  unique_m <- unique(all_m_genes)
  
  out_file <- file.path(OUTDIR, paste0("UpSet_Aggregated_", cluster_name, ".png"))
  col <- if(cluster_name == "Dependent") "#42A5F5" else "#81C784"
  
  draw_aggregated_upset(unique_h, unique_m, gene_theme_map, 
                        paste(cluster_name, " Conservation"), col, out_file)
}

# --- Execute Main Logic ---

# 1. Generate Full TSV (Existing logic preserved for data record)
process_overlap_tsv <- function(cluster_name) {
  h_sub <- h_go %>% filter(Cluster == cluster_name)
  m_sub <- m_go %>% filter(Cluster == cluster_name)
  common_terms <- intersect(h_sub$Description, m_sub$Description)
  if(length(common_terms) == 0) return(NULL)
  results <- list()
  for (term in common_terms) {
    h_ids <- h_sub %>% filter(Description == term) %>% pull(geneID) %>% head(1)
    m_ids <- m_sub %>% filter(Description == term) %>% pull(geneID) %>% head(1)
    common_genes <- intersect(map_ids(h_ids, "human"), map_ids(m_ids, "mouse"))
    results[[term]] <- data.frame(Category=cluster_name, GO_Term=term, Conserved_Count=length(common_genes), Conserved_Genes=paste(common_genes, collapse=", "))
  }
  return(bind_rows(results))
}

dep_res <- process_overlap_tsv("Dependent")
ind_res <- process_overlap_tsv("Independent")
final_res <- bind_rows(dep_res, ind_res)
if (nrow(final_res) > 0) {
  write.table(final_res, file.path(OUTDIR, "Functional_Conservation_Gene_Breakdown.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  
  # Summary Bar Plot
  high_res <- final_res %>% filter(tolower(GO_Term) %in% tolower(PRIORITY_THEMES))
  if(nrow(high_res) > 0) {
     p <- ggplot(high_res, aes(x = reorder(GO_Term, Conserved_Count), y = Conserved_Count, fill = Category)) +
      geom_bar(stat="identity", position="dodge") + coord_flip() + theme_bw() +
      scale_fill_manual(values = c("Dependent" = "#42A5F5", "Independent" = "#81C784"))
     ggsave(file.path(OUTDIR, "Conserved_Themes_Summary.png"), plot = p, width = 8, height = 6)
  }
}

# 2. Generate Aggregated UpSet Plots (New Request)
message("[INFO] Generating Aggregated UpSet Plot for Dependent...")
generate_aggregated_plot("Dependent")
message("[INFO] Generating Aggregated UpSet Plot for Independent...")
generate_aggregated_plot("Independent")

message("Step 16: Functional conservation analysis (Aggregated UpSets) complete.")
