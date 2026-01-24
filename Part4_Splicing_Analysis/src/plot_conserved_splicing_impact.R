#!/usr/bin/env Rscript

# ==============================================================================
# CONSERVED SPLICING IMPACT VISUALIZATION (Step 15) - Enhanced
# ==============================================================================
# Extracts mean PSI for conserved driver genes (Dependent & Independent)
# across Human and Mouse. Adds gene descriptions, IDs, and rounds PSI values.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(VennDiagram)
  library(grid)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript plot_conserved_splicing_impact.R --h_dep <path> --h_ind <path> --m_dep <path> --m_ind <path> --outdir <path>")
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

H_DEP <- params$h_dep
H_IND <- params$h_ind
M_DEP <- params$m_dep
M_IND <- params$m_ind
OUTDIR <- params$outdir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# --- Helpers ---
load_and_label <- function(path, species_name, category) {
  if(!file.exists(path)) return(NULL)
  df <- read.table(path, sep="\t", header=TRUE, comment.char = "", quote = "")
  if(nrow(df) == 0) return(NULL)
  df$Species <- species_name
  df$Category <- category
  return(df)
}

# --- 1. Load All Data ---
message("[INFO] Loading datasets...")
h_dep_df <- load_and_label(H_DEP, "Human", "Dependent")
h_ind_df <- load_and_label(H_IND, "Human", "Independent")
m_dep_df <- load_and_label(M_DEP, "Mouse", "Dependent")
m_ind_df <- load_and_label(M_IND, "Mouse", "Independent")

# --- 2. Venn Diagrams ---
message("[INFO] Generating Venn Diagrams...")
generate_venn <- function(h_genes, m_genes, category, filename, color1="#42A5F5", color2="#EF5350") {
  h_set <- unique(toupper(h_genes))
  m_set <- unique(toupper(m_genes))
  common <- intersect(h_set, m_set)
  
  pdf(file.path(OUTDIR, filename), width=6, height=6)
  grid.newpage()
  v <- draw.pairwise.venn(
    area1 = length(h_set),
    area2 = length(m_set),
    cross.area = length(common),
    category = c("Human", "Mouse"),
    fill = c(color1, color2),
    alpha = 0.5,
    main = paste("Conserved", category, "Genes"),
    main.cex = 1.2
  )
  grid.draw(v)
  dev.off()
  return(common)
}

dep_cons <- if(!is.null(h_dep_df) && !is.null(m_dep_df)) {
  generate_venn(h_dep_df$geneSymbol, m_dep_df$geneSymbol, "Dependent", "Venn_Conserved_Dependent.pdf")
} else character(0)

ind_cons <- if(!is.null(h_ind_df) && !is.null(m_ind_df)) {
  generate_venn(h_ind_df$geneSymbol, m_ind_df$geneSymbol, "Independent", "Venn_Conserved_Independent.pdf", color1="#81C784", color2="#FFB74D")
} else character(0)

# --- 3. Gene Descriptions ---
message("[INFO] Mapping gene descriptions...")
get_descriptions <- function(symbols, species_name) {
  db <- if(species_name == "Human") org.Hs.eg.db else org.Mm.eg.db
  syms_to_query <- if(species_name == "Human") toupper(symbols) else str_to_title(tolower(symbols))
  
  tryCatch({
    mapping <- AnnotationDbi::select(db, keys = syms_to_query, columns = c("GENENAME"), keytype = "SYMBOL")
    # Consolidate duplicates
    mapping <- mapping %>% 
      dplyr::group_by(SYMBOL) %>% 
      dplyr::summarize(Description = paste(unique(GENENAME), collapse="; "))
    # Return as named vector for easy lookup
    res <- setNames(mapping$Description, toupper(mapping$SYMBOL))
    return(res)
  }, error = function(e) {
    message("[WARN] Failed to map descriptions for ", species_name)
    u_syms <- unique(toupper(symbols))
    return(setNames(rep("NA", length(u_syms)), u_syms))
  })
}

all_symbols <- unique(toupper(c(h_dep_df$geneSymbol, h_ind_df$geneSymbol, m_dep_df$geneSymbol, m_ind_df$geneSymbol)))
h_desc <- get_descriptions(all_symbols, "Human")
m_desc <- get_descriptions(all_symbols, "Mouse")

# --- 4. Consolidate and Format Plot Data ---
message("[INFO] Consolidating data and rounding PSI values...")
all_data <- bind_rows(h_dep_df, h_ind_df, m_dep_df, m_ind_df)
relevant_data <- all_data %>%
  mutate(geneSymbol_Upper = toupper(geneSymbol)) %>%
  filter(EventType == "RI", 
         (Category == "Dependent" & geneSymbol_Upper %in% dep_cons) |
         (Category == "Independent" & geneSymbol_Upper %in% ind_cons)) %>%
  # Rounding PSI to 3 decimal places
  mutate(across(c(mean_psi_ctrl.WT, mean_psi_case.WT, mean_psi_ctrl.UFM, mean_psi_case.UFM), ~round(as.numeric(.), 3))) %>%
  dplyr::select(GeneID, geneSymbol, geneSymbol_Upper, Category, Species,
                WT_Ctrl = mean_psi_ctrl.WT, WT_Case = mean_psi_case.WT,
                UFM_Ctrl = mean_psi_ctrl.UFM, UFM_Case = mean_psi_case.UFM)

# Add Descriptions
relevant_data$Description <- ifelse(relevant_data$Species == "Human", 
                                   h_desc[relevant_data$geneSymbol_Upper], 
                                   m_desc[relevant_data$geneSymbol_Upper])

if (nrow(relevant_data) == 0) {
  stop("No matching RI events for conserved genes.")
}

# Export Data
write.table(relevant_data, file.path(OUTDIR, "Conserved_Mammalian_Splicing_Summary.tsv"), 
            sep="\t", row.names=FALSE, quote=FALSE)
message("[INFO] Published Consolidated TSV with GeneID, Descriptions, and 3-decimal PSI.")

# --- 5. Splicing Impact Plots ---
plot_df <- relevant_data %>%
  pivot_longer(cols = c(WT_Ctrl, WT_Case, UFM_Ctrl, UFM_Case), 
               names_to = "Condition", values_to = "PSI") %>%
  mutate(Genotype = factor(ifelse(grepl("WT", Condition), "WT", "UFM1"), levels=c("WT", "UFM1")),
         Treatment = factor(ifelse(grepl("Ctrl", Condition), "Ctrl", "Case"), levels=c("Ctrl", "Case")))

# Dependent Plot
dep_plot_df <- plot_df %>% filter(Category == "Dependent")
p_dep <- ggplot(dep_plot_df, aes(x = Treatment, y = PSI, color = Genotype, group = interaction(geneSymbol, Genotype, Species))) +
  geom_point(size = 2, alpha = 0.8) +
  geom_line(alpha = 0.5) +
  facet_grid(Species ~ geneSymbol, scales = "free_y") +
  theme_bw() +
  scale_color_manual(values = c("WT" = "#42A5F5", "UFM1" = "#EF5350")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text.x = element_text(size = 7, angle = 90),
        legend.position = "top") +
  labs(title = "Conserved UFM1-Dependent Splicing Targets", x = "", y = "PSI (Rounded)")

w_dep <- max(12, length(unique(dep_plot_df$geneSymbol)) * 0.7)
ggsave(file.path(OUTDIR, "Conserved_Dependent_PSI_Facets.png"), plot = p_dep, width = w_dep, height = 7, dpi = 200)

message("Step 15: Cross-species enhanced visualization complete.")
