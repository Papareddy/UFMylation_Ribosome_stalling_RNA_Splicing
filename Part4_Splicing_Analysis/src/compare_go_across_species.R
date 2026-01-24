#!/usr/bin/env Rscript

# ==============================================================================
# CROSS-SPECIES GO ENRICHMENT COMPARISON (Pathway Drivers) - Robust Version
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(VennDiagram)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript compare_go_across_species.R --human_drivers <path> --mouse_drivers <path> --outdir <path>")
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

HUMAN_PATH <- params$human_drivers
MOUSE_PATH <- params$mouse_drivers
OUTDIR <- params$outdir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# --- Load Data ---
load_drivers <- function(path, species_name) {
  if(!file.exists(path)) return(NULL)
  df <- read.csv(path, stringsAsFactors = FALSE)
  if(nrow(df) > 0) {
    df$Species <- species_name
    return(df)
  } else {
    message("[INFO] No drivers found for ", species_name)
    return(data.frame())
  }
}

h_df <- load_drivers(HUMAN_PATH, "Human")
m_df <- load_drivers(MOUSE_PATH, "Mouse")

if(is.null(h_df) || is.null(m_df)) stop("Critical driver files missing.")

# Check if both have data
if (nrow(h_df) > 0 && nrow(m_df) > 0) {
  common_terms <- intersect(h_df$ID, m_df$ID)
  message("Found ", length(common_terms), " overlapping pathways.")

  merged_common <- h_df %>% 
    filter(ID %in% common_terms) %>%
    select(ID, Description, p.adjust_h = p.adjust, GeneRatio_h = GeneRatio) %>%
    inner_join(m_df %>% filter(ID %in% common_terms) %>% select(ID, p.adjust_m = p.adjust, GeneRatio_m = GeneRatio), by="ID") %>%
    mutate(Conserved_Score = ((-log10(p.adjust_h)) + (-log10(p.adjust_m))) / 2)

  write.table(merged_common[order(-merged_common$Conserved_Score), ], 
              file.path(OUTDIR, "Conserved_Dependent_Pathways.tsv"), 
              sep="\t", row.names=FALSE, quote=FALSE)

  # Venn
  pdf(file.path(OUTDIR, "Venn_Overlap_Dependent.pdf"), width=6, height=6)
  grid.newpage()
  draw.pairwise.venn(area1 = nrow(h_df), area2 = nrow(m_df), cross.area = length(common_terms),
                     category = c("Human Dep", "Mouse Dep"), fill = c("#E57373", "#64B5F6"), alpha = 0.5)
  dev.off()

  # Dot Plot
  if(nrow(merged_common) > 0) {
    top_common <- merged_common %>% arrange(-Conserved_Score) %>% head(25)
    plot_df <- top_common %>% 
      pivot_longer(cols = c(p.adjust_h, p.adjust_m), names_to = "Species", values_to = "p_adj") %>%
      mutate(Species = ifelse(Species == "p.adjust_h", "Human", "Mouse"))
    
    p <- ggplot(plot_df, aes(x = reorder(Description, Conserved_Score), y = -log10(p_adj), color = Species)) +
      geom_point(size = 5, alpha = 0.8) + coord_flip() + theme_bw() +
      scale_color_manual(values = c("Human" = "#E57373", "Mouse" = "#64B5F6")) +
      labs(title = "Conserved UFM1-Dependent Pathways", x = "", y = "-log10(adjusted p-value)")
    ggsave(file.path(OUTDIR, "Conserved_Dependent_DotPlot.png"), plot = p, width = 11, height = 8, dpi = 150)
  }
} else {
  message("[WARN] Cross-species overlap not possible: one or both gene sets yielded no significant pathways.")
}

message("Cross-species comparison task finished.")
