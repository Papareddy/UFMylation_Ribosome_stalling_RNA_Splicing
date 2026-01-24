#!/usr/bin/env Rscript
#' plot_ese_enrichment.R
#' Step C: ESE Motif Visualization
#' 
#' Generates:
#' 1. Comparative scatter plot (Log2 enrichment)
#' 2. CAG stats visualization

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
})

# === ARGUMENTS ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_ese_enrichment.R <motif_enrichment_csv> <cag_stats_csv> <output_dir>")
}

enrichment_file <- args[1]
cag_stats_file <- args[2]
output_dir <- ifelse(length(args) >= 3, args[3], dirname(enrichment_file))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define RBP families to highlight
srsf_family <- c("SRSF1", "SRSF2", "SRSF3", "SRSF4", "SRSF5", "SRSF6", "SRSF7", 
                 "SRSF9", "SRSF10", "SRp20", "SRp40", "SRp55", "SC35", "ASF", "SF2")
pcbp_family <- c("PCBP1", "PCBP2", "PCBP3", "PCBP4", "hnRNPK", "hnRNPE1", "hnRNPE2")
rbm39_family <- c("RBM39", "CAPER", "CAPERalpha")

highlight_rbps <- c(srsf_family, pcbp_family, rbm39_family)

# === PLOT 1: SCATTER PLOT ===
message("[INFO] Generating scatter plot...")

if (file.exists(enrichment_file)) {
  enrichment <- read_csv(enrichment_file, show_col_types = FALSE)
  
  if (nrow(enrichment) > 0) {
    # Classify motifs
    enrichment <- enrichment %>%
      mutate(
        Family = case_when(
          str_detect(Motif, regex(paste(srsf_family, collapse = "|"), ignore_case = TRUE)) ~ "SRSF",
          str_detect(Motif, regex(paste(pcbp_family, collapse = "|"), ignore_case = TRUE)) ~ "PCBP",
          str_detect(Motif, regex(paste(rbm39_family, collapse = "|"), ignore_case = TRUE)) ~ "RBM39",
          TRUE ~ "Other"
        ),
        Label = ifelse(Family != "Other", Motif, NA_character_)
      )
    
    # Create scatter plot
    p_scatter <- ggplot(enrichment, aes(x = Log2_Enrichment_Indep_vs_Control, 
                                         y = Log2_Enrichment_Dep_vs_Control)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray70") +
      geom_point(aes(color = Family, size = Family), alpha = 0.7) +
      geom_text_repel(aes(label = Label), size = 3, max.overlaps = 20, 
                      box.padding = 0.5, point.padding = 0.3) +
      scale_color_manual(values = c("SRSF" = "#E41A1C", "PCBP" = "#377EB8", 
                                     "RBM39" = "#4DAF4A", "Other" = "gray60")) +
      scale_size_manual(values = c("SRSF" = 3.5, "PCBP" = 3.5, "RBM39" = 3.5, "Other" = 1.5)) +
      labs(
        title = "RBP Motif Enrichment: UFM1-Dependent vs Independent",
        subtitle = "Relative to Constitutive Introns (Control)",
        x = "Log2 Enrichment (Independent / Control)",
        y = "Log2 Enrichment (Dependent / Control)",
        color = "RBP Family"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      ) +
      guides(size = "none")
    
    ggsave(file.path(output_dir, "RBP_Enrichment_Scatter.pdf"), p_scatter, 
           width = 10, height = 8)
    message("[INFO] Scatter plot saved.")
  }
} else {
  message("[WARN] Enrichment file not found, skipping scatter plot.")
}

# === PLOT 2: CAG STATS VISUALIZATION ===
message("[INFO] Generating CAG stats visualization...")

if (file.exists(cag_stats_file)) {
  cag_stats <- read_csv(cag_stats_file, show_col_types = FALSE)
  
  if (nrow(cag_stats) > 0) {
    # Prepare data for plotting
    cag_long <- cag_stats %>%
      dplyr::select(Group, Pct_Exonic_CAG, Pct_Terminal_CAG) %>%
      pivot_longer(cols = c(Pct_Exonic_CAG, Pct_Terminal_CAG), 
                   names_to = "Metric", values_to = "Percentage") %>%
      mutate(
        Metric = case_when(
          Metric == "Pct_Exonic_CAG" ~ "Any Exonic CAG",
          Metric == "Pct_Terminal_CAG" ~ "Terminal CAG (-3)",
          TRUE ~ Metric
        ),
        Group = factor(Group, levels = c("UFM1_dependent", "UFM1_independent", "Constitutive"))
      )
    
    # Bar plot
    p_cag <- ggplot(cag_long, aes(x = Group, y = Percentage, fill = Metric)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
                position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
      scale_fill_manual(values = c("Any Exonic CAG" = "#2171B5", "Terminal CAG (-3)" = "#6BAED6")) +
      labs(
        title = "CAG Trinucleotide Frequency in 5' Splice Site Exons",
        subtitle = "Percentage of introns with CAG in upstream exon (-30 to -1 bp)",
        x = "Group",
        y = "Percentage of Sequences (%)",
        fill = "CAG Position"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "top",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      ylim(0, max(cag_long$Percentage, na.rm = TRUE) * 1.2)
    
    ggsave(file.path(output_dir, "CAG_Stats_Barplot.pdf"), p_cag, 
           width = 8, height = 6)
    message("[INFO] CAG bar plot saved.")
    
    # Print summary table
    message("\n", paste(rep("=", 60), collapse = ""))
    message("CAG ANALYSIS SUMMARY")
    message(paste(rep("=", 60), collapse = ""))
    print(cag_stats)
    message(paste(rep("=", 60), collapse = ""), "\n")
  }
} else {
  message("[WARN] CAG stats file not found, skipping CAG visualization.")
}

message("[INFO] Step C complete. Visualizations generated.")
