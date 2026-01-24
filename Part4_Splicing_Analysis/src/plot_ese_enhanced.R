#!/usr/bin/env Rscript
#' plot_ese_enhanced.R
#' Enhanced ESE Visualization with all analysis results

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_ese_enhanced.R <results_dir> <output_dir>")
}

RESULTS_DIR <- args[1]
OUTPUT_DIR <- args[2]
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# === 1. CAG STATS BAR PLOT ===
cag_file <- file.path(RESULTS_DIR, "consensus_scan/CAG_Enhanced_Stats.csv")
if (file.exists(cag_file)) {
  cag <- read_csv(cag_file, show_col_types = FALSE)
  
  cag_long <- cag %>%
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
  
  p_cag <- ggplot(cag_long, aes(x = Group, y = Percentage, fill = Metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
              position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("Any Exonic CAG" = "#2171B5", "Terminal CAG (-3)" = "#6BAED6")) +
    labs(title = "CAG Trinucleotide Frequency in 5'SS Exons (±50bp)",
         x = "Group", y = "Percentage (%)", fill = "CAG Position") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, max(cag_long$Percentage, na.rm = TRUE) * 1.25)
  
  ggsave(file.path(OUTPUT_DIR, "CAG_Stats_Barplot.pdf"), p_cag, width = 8, height = 6)
  message("[INFO] CAG bar plot saved.")
}

# === 2. ESEFINDER HEATMAP ===
esefinder_file <- file.path(RESULTS_DIR, "consensus_scan/ESEfinder_Stats.csv")
if (file.exists(esefinder_file)) {
  esefinder <- read_csv(esefinder_file, show_col_types = FALSE)
  
  ese_summary <- esefinder %>%
    group_by(SR_Protein, Group) %>%
    summarize(Pct = mean(Pct_With_Motif, na.rm = TRUE), .groups = "drop") %>%
    mutate(Group = factor(Group, levels = c("UFM1_dependent", "UFM1_independent", "Constitutive")))
  
  p_ese <- ggplot(ese_summary, aes(x = Group, y = SR_Protein, fill = Pct)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = paste0(round(Pct, 1), "%")), size = 4) +
    scale_fill_gradient2(low = "white", mid = "#9ECAE1", high = "#08519C", midpoint = 10) +
    labs(title = "ESEfinder Motif Frequency by Group", x = "", y = "SR Protein", fill = "% Sequences") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(OUTPUT_DIR, "ESEfinder_Heatmap.pdf"), p_ese, width = 8, height = 5)
  message("[INFO] ESEfinder heatmap saved.")
}

# === 3. RBP ENRICHMENT SCATTER (TOP 10 HIGHLIGHTED) ===
rbp_file <- file.path(RESULTS_DIR, "RBP_Enrichment_Human.csv")
if (file.exists(rbp_file)) {
  rbp <- read_csv(rbp_file, show_col_types = FALSE)
  
  # Top 10 by absolute enrichment
  rbp <- rbp %>%
    mutate(
      Diff = abs(Log2_Dep - Log2_Ind),
      Rank = row_number(desc(Log2_Dep)),
      Label = ifelse(Rank <= 10, Motif_Name, NA_character_)
    )
  
  p_rbp <- ggplot(rbp, aes(x = Log2_Ind, y = Log2_Dep)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray70") +
    geom_point(aes(size = Hits_Dep, color = Rank <= 10), alpha = 0.6) +
    geom_text_repel(aes(label = Label), size = 3.5, max.overlaps = 15, 
                    box.padding = 0.5, point.padding = 0.3, fontface = "bold") +
    scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "gray60"), guide = "none") +
    scale_size_continuous(range = c(1, 6), name = "Hits") +
    labs(
      title = "RBP Motif Enrichment: UFM1-Dependent vs Independent",
      subtitle = "Top 10 enriched in Dependent highlighted (relative to Constitutive)",
      x = "Log2 Enrichment (Independent / Control)",
      y = "Log2 Enrichment (Dependent / Control)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
  
  ggsave(file.path(OUTPUT_DIR, "RBP_Enrichment_Scatter_Top10.pdf"), p_rbp, width = 10, height = 8)
  message("[INFO] RBP scatter plot saved.")
}

# === 4. RESCUE-ESE BAR ===
rescue_file <- file.path(RESULTS_DIR, "consensus_scan/RESCUE_ESE_Stats.csv")
if (file.exists(rescue_file)) {
  rescue <- read_csv(rescue_file, show_col_types = FALSE)
  
  p_rescue <- ggplot(rescue, aes(x = factor(Group, levels = c("UFM1_dependent", "UFM1_independent", "Constitutive")), 
                                  y = Pct_With_ESE, fill = Group)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = paste0(round(Pct_With_ESE, 1), "%")), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("UFM1_dependent" = "#E41A1C", "UFM1_independent" = "#377EB8", "Constitutive" = "gray60")) +
    labs(title = "RESCUE-ESE Hexamer Presence", x = "", y = "% Sequences with ESE") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    ylim(0, max(rescue$Pct_With_ESE, na.rm = TRUE) * 1.2)
  
  ggsave(file.path(OUTPUT_DIR, "RESCUE_ESE_Barplot.pdf"), p_rescue, width = 6, height = 5)
  message("[INFO] RESCUE-ESE bar plot saved.")
}

# === 5. BRANCH POINT BAR ===
bp_file <- file.path(RESULTS_DIR, "consensus_scan/BranchPoint_Stats.csv")
if (file.exists(bp_file)) {
  bp <- read_csv(bp_file, show_col_types = FALSE)
  
  p_bp <- ggplot(bp, aes(x = factor(Group, levels = c("UFM1_dependent", "UFM1_independent", "Constitutive")), 
                          y = Pct_With_BP, fill = Group)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = paste0(round(Pct_With_BP, 1), "%")), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("UFM1_dependent" = "#E41A1C", "UFM1_independent" = "#377EB8", "Constitutive" = "gray60")) +
    labs(title = "Branch Point Consensus (YNYURAY)", x = "", y = "% Sequences with BP") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    ylim(0, max(bp$Pct_With_BP, na.rm = TRUE) * 1.3)
  
  ggsave(file.path(OUTPUT_DIR, "BranchPoint_Barplot.pdf"), p_bp, width = 6, height = 5)
  message("[INFO] Branch Point bar plot saved.")
}

message("[INFO] All plots generated in: ", OUTPUT_DIR)
