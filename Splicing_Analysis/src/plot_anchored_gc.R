#!/usr/bin/env Rscript

# plot_anchored_gc.R
# Usage: Rscript plot_anchored_gc.R <input_tsv> <window_size> <output_pdf>

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript plot_anchored_gc.R <input_tsv> <window_size> <output_pdf>\nExample: Rscript plot_anchored_gc.R results/data.tsv 500 my_plot.pdf")
}

input_file <- args[1]
window_size <- as.numeric(args[2])
output_file <- args[3]

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

message("Loading libraries...")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(gridExtra)
})

message(paste("Reading", input_file))
df <- read_tsv(input_file, show_col_types = FALSE)

# Filter by window
message(paste("Filtering data to +/-", window_size, "bp"))
df_filtered <- df %>%
  filter(Position >= -window_size & Position <= window_size)

# Set up colors
# Match Python: Lost=C1 (Orange), Preserved=C0 (Blue), Constitutive=C2 (Green)
# R default palette is differnt, so let's set manual colors close to matplotlib's C0/C1/C2
custom_colors <- c(
  "Lost" = "#ff7f0e",       # C1 Orange
  "Preserved" = "#1f77b4",  # C0 Blue
  "Constitutive" = "#2ca02c", # C2 Green
  "UFM1_dependent" = "#ff7f0e",
  "UFM1_independent" = "#1f77b4"
)

message("Plotting...")

# Start Codon Plot
p_start <- ggplot(df_filtered %>% filter(Feature == "Start_Codon"), 
                  aes(x = Position, y = GC_Percent, color = Group)) +
  geom_line(size = 1) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(title = paste("Start Codon Anchored GC Content (+/-", window_size, "bp)"),
       y = "GC Fraction", x = "Position relative to Start Codon") +
  theme(legend.position = "bottom")

# Stop Codon Plot
p_stop <- ggplot(df_filtered %>% filter(Feature == "Stop_Codon"), 
                 aes(x = Position, y = GC_Percent, color = Group)) +
  geom_line(size = 1) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(title = paste("Stop Codon Anchored GC Content (+/-", window_size, "bp)"),
       y = "GC Fraction", x = "Position relative to Stop Codon") +
  theme(legend.position = "bottom")

# Combine and Save
message(paste("Saving to", output_file))
p_combined <- grid.arrange(p_start, p_stop, ncol = 2)
ggsave(output_file, p_combined, width = 12, height = 6)

message("Done.")
