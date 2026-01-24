#!/usr/bin/env Rscript

# Usage: Rscript plot_custom_logo.R <meme_file> <out_file>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_custom_logo.R <meme_file> [out_file]")
}

meme_file <- args[1]
out_file <- if (length(args) >= 2) args[2] else "custom_logo.pdf"

library(ggplot2)
library(ggseqlogo)

# Function to parse MEME file (Simplified for this custom file)
parse_meme_matrix <- function(file, motif_name_pattern) {
  lines <- readLines(file)
  
  # Find motif start
  motif_idx <- grep(paste0("MOTIF ", motif_name_pattern), lines)
  if (length(motif_idx) == 0) return(NULL)
  
  # Find matrix start
  # Look for "letter-probability matrix" after motif line
  matrix_start <- -1
  for (i in (motif_idx + 1):length(lines)) {
    if (grepl("letter-probability matrix", lines[i])) {
      matrix_start <- i + 1
      break
    }
  }
  
  if (matrix_start == -1) return(NULL)
  
  # Read matrix lines until empty line or next motif
  matrix_lines <- c()
  for (i in matrix_start:length(lines)) {
    if (trimws(lines[i]) == "" || grepl("MOTIF", lines[i])) break
    matrix_lines <- c(matrix_lines, lines[i])
  }
  
  # Convert to numeric matrix
  mat_data <- do.call(rbind, lapply(matrix_lines, function(x) {
    as.numeric(strsplit(trimws(x), "\\s+")[[1]])
  }))
  
  # Transpose for ggseqlogo (rows = bases, cols = positions)
  # MEME is usually ACGT columns, positions rows -> Transpose to Bases x Rows
  return(t(mat_data))
}

# 1. Parse SRSF3
pwm_srsf3 <- parse_meme_matrix(meme_file, "SRSF3_SRp20")

if (!is.null(pwm_srsf3)) {
  rownames(pwm_srsf3) <- c("A", "C", "G", "T") # Assuming ACGT alphabet
  
  p <- ggseqlogo(pwm_srsf3) +
    ggtitle("SRSF3_SRp20 (Custom)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(out_file, p, width = 6, height = 4)
  print(paste("Saved logo to", out_file))
} else {
  stop("Motif SRSF3_SRp20 not found in file.")
}
