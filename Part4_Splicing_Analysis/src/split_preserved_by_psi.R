
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
out_dir <- args[2]

if (is.na(input_file) || is.na(out_dir)) {
  stop("Usage: Rscript split_preserved_by_psi.R <preserved.tsv> <out_dir>")
}

if (!file.exists(input_file)) {
  stop(paste("File not found:", input_file))
}

df <- read_tsv(input_file, show_col_types = FALSE)

# Check for IncLevel1 (WT PSI)
# rMATS cols: IncLevel1, IncLevel2. usually comma separated. 
# We need mean PSI.
parse_psi <- function(x) {
  vals <- as.numeric(unlist(strsplit(as.character(x), ",")))
  return(mean(vals, na.rm=TRUE))
}

# Add mean PSI column
df$Mean_WT_PSI <- sapply(df$IncLevel1, parse_psi)

# Split
# Included: PSI > 0.5
# Excluded: PSI <= 0.5
included <- df %>% filter(Mean_WT_PSI > 0.5)
excluded <- df %>% filter(Mean_WT_PSI <= 0.5)

write_tsv(included, file.path(out_dir, "preserved_included.tsv"))
write_tsv(excluded, file.path(out_dir, "preserved_excluded.tsv"))

cat(sprintf("Split Preserved: %d Included (PSI>0.5), %d Excluded (PSI<=0.5)\n", nrow(included), nrow(excluded)))
