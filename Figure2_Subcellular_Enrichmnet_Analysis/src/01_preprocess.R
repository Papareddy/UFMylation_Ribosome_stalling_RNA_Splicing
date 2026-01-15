# src/01_preprocess.R

# Load required libraries
library(optparse)

# Define command-line options
option_list <- list(
  make_option("--input", type="character", default=NULL, help="Input data file path", metavar="character"),
  make_option("--output", type="character", default="results/preprocessed_data.csv", help="Output data file path", metavar="character"),
  make_option("--pval", type="double", default=0.05, help="P-value threshold", metavar="number"),
  make_option("--fc_threshold", type="double", default=1.0, help="Log2 Fold-change threshold", metavar="number")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if input file is provided
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Input file must be specified.", call.=FALSE)
}

# Print parameters
cat("Input file:", opt$input, "\n")
cat("Output file:", opt$output, "\n")
cat("P-value threshold:", opt$pval, "\n")
cat("Log2 FC threshold:", opt$fc_threshold, "\n")

# --- Dummy Pre-processing ---
# In a real script, you would load the data, filter it, and save it.
# For this placeholder, we'll just create a dummy output file.
cat("Simulating data pre-processing...\n")
dummy_data <- data.frame(
  GeneID = c("GENE1", "GENE2", "GENE3"),
  log2FC = c(1.5, -1.2, 2.0),
  pvalue = c(0.01, 0.04, 0.001)
)

# Write dummy output
write.csv(dummy_data, file = opt$output, row.names = FALSE)

cat("Pre-processing complete. Output written to", opt$output, "\n")
