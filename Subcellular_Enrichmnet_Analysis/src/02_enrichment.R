# src/02_enrichment.R

# Load required libraries
library(optparse)

# Define command-line options
option_list <- list(
  make_option("--input", type="character", default=NULL, help="Input pre-processed data file path", metavar="character"),
  make_option("--output", type="character", default="results/enrichment_results.csv", help="Output enrichment file path", metavar="character")
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

# --- Dummy Enrichment Analysis ---
# In a real script, you would perform Fisher's Exact Test or similar.
# For this placeholder, we'll create a dummy output file.
cat("Simulating enrichment analysis...\n")
dummy_enrichment <- data.frame(
  Compartment = c("Nucleus", "Cytosol", "Ribosome"),
  p.value = c(0.001, 0.045, 0.1),
  enriched_genes = c("GENE3", "GENE1", "GENE2")
)

# Write dummy output
write.csv(dummy_enrichment, file = opt$output, row.names = FALSE)

cat("Enrichment analysis complete. Output written to", opt$output, "\n")
