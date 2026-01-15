# src/03_visualize.R

# Load required libraries
library(optparse)
library(ggplot2)

# Define command-line options
option_list <- list(
  make_option("--input", type="character", default=NULL, help="Input enrichment results file path", metavar="character"),
  make_option("--output", type="character", default="results/enrichment_plot.png", help="Output plot file path", metavar="character"),
  make_option("--top_n", type="integer", default=20, help="Number of top items to plot", metavar="number")
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
cat("Top N:", opt$top_n, "\n")

# --- Dummy Visualization ---
# In a real script, you would create a meaningful plot from the enrichment data.
# For this placeholder, we'll create a dummy plot.
cat("Simulating visualization...\n")

# Read dummy data (or real data in a real script)
enrichment_data <- read.csv(opt$input)

# Create a simple bar plot
p <- ggplot(enrichment_data, aes(x = reorder(Compartment, -log10(p.value)), y = -log10(p.value))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(
    title = "Compartment Enrichment Analysis",
    subtitle = paste("Top", opt$top_n, "compartments"),
    x = "Compartment",
    y = "-log10(p-value)"
  ) +
  theme_minimal()

# Save dummy plot
ggsave(opt$output, plot = p, width = 8, height = 6)

cat("Visualization complete. Plot saved to", opt$output, "\n")
