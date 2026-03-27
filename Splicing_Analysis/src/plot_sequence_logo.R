#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggseqlogo)
  library(ggplot2)
  library(Biostrings)
  library(argparse)
  library(cowplot)
})

# Argument Parser
parser <- ArgumentParser(description="Generate Positional Sequence Logos")
parser$add_argument("--fasta5", required=TRUE, help="FASTA file for 5'SS (Donor)")
parser$add_argument("--fasta3", required=TRUE, help="FASTA file for 3'SS (Acceptor)")
parser$add_argument("--outdir", required=TRUE, help="Output directory")
parser$add_argument("--prefix", default="Logo", help="Output filename prefix")
args <- parser$parse_args()

# Function to read and prep sequences
get_seqs <- function(fasta_file) {
  if (!file.exists(fasta_file)) return(NULL)
  seqs <- readDNAStringSet(fasta_file)
  # Convert to character vector
  as.character(seqs)
}

# Main Plotting Logic
do_logo <- function(seqs, title, x_start=-100) {
  if (is.null(seqs) || length(seqs) == 0) return(NULL)
  
  # Validate Length (assume all same, check first)
  len <- nchar(seqs[1])
  
  # Generate Logo
  # ggseqlogo method='bits' (Information Content) or 'prob' (Probability)
  # User asked for "probability (height = information content or frequency)"
  # Usually 'bits' is standard for "height = information". 'prob' is frequency.
  # Let's use 'bits' as it highlights conservation better vs background.
  
  p <- ggseqlogo(seqs, method='bits') +
    scale_x_continuous(breaks = seq(1, len, by=20), 
                       labels = seq(x_start, x_start + len - 1, by=20)) +
    theme_classic() +
    labs(title=title, x="Distance from Splice Site (bp)", y="Bits") +
    geom_vline(xintercept = -x_start + 0.5, linetype="dashed", color="black") # Center line
    
  return(p)
}

# 1. Load Data
s5 <- get_seqs(args$fasta5)
s3 <- get_seqs(args$fasta3)

# 2. Plot 5'SS
# Assuming centered at window start? No.
# plot_rna_map used +/- 100 window.
# So center is exactly middle.
# Length should be 200.
# Left is -100 (Exon). Right is +100 (Intron).
# Junction at 100|101.
# x-axis logic:
# Index 1 = -100.
# Index 100 = -1.
# Index 101 = +1 (skip 0).
# Index 200 = +100.

# Function to generate custom labels skipping 0
make_labels <- function(len, half) {
  # Indices 1..len
  # 1 maps to -half
  # half maps to -1
  # half+1 maps to +1
  # len maps to +half
  
  # Just generate the vector
  lab <- c(seq(-half, -1), seq(1, half))
  return(lab)
}

plot_logo_custom <- function(seqs, title, half_window=100) {
  if (is.null(seqs)) return(ggplot() + theme_void() + labs(title=paste(title, "(No Data)")))
  
  p <- ggseqlogo(seqs, method='bits')
  
  # Custom Axis
  # Label every 20 bp
  breaks <- seq(1, 2*half_window, by=20)
  # Map breaks to genomic relative coords
  # 1 -> -100
  # 21 -> -80
  # ...
  # 101 -> +1
  
  labels <- c()
  for (b in breaks) {
    if (b <= half_window) {
      labels <- c(labels, b - half_window - 1) # 1 - 101 = -100
    } else {
      labels <- c(labels, b - half_window)   # 101 - 100 = 1
    }
  }
  
  p <- p + 
    scale_x_continuous(breaks=breaks, labels=labels) +
    geom_vline(xintercept = half_window + 0.5, linetype="dashed", color="black") +
    theme_cowplot() +
    labs(title=title, x="Distance to Splice Site (bp)", y="Information (Bits)")
  
  return(p)
}


# Auto-detect window size from first sequence
# Length is 2 * half_window
len_seq <- nchar(s5[1])
half_window <- len_seq / 2
cat(paste("[INFO] Detected Sequence Length:", len_seq, "-> Half Window:", half_window, "bp\n"))

p5 <- plot_logo_custom(s5, "5'SS (Donor) - UFM1-dependent (Lost)", half_window)
p3 <- plot_logo_custom(s3, "3'SS (Acceptor) - UFM1-dependent (Lost)", half_window)

# Combine
final_plot <- plot_grid(p5, p3, ncol=1, align='v')

outfile <- file.path(args$outdir, paste0(args$prefix, "_SequenceLogo.pdf"))
ggsave(outfile, final_plot, width=10, height=8)
cat(paste("[INFO] Logo saved to", outfile, "\n"))
