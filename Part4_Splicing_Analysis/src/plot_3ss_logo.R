#!/usr/bin/env Rscript
library(ggseqlogo)
library(ggplot2)
library(Biostrings)
library(optparse)

option_list <- list(
  make_option(c("--dep"), type="character", help="Path to UFM1_dependent FASTA"),
  make_option(c("--indep"), type="character", help="Path to UFM1_independent FASTA"),
  make_option(c("--const"), type="character", help="Path to Constitutive FASTA"),
  make_option(c("--out"), type="character", help="Output PDF path")
)

opt <- parse_args(OptionParser(option_list=option_list))

load_seqs <- function(fasta_path) {
  if (!file.exists(fasta_path)) return(NULL)
  s <- readDNAStringSet(fasta_path)
  return(as.character(s))
}

data_list <- list()

seq_dep <- load_seqs(opt$dep)
if (!is.null(seq_dep)) data_list[["UFM1-dependent"]] <- seq_dep

seq_indep <- load_seqs(opt$indep)
if (!is.null(seq_indep)) data_list[["UFM1-independent"]] <- seq_indep

seq_const <- load_seqs(opt$const)
if (!is.null(seq_const)) data_list[["Constitutive"]] <- seq_const

if (length(data_list) == 0) {
  cat("No sequences loaded. Creating empty output file.\n")
  pdf(opt$out)
  plot.new()
  text(0.5, 0.5, "No sequences available for analysis")
  dev.off()
  quit(save="no", status=0)
}

# Ensure all groups have at least one sequence to avoid ggseqlogo error
data_list <- data_list[sapply(data_list, length) > 0]

if (length(data_list) == 0) {
  cat("All groups are empty. Creating empty output file.\n")
  pdf(opt$out)
  plot.new()
  text(0.5, 0.5, "No sequences available for analysis")
  dev.off()
  quit(save="no", status=0)
}

# Generate Plot
p <- ggseqlogo(data_list, ncol=1) +
  theme_logo() +
  scale_x_continuous(breaks=c(1, 11, 21, 31, 40), labels=c("-20", "-10", "0", "10", "20")) +
  annotate('rect', xmin = 20.5, xmax = 20.5, ymin = 0, ymax = 2, alpha = .5, col='black', linetype='dashed') +
  labs(title="Sequence Logo around 3' Splice Site (+/- 20bp)")

ggsave(opt$out, p, width=8, height=3 * length(data_list))
