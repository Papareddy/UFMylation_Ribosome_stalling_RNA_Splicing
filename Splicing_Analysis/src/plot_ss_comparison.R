#!/usr/bin/env Rscript

# plot_ss_comparison.R
# Visualizes Splice Site Logos: Constitutive vs Preserved vs Lost
# Comparison of 3 Groups to show "Weakness" relative to Constitutive.

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 7) {
  stop("Usage: Rscript plot_ss_comparison.R <const_5ss> <const_3ss> <pres_5ss> <pres_3ss> <lost_5ss> <lost_3ss> <out_pdf>")
}

c_5ss_file <- args[1]
c_3ss_file <- args[2]
p_5ss_file <- args[3]
p_3ss_file <- args[4]
l_5ss_file <- args[5]
l_3ss_file <- args[6]
out_pdf <- args[7]

# Load Libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggseqlogo)
    library(Biostrings)
    library(gridExtra)
    library(grid)
})

# Helper
get_sequences <- function(fasta_file) {
    if (!file.exists(fasta_file)) return(NULL)
    seqs <- readDNAStringSet(fasta_file)
    as.character(seqs)
}

# Load Data
c_5ss <- get_sequences(c_5ss_file)
c_3ss <- get_sequences(c_3ss_file)
p_5ss <- get_sequences(p_5ss_file)
p_3ss <- get_sequences(p_3ss_file)
l_5ss <- get_sequences(l_5ss_file)
l_3ss <- get_sequences(l_3ss_file)

# Helper: Slice
slice_seqs <- function(seqs, start, end) {
    sapply(seqs, function(x) substr(x, start, end), USE.NAMES=FALSE)
}

# 5'SS: Pos 2-8 (GT + Flank)
c_5ss_plot <- slice_seqs(c_5ss, 2, 8)
p_5ss_plot <- slice_seqs(p_5ss, 2, 8)
l_5ss_plot <- slice_seqs(l_5ss, 2, 8)

# 3'SS: Pos 10-23 (Py Tract End + AG + Exon)
c_3ss_plot <- slice_seqs(c_3ss, 10, 23)
p_3ss_plot <- slice_seqs(p_3ss, 10, 23)
l_3ss_plot <- slice_seqs(l_3ss, 10, 23)

# Generate Plots
# Row 1: Constitutive
p1_c <- ggseqlogo(c_5ss_plot) + ggtitle("Constitutive: Strong 5'SS (Donor)") + theme(axis.text.x=element_blank())
p2_c <- ggseqlogo(c_3ss_plot) + ggtitle("Constitutive: Strong 3'SS (Acceptor)") + theme(axis.text.x=element_blank())

# Row 2: Preserved
p1_p <- ggseqlogo(p_5ss_plot) + ggtitle("Preserved: Weak 5'SS (Donor)") + theme(axis.text.x=element_blank())
p2_p <- ggseqlogo(p_3ss_plot) + ggtitle("Preserved: Weak 3'SS (Acceptor)") + theme(axis.text.x=element_blank())

# Row 3: Lost
p1_l <- ggseqlogo(l_5ss_plot) + ggtitle("Lost: Weak 5'SS (Donor)") + theme(axis.text.x=element_blank())
p2_l <- ggseqlogo(l_3ss_plot) + ggtitle("Lost: Weak 3'SS (Acceptor)") + theme(axis.text.x=element_blank())

# Arrange
pdf(out_pdf, width=10, height=8)
grid.arrange(
    p1_c, p2_c,
    p1_p, p2_p,
    p1_l, p2_l,
    ncol=2, nrow=3,
    top = textGrob("Splice Site Strength Comparison: Constitutive vs Preserved vs Lost", gp=gpar(fontsize=15, font=2))
)
dev.off()

print(paste("Generated:", out_pdf))
