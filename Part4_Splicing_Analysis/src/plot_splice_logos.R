#!/usr/bin/env Rscript

# plot_splice_logos.R
# Visualizes Splice Site Logos (Lost vs Constitutive) and SRSF5 Motif Hits
# Usage: Rscript plot_splice_logos.R <lost_5ss> <lost_3ss> <lost_srsf5> <pres_5ss> <pres_3ss> <pres_srsf5> <out_pdf>

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 7) {
  stop("Usage: Rscript plot_splice_logos.R <lost_5ss> <lost_3ss> <lost_srsf5> <pres_5ss> <pres_3ss> <pres_srsf5> <out_pdf>")
}

# Lost Inputs
l_5ss_file <- args[1]
l_3ss_file <- args[2]
l_srsf5_file <- args[3]

# Preserved Inputs
p_5ss_file <- args[4]
p_3ss_file <- args[5]
p_srsf5_file <- args[6]

out_pdf <- args[7]

# Load Libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggseqlogo)
    library(Biostrings)
    library(gridExtra)
    library(grid) # Required for textGrob/gpar
})

# Helper: Read FASTA
get_sequences <- function(fasta_file) {
    if (!file.exists(fasta_file)) return(NULL)
    seqs <- readDNAStringSet(fasta_file)
    as.character(seqs)
}

# Load Data
l_5ss_seq <- get_sequences(l_5ss_file)
l_3ss_seq <- get_sequences(l_3ss_file)
l_srsf5_seq <- get_sequences(l_srsf5_file)

p_5ss_seq <- get_sequences(p_5ss_file)
p_3ss_seq <- get_sequences(p_3ss_file)
p_srsf5_seq <- get_sequences(p_srsf5_file)

# Helper: Slice Strings
slice_seqs <- function(seqs, start, end) {
    sapply(seqs, function(x) substr(x, start, end), USE.NAMES=FALSE)
}

# Process 5'SS (User wants positions 2-8)
# MaxEntScan 5'SS is 9-mer. 2-8 corresponds to indices 2:8.
l_5ss_plot <- slice_seqs(l_5ss_seq, 2, 8)
p_5ss_plot <- slice_seqs(p_5ss_seq, 2, 8)

# Process 3'SS (User wants "AG" region)
# MaxEntScan 3'SS is 23-mer. AG is at 19-20.
# Let's show 10-23 (Py tract end + AG + 3 intron/exon jxn bases)
l_3ss_plot <- slice_seqs(l_3ss_seq, 10, 23)
p_3ss_plot <- slice_seqs(p_3ss_seq, 10, 23)

# Process SRSF5 (Full 5-mer hits)
l_srsf5_plot <- l_srsf5_seq
p_srsf5_plot <- p_srsf5_seq

# Generate Plots

# Row 1: Lost
pl1 <- ggseqlogo(l_5ss_plot) + ggtitle("Lost: 5'SS (Pos 2-8)") + theme(axis.text.x=element_blank())
pl2 <- ggseqlogo(l_srsf5_plot) + ggtitle("Lost: SRSF5 Hit (Best Match)") + theme(axis.text.x=element_blank())
pl3 <- ggseqlogo(l_3ss_plot) + ggtitle("Lost: 3'SS (AG Region)") + theme(axis.text.x=element_blank())

# Row 2: Preserved
pp1 <- ggseqlogo(p_5ss_plot) + ggtitle("Preserved: 5'SS (Pos 2-8)") + theme(axis.text.x=element_blank())
pp2 <- ggseqlogo(p_srsf5_plot) + ggtitle("Preserved: SRSF5 Hit (Best Match)") + theme(axis.text.x=element_blank())
pp3 <- ggseqlogo(p_3ss_plot) + ggtitle("Preserved: 3'SS (AG Region)") + theme(axis.text.x=element_blank())

# Arrange
pdf(out_pdf, width=12, height=6)
grid.arrange(
    pl1, pl2, pl3,
    pp1, pp2, pp3,
    ncol=3, nrow=2,
    top = textGrob("Splicing Signal Comparison: Lost vs Preserved (Human)", gp=gpar(fontsize=15, font=2))
)
dev.off()

print(paste("Generated:", out_pdf))
