#!/usr/bin/env Rscript
#' ese_motif_extraction.R
#' Step A: Sequence Extraction using GRanges
#' 
#' Uses GRanges directly for ROI extraction from RDS files.
#' Generates BED and FASTA files for motif scanning.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(Biostrings)
  library(rtracklayer)
})

# === ARGUMENTS ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript ese_motif_extraction.R <events_rds> <genome_fasta> <gtf> <outdir>")
}

RDS_PATH <- args[1]
GENOME_FASTA <- args[2]
GTF_PATH <- args[3]
OUTDIR <- args[4]

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

message("[INFO] Loading events from RDS...")
events <- readRDS(RDS_PATH)

# Filter for RI events only
ri_events <- events[events$EventType == "RI"]
message(paste("[INFO] Total RI events:", length(ri_events)))

# Split by group
dep_ri <- ri_events[ri_events$Group == "UFM1_dependent"]
indep_ri <- ri_events[ri_events$Group == "UFM1_independent"]
message(paste("[INFO] UFM1-Dependent:", length(dep_ri), ", UFM1-Independent:", length(indep_ri)))

# === HELPER FUNCTIONS ===

extract_rois <- function(gr, roi_name) {
  #' Extract ROIs relative to intron boundaries
  #' For RI events: the GRanges represents the retained intron
  #' 5'SS (donor) is at start, 3'SS (acceptor) is at end
  #' Using ±50bp windows for better ESE context
  
  if (length(gr) == 0) return(list())
  
  # ROI_1: 5'SS ±50bp (centered on donor = start of intron)
  roi1 <- GRanges(
    seqnames = seqnames(gr),
    ranges = IRanges(start = start(gr) - 50, end = start(gr) + 50),
    strand = strand(gr)
  )
  mcols(roi1)$name <- paste0(roi_name, "_", seq_along(gr), "_roi1_5ss")
  
  # ROI_2: 3'SS ±50bp (centered on acceptor = end of intron)
  roi2 <- GRanges(
    seqnames = seqnames(gr),
    ranges = IRanges(start = end(gr) - 50, end = end(gr) + 50),
    strand = strand(gr)
  )
  mcols(roi2)$name <- paste0(roi_name, "_", seq_along(gr), "_roi2_3ss")

  
  # ROI_3: Branch Point -15 to -45 upstream of 3'SS
  # For + strand: [end - 45, end - 15] (upstream of acceptor)
  # For - strand: [end + 15, end + 45] (upstream of acceptor in transcript coords)
  plus_idx <- which(as.character(strand(gr)) == "+")
  minus_idx <- which(as.character(strand(gr)) == "-")
  
  roi3_list <- list()
  if (length(plus_idx) > 0) {
    roi3_plus <- GRanges(
      seqnames = seqnames(gr)[plus_idx],
      ranges = IRanges(start = end(gr)[plus_idx] - 45, end = end(gr)[plus_idx] - 15),
      strand = strand(gr)[plus_idx]
    )
    mcols(roi3_plus)$name <- paste0(roi_name, "_", plus_idx, "_roi3_bp")
    roi3_list <- c(roi3_list, roi3_plus)
  }
  
  if (length(minus_idx) > 0) {
    roi3_minus <- GRanges(
      seqnames = seqnames(gr)[minus_idx],
      ranges = IRanges(start = end(gr)[minus_idx] + 15, end = end(gr)[minus_idx] + 45),
      strand = strand(gr)[minus_idx]
    )
    mcols(roi3_minus)$name <- paste0(roi_name, "_", minus_idx, "_roi3_bp")
    roi3_list <- c(roi3_list, roi3_minus)
  }
  
  roi3 <- do.call(c, roi3_list)
  if (is.null(roi3)) roi3 <- GRanges()

  
  return(list(roi1_5ss = roi1, roi2_3ss = roi2, roi3_bp = roi3))
}

write_bed <- function(gr, filepath) {
  if (length(gr) == 0) return()
  
  # Normalize chromosome names for Ensembl FASTA (strip chr prefix)
  chrom_names <- as.character(seqnames(gr))
  chrom_names <- sub("^chr", "", chrom_names)
  chrom_names <- ifelse(chrom_names == "M", "MT", chrom_names)
  
  bed_df <- data.frame(
    chrom = chrom_names,
    start = start(gr) - 1,  # BED is 0-based
    end = end(gr),
    name = mcols(gr)$name,
    score = 0,
    strand = as.character(strand(gr))
  )
  write.table(bed_df, filepath, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  message(paste("[INFO] Wrote", nrow(bed_df), "entries to", basename(filepath)))
}


extract_fasta <- function(bed_file, genome_fasta, output_fasta) {
  cmd <- paste("bedtools getfasta -s -fi", genome_fasta, "-bed", bed_file, "-fo", output_fasta, "-name")
  system(cmd)
  message(paste("[INFO] Extracted sequences to", basename(output_fasta)))
}

# === LOAD CONSTITUTIVE INTRONS ===
message("[INFO] Loading constitutive introns from GTF...")
txdb <- makeTxDbFromGFF(GTF_PATH, format = "gtf")

# Get all introns from TxDb
all_introns <- intronsByTranscript(txdb, use.names = TRUE)
all_introns_gr <- unlist(all_introns)

# Remove any introns that overlap with RI events (to get truly constitutive)
ri_all <- c(dep_ri, indep_ri)
constitutive <- all_introns_gr[!overlapsAny(all_introns_gr, ri_all)]

# Sample to reasonable size
if (length(constitutive) > 500) {
  set.seed(42)
  constitutive <- constitutive[sample(length(constitutive), 500)]
}
message(paste("[INFO] Constitutive introns (control):", length(constitutive)))

# === PROCESS EACH GROUP ===
groups <- list(
  UFM1_dependent = dep_ri,
  UFM1_independent = indep_ri,
  Constitutive = constitutive
)

for (group_name in names(groups)) {
  gr <- groups[[group_name]]
  
  if (length(gr) == 0) {
    message(paste("[WARN] Skipping", group_name, "(no events)"))
    next
  }
  
  message(paste("[INFO] Processing", group_name, "(n =", length(gr), ")..."))
  
  # Extract ROIs
  rois <- extract_rois(gr, group_name)
  
  # Write BED files and extract FASTA
  for (roi_name in names(rois)) {
    roi_gr <- rois[[roi_name]]
    if (length(roi_gr) == 0) next
    
    bed_path <- file.path(OUTDIR, paste0(group_name, ".", roi_name, ".bed"))
    fa_path <- file.path(OUTDIR, paste0(group_name, ".", roi_name, ".fa"))
    
    write_bed(roi_gr, bed_path)
    extract_fasta(bed_path, GENOME_FASTA, fa_path)
  }
}

message("[INFO] Step A complete. FASTA files generated in:", OUTDIR)
