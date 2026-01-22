#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT: Anchored Stop Codon Intron Density Analysis
# PURPOSE: Count normalized intron density anchored at stop codons
#          for UFM1-dependent, UFM1-independent, and constitutive introns
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(tidyverse)
  library(optparse)
})

# ==============================================================================
# 1. ARGUMENT PARSING
# ==============================================================================
option_list <- list(
  make_option(c("--gtf"), type="character", help="Path to GTF file"),
  make_option(c("--events"), type="character", help="Path to UFM1 events RDS file"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory"),
  make_option(c("--window"), type="integer", default=1000, help="Window size (+/- bp) around stop codon")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$gtf) || is.null(opt$events)) {
  stop("Missing required arguments: --gtf, --events")
}

if(!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

WIN <- opt$window

message("--- Loading Data ---")
txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")
events <- readRDS(opt$events)

# Harmonize chromosome styles
tryCatch({
  tx_gr <- transcripts(txdb)
  if (!any(seqlevels(events) %in% seqlevels(tx_gr))) {
    message("Harmonizing chromosome naming styles...")
    seqlevelsStyle(events) <- seqlevelsStyle(tx_gr)[1]
  }
}, error=function(e) message("Style harmonization skipped"))

# ==============================================================================
# 2. DEFINE GROUPS
# ==============================================================================
message("--- Defining Groups ---")

# A. Experimental Groups (RI events only)
ri_events <- events[events$EventType == "RI", ]
dep_ri <- ri_events[ri_events$Group == "UFM1_dependent", ]
indep_ri <- ri_events[ri_events$Group == "UFM1_independent", ]

# B. Constitutive Introns (all introns in 3'UTR, excluding experimental)
all_introns_grl <- intronsByTranscript(txdb, use.names=TRUE)
all_introns_flat <- unlist(all_introns_grl, use.names=TRUE)

utr3_grl <- threeUTRsByTranscript(txdb, use.names=TRUE)
utr3_ranges <- unlist(range(utr3_grl), use.names=TRUE)

utr3_introns <- subsetByOverlaps(all_introns_flat, utr3_ranges, type="within")

# Remove experimental introns from control
targets <- c(dep_ri, indep_ri)
hits <- findOverlaps(utr3_introns, targets)
control_introns <- if(length(hits) > 0) utr3_introns[-unique(queryHits(hits))] else utr3_introns

message(paste("N Dependent:", length(dep_ri)))
message(paste("N Independent:", length(indep_ri)))
message(paste("N Control:", length(control_introns)))

# ==============================================================================
# 3. GET STOP CODON POSITIONS
# ==============================================================================
message("--- Extracting Stop Codon Positions ---")

# Get CDS by transcript
cds_grl <- cdsBy(txdb, by="tx", use.names=TRUE)

# Calculate stop codon position for each transcript
get_stop_codon_positions <- function(cds_grl) {
  tx_ids <- names(cds_grl)
  stop_pos <- sapply(tx_ids, function(tx) {
    cds <- cds_grl[[tx]]
    if (length(cds) == 0) return(NA)
    strand_val <- as.character(strand(cds)[1])
    if (strand_val == "+") {
      return(max(end(cds)))
    } else {
      return(min(start(cds)))
    }
  })
  
  # Get chromosome and strand for each transcript
  chr_val <- sapply(tx_ids, function(tx) as.character(seqnames(cds_grl[[tx]])[1]))
  strand_val <- sapply(tx_ids, function(tx) as.character(strand(cds_grl[[tx]])[1]))
  
  data.frame(
    tx_id = tx_ids,
    chr = chr_val,
    stop_pos = as.numeric(stop_pos),
    strand = strand_val,
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(stop_pos))
}

stop_codon_df <- get_stop_codon_positions(cds_grl)
message(paste("Found stop codons for", nrow(stop_codon_df), "transcripts"))

# ==============================================================================
# 4. CALCULATE INTRON DENSITY RELATIVE TO STOP CODON
# ==============================================================================
message("--- Calculating Intron Density ---")

# OPTIMIZATION: Get transcripts once instead of inside function
message("Loading transcript annotations (this may take a moment)...")
tx_gr <- transcripts(txdb)
message(paste("Loaded", length(tx_gr), "transcripts"))

calculate_density <- function(gr, group_name, stop_df, window, tx_granges) {
  if (length(gr) == 0) return(data.frame())
  
  # Map introns to transcripts (using preloaded tx_granges)
  hits <- findOverlaps(gr, tx_granges)
  if (length(hits) == 0) return(data.frame())
  
  intron_tx_map <- data.frame(
    intron_idx = queryHits(hits),
    tx_id = tx_granges$tx_name[subjectHits(hits)]
  )
  
  # Join with stop codon positions
  merged <- inner_join(intron_tx_map, stop_df, by="tx_id")
  if (nrow(merged) == 0) return(data.frame())
  
  # Calculate intron center relative to stop codon
  intron_centers <- (start(gr)[merged$intron_idx] + end(gr)[merged$intron_idx]) / 2
  
  # Calculate distance: positive = downstream of stop, negative = upstream
  merged$distance <- ifelse(merged$strand == "+",
                            intron_centers - merged$stop_pos,
                            merged$stop_pos - intron_centers)
  
  # Filter to window
  merged_filt <- merged %>% filter(abs(distance) <= window)
  
  if (nrow(merged_filt) == 0) return(data.frame())
  
  # Create density using histogram bins (25bp bins as requested)
  breaks <- seq(-window, window, by=25)
  counts <- hist(merged_filt$distance, breaks=breaks, plot=FALSE)$counts
  mids <- breaks[-length(breaks)] + 12.5  # Center of 25bp bins
  
  # Return both raw counts and normalized density
  data.frame(
    Position = mids,
    RawCounts = counts,
    Density = counts / length(gr),
    Group = group_name
  )
}

# Calculate for each group (using preloaded transcripts)
density_dep <- calculate_density(dep_ri, "UFM1_dependent", stop_codon_df, WIN, tx_gr)
density_indep <- calculate_density(indep_ri, "UFM1_independent", stop_codon_df, WIN, tx_gr)
density_control <- calculate_density(control_introns, "Constitutive", stop_codon_df, WIN, tx_gr)

all_density <- bind_rows(density_dep, density_indep, density_control)

# ==============================================================================
# 5. PLOTTING
# ==============================================================================
message("--- Generating Plot ---")

if (nrow(all_density) > 0) {
  # Define colors
  cols <- c("UFM1_dependent" = "#E57373", 
            "UFM1_independent" = "#64B5F6", 
            "Constitutive" = "#9E9E9E")
  
  # Create counts for legend
  leg_labels <- c(
    paste0("UFM1 Dep (n=", length(dep_ri), ")"),
    paste0("UFM1 Indep (n=", length(indep_ri), ")"),
    paste0("Constitutive (n=", length(control_introns), ")")
  )
  
  # Calculate max y for plot
  max_y <- max(all_density$Density, na.rm=TRUE) * 1.1
  
  # Open PDF
  pdf(file.path(opt$outdir, "Anchored_StopCodon_IntronDensity.pdf"), width=8, height=6)
  
  par(mar=c(5, 5, 4, 2))
  
  plot(1, type="n", xlim=c(-WIN, WIN), ylim=c(0, max_y),
       xlab="Distance from Stop Codon (bp)", 
       ylab="Normalized Intron Density",
       main="Intron Density Anchored at Stop Codon",
       bty="l", las=1)
  
  # Shade upstream (CDS) and downstream (3'UTR) regions
  rect(-WIN, 0, 0, max_y*1.2, col=rgb(0.9,0.9,1,0.3), border=NA)
  rect(0, 0, WIN, max_y*1.2, col=rgb(1,0.9,0.9,0.3), border=NA)
  
  abline(v=0, lty=2, lwd=2)
  text(-WIN/2, max_y*0.95, "CDS", col="blue", cex=0.9)
  text(WIN/2, max_y*0.95, "3'UTR", col="red", cex=0.9)
  
  # Plot lines
  for (g in c("Constitutive", "UFM1_independent", "UFM1_dependent")) {
    dat <- all_density[all_density$Group == g, ]
    if (nrow(dat) > 0) {
      lines(dat$Position, dat$Density, col=cols[g], lwd=2)
    }
  }
  
  legend("topright", legend=leg_labels, col=cols, lwd=2, bty="n", cex=0.9)
  
  dev.off()
  
  # Save normalized density data
  write_tsv(all_density, file.path(opt$outdir, "Anchored_StopCodon_IntronDensity.tsv"))
  
  # Save raw counts only (not normalized)
  raw_counts_only <- all_density %>% select(Position, RawCounts, Group)
  write_tsv(raw_counts_only, file.path(opt$outdir, "Anchored_StopCodon_RawCounts.tsv"))
  
  message("Done. Saved to: ", opt$outdir)
} else {
  message("[WARN] No density data to plot.")
}
