#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT: Combined NMD Architecture Pipeline (EJC + PTC) - Base R
# PURPOSE: 1. Calculate Neighbor EJC Density (Self-Subtracted)
#          2. Calculate Stop Codon Density at Splice Sites
#          3. Plot all 4 views (5'/3' EJC, 5'/3' PTC) with Sample Counts
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
  library(rtracklayer)
  library(tidyverse) # Used for efficient joins/data wrangling
  library(optparse)
})

# ==============================================================================
# 1. ARGUMENT PARSING
# ==============================================================================
option_list <- list(
  make_option(c("--gtf"), type="character", help="Path to GTF file"),
  make_option(c("--fasta"), type="character", help="Path to genome FASTA file"),
  make_option(c("--events"), type="character", help="Path to UFM1 events RDS file"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory"),
  make_option(c("--ejc_window"), type="integer", default=5000, help="EJC Window (+/- bp)"),
  make_option(c("--stop_window"), type="integer", default=25, help="Stop Codon Window (+/- bp)")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$gtf) || is.null(opt$fasta) || is.null(opt$events)) {
  stop("Missing required arguments: --gtf, --fasta, --events")
}

if(!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

# Constants
WIN_EJC  <- opt$ejc_window
WIN_STOP <- opt$stop_window

message("--- Loading Data ---")
txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")
events <- readRDS(opt$events)
genome_fa <- FaFile(opt$fasta)

# ==============================================================================
# 2. GROUP DEFINITION & CONTROL SUBTRACTION
# ==============================================================================
message("--- Defining Groups ---")

# A. Experimental Groups
ri_events <- events[events$EventType == "RI", ]
dep_ri    <- ri_events[ri_events$Group == "UFM1_dependent", ]
indep_ri  <- ri_events[ri_events$Group == "UFM1_independent", ]

# B. Control Group (Subtraction Logic)
all_introns_grl  <- intronsByTranscript(txdb, use.names=TRUE)
all_introns_flat <- unlist(all_introns_grl, use.names=TRUE)
mcols(all_introns_flat)$intron_uid <- paste0("intron_", seq_along(all_introns_flat))

utr3_grl     <- threeUTRsByTranscript(txdb, use.names=TRUE)
utr3_ranges  <- unlist(range(utr3_grl), use.names=TRUE)

# Harmonize chromosome naming styles
tryCatch({
  if (!any(seqlevels(all_introns_flat) %in% seqlevels(utr3_ranges))) {
    seqlevelsStyle(all_introns_flat) <- seqlevelsStyle(utr3_ranges)[1]
  }
}, error=function(e) message("Style harmonization skipped"))

utr3_introns <- subsetByOverlaps(all_introns_flat, utr3_ranges, type="within")

# Remove overlaps with targets
targets <- c(dep_ri, indep_ri)
hits <- findOverlaps(utr3_introns, targets)
control_introns <- if(length(hits) > 0) utr3_introns[-unique(queryHits(hits))] else utr3_introns

# C. Calculate Counts for Legend
counts_info <- list(
  "UFM1_dependent"       = length(dep_ri),
  "UFM1_independent"     = length(indep_ri),
  "Control 3'UTR Intron" = length(control_introns)
)

message(paste("N Dependent:", counts_info[["UFM1_dependent"]]))
message(paste("N Independent:", counts_info[["UFM1_independent"]]))
message(paste("N Control:", counts_info[["Control 3'UTR Intron"]]))

# ==============================================================================
# 3. ANALYSIS 1: EJC DENSITY (Self-Subtracted)
# ==============================================================================
get_ejc_density <- function(gr, grp_name, txdb_obj, anchor_type="3p") {
  
  if (length(gr) == 0) return(data.frame())
  
  # 1. Map RIs to Transcripts
  if (is.null(names(gr))) {
    tx_gr <- transcripts(txdb_obj)
    tryCatch({seqlevelsStyle(gr) <- seqlevelsStyle(tx_gr)[1]}, error=function(e){})
    hits <- findOverlaps(gr, tx_gr)
    if (length(hits) == 0) return(data.frame())
    map_df <- data.frame(RI_Index = queryHits(hits), tx_id = tx_gr$tx_name[subjectHits(hits)])
    expanded_gr <- gr[map_df$RI_Index]
    mcols(expanded_gr)$tx_id <- map_df$tx_id
    gr_to_use <- expanded_gr
  } else {
    gr_to_use <- gr
    mcols(gr_to_use)$tx_id <- names(gr)
  }
  
  # 2. Identify "Self" Introns to EXCLUDE
  self_hits <- findOverlaps(gr_to_use, all_introns_flat)
  forbidden_df <- data.frame(
    tx_id = mcols(gr_to_use)$tx_id[queryHits(self_hits)],
    intron_uid = mcols(all_introns_flat)$intron_uid[subjectHits(self_hits)]
  )
  forbidden_df$block_key <- paste(forbidden_df$tx_id, forbidden_df$intron_uid, sep="|")
  
  # 3. Define Anchor
  is_plus <- as.character(strand(gr_to_use)) == "+"
  if (anchor_type == "5p") {
    anchor_pos <- ifelse(is_plus, start(gr_to_use), end(gr_to_use)) # Start of Intron
  } else {
    anchor_pos <- ifelse(is_plus, end(gr_to_use), start(gr_to_use))   # End of Intron
  }
  
  # 4. Get Neighbors
  host_tx_ids <- unique(mcols(gr_to_use)$tx_id)
  relevant_introns <- all_introns_flat[names(all_introns_flat) %in% host_tx_ids]
  
  neighbors_df <- data.frame(
    tx_id = names(relevant_introns),
    intron_uid = mcols(relevant_introns)$intron_uid,
    intron_center = (start(relevant_introns) + end(relevant_introns)) / 2
  )
  
  # 5. Join & Calc Distance
  ri_df <- data.frame(
    Group = grp_name, Anchor = anchor_pos, Strand = as.character(strand(gr_to_use)),
    tx_id = mcols(gr_to_use)$tx_id
  )
  merged <- inner_join(ri_df, neighbors_df, by="tx_id")
  
  if (nrow(merged) == 0) return(data.frame())
  
  # 6. Remove Self & Finalize
  merged$key <- paste(merged$tx_id, merged$intron_uid, sep="|")
  merged <- merged[!(merged$key %in% forbidden_df$block_key), ]
  
  merged <- merged %>%
    mutate(raw_dist = intron_center - Anchor) %>%
    mutate(distance = ifelse(Strand == "-", -raw_dist, raw_dist))
  
  return(merged)
}

message("--- Calculating EJC Densities ---")
# 5' EJC Data
ejc_5p <- bind_rows(
  get_ejc_density(dep_ri, "UFM1_dependent", txdb, "5p"),
  get_ejc_density(indep_ri, "UFM1_independent", txdb, "5p"),
  get_ejc_density(control_introns, "Control 3'UTR Intron", txdb, "5p")
) %>% filter(distance >= -WIN_EJC & distance <= WIN_EJC)

# 3' EJC Data
ejc_3p <- bind_rows(
  get_ejc_density(dep_ri, "UFM1_dependent", txdb, "3p"),
  get_ejc_density(indep_ri, "UFM1_independent", txdb, "3p"),
  get_ejc_density(control_introns, "Control 3'UTR Intron", txdb, "3p")
) %>% filter(distance >= -WIN_EJC & distance <= WIN_EJC)


# ==============================================================================
# 4. ANALYSIS 2: STOP CODON DENSITY
# ==============================================================================
analyze_stops_robust <- function(gr, grp_name, fasta_file, window_size=25) {
  
  if (length(gr) == 0) return(NULL)
  
  # A. Fix Chromosome Names
  fasta_levels <- seqlevels(seqinfo(fasta_file))
  if (!any(grepl("^chr", fasta_levels)) && any(grepl("^chr", seqlevels(gr)))) {
    seqlevelsStyle(gr) <- "NCBI"
  } else if (any(grepl("^chr", fasta_levels)) && !any(grepl("^chr", seqlevels(gr)))) {
    seqlevelsStyle(gr) <- "UCSC"
  }
  
  # B. Define Windows (5' and 3')
  is_plus <- as.character(strand(gr)) == "+"
  
  # 5' (Start) and 3' (End)
  c_5p <- ifelse(is_plus, start(gr), end(gr))
  c_3p <- ifelse(is_plus, end(gr), start(gr))
  
  gr_5p <- GRanges(seqnames(gr), IRanges(c_5p - window_size, c_5p + window_size), strand=strand(gr))
  gr_3p <- GRanges(seqnames(gr), IRanges(c_3p - window_size, c_3p + window_size), strand=strand(gr))
  
  # Prune
  full_width <- 2 * window_size + 1
  gr_5p <- gr_5p[width(gr_5p) == full_width]
  gr_3p <- gr_3p[width(gr_3p) == full_width]
  
  if (length(gr_5p) == 0) return(NULL)
  
  # C. Get Seqs & Count
  seqs_5p <- getSeq(fasta_file, gr_5p)
  seqs_3p <- getSeq(fasta_file, gr_3p)
  
  count_motifs <- function(dna_set) {
    counts <- rep(0, full_width)
    for (p in c("TAA", "TAG", "TGA")) {
      m <- vmatchPattern(p, dna_set)
      s <- unlist(startIndex(m))
      if (length(s) > 0) {
        t <- table(s)
        i <- as.numeric(names(t))
        valid <- i <= full_width
        counts[i[valid]] <- counts[i[valid]] + as.numeric(t[valid])
      }
    }
    return(counts / length(dna_set))
  }
  
  # D. Format
  x_axis <- seq(-window_size, window_size)
  df_5p <- data.frame(Position = x_axis, Density = count_motifs(seqs_5p), Site = "5p", Group = grp_name)
  df_3p <- data.frame(Position = x_axis, Density = count_motifs(seqs_3p), Site = "3p", Group = grp_name)
  
  return(rbind(df_5p, df_3p))
}

message("--- Calculating Stop Densities ---")
stop_df <- rbind(
  analyze_stops_robust(dep_ri, "UFM1_dependent", genome_fa, WIN_STOP),
  analyze_stops_robust(indep_ri, "UFM1_independent", genome_fa, WIN_STOP),
  analyze_stops_robust(control_introns, "Control 3'UTR Intron", genome_fa, WIN_STOP)
)

# ==============================================================================
# 5. BASE R PLOTTING
# ==============================================================================
message("--- Generating Combined Plots ---")

# Define Colors & Transparency
make_alpha <- function(col, alpha=50) {
  rgb(t(col2rgb(col)), alpha=alpha, maxColorValue=255)
}
cols_solid <- c("UFM1_dependent"="#E57373", "UFM1_independent"="#64B5F6", "Control 3'UTR Intron"="#9E9E9E")
cols_fill  <- c("UFM1_dependent"=make_alpha("#E57373", 60), "UFM1_independent"=make_alpha("#64B5F6", 60), "Control 3'UTR Intron"=make_alpha("#9E9E9E", 60))

# Open PDF Device
pdf(file.path(opt$outdir, "EJC_PTC_Analysis.pdf"), width=10, height=8)

# Layout: 4 Panels (2x2)
# Row 1: EJC Density (Metagene)
# Row 2: Stop Codon Density (PTC)
par(mfrow=c(2,2), mar=c(4, 4, 3, 1), oma=c(0,0,2,0))

# --- PLOT 1: EJC 5' ANCHOR ---
d_5p <- ejc_5p
d_3p <- ejc_3p

# Calculate dynamic y-max from all EJC densities
get_max_density <- function(data, groups, win) {
  max_vals <- c()
  for(g in groups) {
    vals <- data$distance[data$Group == g]
    if(length(vals) > 10) {
      den <- density(vals, from=-win, to=win, bw=200)
      max_vals <- c(max_vals, max(den$y))
    }
  }
  if(length(max_vals) == 0) return(0.0002)
  return(max(max_vals) * 1.1)
}

groups_list <- c("Control 3'UTR Intron", "UFM1_independent", "UFM1_dependent")
max_y <- max(get_max_density(d_5p, groups_list, WIN_EJC), 
             get_max_density(d_3p, groups_list, WIN_EJC))

plot(1, type="n", xlim=c(-WIN_EJC, WIN_EJC), ylim=c(0, max_y),
     xlab="Distance (bp)", ylab="Splice Junction Density", main="A. Splice Junction Density: 5' Anchored", bty="l", las=1)
rect(0, 0, WIN_EJC, 1, col=rgb(0.9,0.9,0.9,0.5), border=NA) # Shade Intron Body (Right)
abline(v=0, lty=2); text(0, max_y, "Start", pos=3, cex=0.8)

for(g in c("Control 3'UTR Intron", "UFM1_independent", "UFM1_dependent")) {
  vals <- d_5p$distance[d_5p$Group == g]
  if(length(vals)>10) {
    den <- density(vals, from=-WIN_EJC, to=WIN_EJC, bw=200)
    polygon(c(den$x, rev(den$x)), c(den$y, rep(0, length(den$y))), col=cols_fill[g], border=NA)
    lines(den$x, den$y, col=cols_solid[g], lwd=2)
  }
}

# --- PLOT 2: EJC 3' ANCHOR ---
plot(1, type="n", xlim=c(-WIN_EJC, WIN_EJC), ylim=c(0, max_y),
     xlab="Distance (bp)", ylab="Splice Junction Density", main="B. Splice Junction Density: 3' Anchored", bty="l", las=1)
rect(55, 0, WIN_EJC, 1, col=rgb(1,0,0,0.05), border=NA) # Shade NMD Zone
text(2500, max_y*0.8, "NMD Zone", col="darkred", font=3, cex=0.8)
abline(v=0, lty=2); text(0, max_y, "End", pos=3, cex=0.8)

for(g in c("Control 3'UTR Intron", "UFM1_independent", "UFM1_dependent")) {
  vals <- d_3p$distance[d_3p$Group == g]
  if(length(vals)>10) {
    den <- density(vals, from=-WIN_EJC, to=WIN_EJC, bw=200)
    polygon(c(den$x, rev(den$x)), c(den$y, rep(0, length(den$y))), col=cols_fill[g], border=NA)
    lines(den$x, den$y, col=cols_solid[g], lwd=2)
  }
}

# --- PLOT 3: STOP CODON 5' ANCHOR ---
if (!is.null(stop_df) && nrow(stop_df) > 0) {
  s_5p <- stop_df[stop_df$Site == "5p", ]
  max_y_s <- max(s_5p$Density, na.rm=TRUE) * 1.1
  plot(1, type="n", xlim=c(-WIN_STOP, WIN_STOP), ylim=c(0, max_y_s),
       xlab="Position (bp)", ylab="Stop Density", main="C. Stop Codons: 5' Anchored", bty="l", las=1)
  rect(0, 0, WIN_STOP, 1, col=rgb(0.9,0.9,0.9,0.5), border=NA) # Shade Intron
  abline(v=0, lty=2)
  
  for(g in c("Control 3'UTR Intron", "UFM1_independent", "UFM1_dependent")) {
    dat <- s_5p[s_5p$Group == g, ]
    if (nrow(dat) > 0) lines(dat$Position, dat$Density, col=cols_solid[g], lwd=2)
  }
  
  # --- PLOT 4: STOP CODON 3' ANCHOR ---
  s_3p <- stop_df[stop_df$Site == "3p", ]
  plot(1, type="n", xlim=c(-WIN_STOP, WIN_STOP), ylim=c(0, max_y_s),
       xlab="Position (bp)", ylab="Stop Density", main="D. Stop Codons: 3' Anchored", bty="l", las=1)
  rect(-WIN_STOP, 0, 0, 1, col=rgb(0.9,0.9,0.9,0.5), border=NA) # Shade Intron
  abline(v=0, lty=2)
  
  for(g in c("Control 3'UTR Intron", "UFM1_independent", "UFM1_dependent")) {
    dat <- s_3p[s_3p$Group == g, ]
    if (nrow(dat) > 0) lines(dat$Position, dat$Density, col=cols_solid[g], lwd=2)
  }
} else {
  plot.new(); text(0.5, 0.5, "No Stop Codon Data")
  plot.new(); text(0.5, 0.5, "No Stop Codon Data")
}

# --- LEGEND WITH COUNTS ---
# Create labels with counts
leg_labels <- c(
  paste0("UFM1 Dep (n=", counts_info[["UFM1_dependent"]], ")"),
  paste0("UFM1 Indep (n=", counts_info[["UFM1_independent"]], ")"),
  paste0("Control (n=", counts_info[["Control 3'UTR Intron"]], ")")
)

# Reset plotting area to draw legend in center or bottom
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", legend=leg_labels, fill=cols_solid, horiz=TRUE, bty="n", cex=1.0, inset=c(0,0.02))

dev.off()

message("Done. Saved to: ", file.path(opt$outdir, "EJC_PTC_Analysis.pdf"))
