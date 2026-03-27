#!/usr/bin/env Rscript

# plot_nmd_metagene.R
#
# Purpose: Visualize the "NMD Risk Landscape" / Architectural Topology.
# Method: Metagene Splice-Junction Density Plot.
#   - Anchor: Center of Retained Intron (0).
#   - Signal: Location of neighbor splice junctions (intron centers).
#   - Hypothesis: High density > 50bp downstream = NMD Trigger.

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
  library(optparse)
})

option_list <- list(
  make_option(c("--gtf"), type="character", help="Path to GTF file"),
  make_option(c("--events"), type="character", help="Path to UFM1 events RDS file"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory"),
  make_option(c("--event_type"), type="character", default="RI", help="Event type to process (e.g. SE, A3SS, RI)")
)

opt <- parse_args(OptionParser(option_list=option_list))
if(!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

message("Loading Data...")
txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")
events <- readRDS(opt$events)

# Filter UFM1 Groups
target_events <- events # Generalize to all types
# Use pre-calculated Group column
dep_events <- target_events[target_events$Group == "UFM1_dependent"]
indep_events <- target_events[target_events$Group == "UFM1_independent"]

# Assign unique IDs for tracking
if (length(dep_events) > 0) {
  mcols(dep_events)$Event_ID <- paste0("Dep_", seq_along(dep_events))
}
if (length(indep_events) > 0) {
  mcols(indep_events)$Event_ID <- paste0("Indep_", seq_along(indep_events))
}

combined_events <- c(dep_events, indep_events)

if (length(combined_events) > 0) {
  # Only keep necessary cols
  mcols(combined_events) <- mcols(combined_events)[, c("Group", "Event_ID")]
}

# --- 1b. Prepare Control Group ---
if (opt$event_type == "RI") {
  message("Extracting Constitutive 3'UTR Introns (Control for RI)...")
  control_group_label <- "Control 3'UTR Intron"
  
  # 1. Get All Introns (with tx_id names)
  all_features_grl <- intronsByTranscript(txdb, use.names=TRUE)
  all_features_flat <- unlist(all_features_grl, use.names=TRUE)
  
  # 2. Get 3'UTR Ranges (Span from start to end of UTR)
  utr3_grl <- threeUTRsByTranscript(txdb, use.names=TRUE)
  utr3_ranges <- unlist(range(utr3_grl), use.names=TRUE)
  
  # 3. Filter: Keep introns that are fully contained in a 3'UTR
  control_features_flat <- subsetByOverlaps(all_features_flat, utr3_ranges, type="within")
} else {
  message("Extracting Constitutive Exons (Control for SE/A3SS/A5SS/MXE)...")
  control_group_label <- "Control Constitutive Exon"
  
  # 1. Get Exons
  all_features_grl <- exonsBy(txdb, by="tx", use.names=TRUE)
  all_features_flat <- unlist(all_features_grl, use.names=TRUE)
  
  # 2. Heuristic for "Constitutive": exons that appear in many transcripts of a gene, 
  # but for metagene just picking random exons or internal exons works.
  # Assuming any exon is a valid control for metagene.
  control_features_flat <- all_features_flat
}

message(paste("Found", length(control_features_flat), "Control features."))

# Make a DF
control_ri_df <- data.frame(
  RI_ID = paste0("Control_", seq_along(control_features_flat)),
  Group = control_group_label,
  RI_Center = (start(control_features_flat) + end(control_features_flat)) / 2,
  Strand = as.character(strand(control_features_flat)),
  tx_id = names(control_features_flat) 
)

# Combine RI Maps
# We need to adapt the RI logic to match this pre-calculated map
# Current logic maps RIs to txs. Control RIs are ALREADY mapped to txs.

# --- 1. Prepare Transcript Architecture ---
message("Preparing Transcript Architecture...")

# Get All Introns flattened (These are the potential EJC sites)
introns_grl <- intronsByTranscript(txdb, use.names=TRUE)
# Flatten to DF for joining
introns_flat <- unlist(introns_grl, use.names=TRUE)
introns_df <- data.frame(
  tx_id = names(introns_flat),
  intron_center = (start(introns_flat) + end(introns_flat)) / 2,
  width = width(introns_flat)
)

# Get Transcripts (to map RI -> TX)
# We map RI to 'transcripts' to find the Host TX ID
transcripts_gr <- transcripts(txdb)
# Align styles
tryCatch({seqlevelsStyle(combined_ri) <- seqlevelsStyle(transcripts_gr)[1]}, error=function(e){})

# --- 2. Map Events to Host Transcripts ---
message("Mapping Events to Host Transcripts...")

# Find overlap
event_tx_map <- data.frame()

if (length(combined_events) > 0) {
    hits <- findOverlaps(combined_events, transcripts_gr)
    # We create a mapping table: Event_ID -> TX_ID
    # Create Event Map (Experimental Groups)
    event_tx_map <- data.frame(
      Event_Index = queryHits(hits),
      TX_Index = subjectHits(hits)
    )
    if (nrow(event_tx_map) > 0) {
        event_tx_map$Event_ID <- combined_events$Event_ID[event_tx_map$Event_Index]
        event_tx_map$Event_Group <- combined_events$Group[event_tx_map$Event_Index]
        event_tx_map$Event_Center <- (start(combined_events)[event_tx_map$Event_Index] + end(combined_events)[event_tx_map$Event_Index]) / 2
        event_tx_map$Strand <- as.character(strand(combined_events)[event_tx_map$Event_Index])
        event_tx_map$tx_id <- transcripts_gr$tx_name[event_tx_map$TX_Index]
    }
}

# Harmonize Control DF columns
control_event_df <- control_ri_df %>% select(Event_ID=RI_ID, Event_Group=Group, Event_Center=RI_Center, Strand, tx_id)

# Combine Experimental and Control Maps
event_tx_map <- bind_rows(event_tx_map, control_event_df)

# !!! CRITICAL FILTER !!!
# An RI might overlap multiple transcripts (isoforms).
# Ideally we want the isoform where this RI *matches* a known intron or fits the structure.
# But for "Metagene" signal, averaging across valid overlaps is acceptable 
# OR limiting to those where the RI overlaps a known intron in that transcript.
# For simplicity and robustness (showing POTENTIAL risk), we keep valid transcript overlaps.

# --- 3. Calculate Junction Distances (Vectorized) ---
message("Calculating Junction Distances...")

# Join Event Map with All Introns on TX_ID
# Inner join: We only care about Events mapped to transcripts that HAVE introns.
merged_df <- inner_join(event_tx_map, introns_df, by="tx_id")

# Calculate Distance
merged_df$raw_dist <- merged_df$intron_center - merged_df$Event_Center

# Strand Correction
# If - strand: Downstream is smaller coordinate. 
# Event (1000) -> Downstream Intron (500). Dist = 500 - 1000 = -500.
# We want this to be POSITIVE distance (Downstream).
# So for - strand, distinct = -raw_dist.
merged_df$distance <- ifelse(merged_df$Strand == "-", -merged_df$raw_dist, merged_df$raw_dist)

# Filter
# 1. Remove "Self" (Distance approx 0, i.e., the RI itself)
# Use a threshold, e.g., < 10 bp difference in centers
filtered_df <- merged_df %>% 
  filter(abs(distance) > 10) %>%
  filter(distance >= -2000 & distance <= 2000)

message(paste("Junctions plotted:", nrow(filtered_df)))

# --- 4. Plotting ---
message("Generating Plot...")

# Define Colors
cols <- c("UFM1_dependent" = "#E57373", "UFM1_independent" = "#64B5F6", 
          "Control 3'UTR Intron" = "#9E9E9E", "Control Constitutive Exon" = "#9E9E9E")

if (nrow(filtered_df) > 0) {
    p <- ggplot(filtered_df, aes(x = distance, color = Event_Group, fill = Event_Group)) +
      # Density Plot
      geom_density(alpha = 0.2, linewidth = 1) +
      
      # The "Danger Zone" Annotation for NMD
      annotate("rect", xmin = 55, xmax = 2000, ymin = 0, ymax = Inf, 
               alpha = 0.05, fill = "red") +
      annotate("text", x = 1000, y = 0, label = "Downstream EJC\n(NMD Trigger Zone)", 
               vjust = -0.5, color = "darkred", fontface = "italic") +
      
      # Vertical Line for RI Center
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      
      # Aesthetics
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols) +
      theme_classic(base_size = 14) +
      labs(
        title = "Metagene Profile: NMD Signature Density",
        subtitle = "Relative position of neighboring splice junctions",
        x = "Distance from Splicing Event Center (bp)\n< Upstream | Downstream >",
        y = "Density of Splice Junctions"
      ) +
      scale_x_continuous(breaks = seq(-2000, 2000, 500)) +
      theme(legend.position = "top")

    outfile <- file.path(opt$outdir, "NMD_Metagene_Profile.pdf")
    ggsave(outfile, p, width = 8, height = 6)
    message("Saved plot to: ", outfile)
} else {
    message("No junctions to plot.")
}
message("Done.")
