#!/usr/bin/env Rscript

# analyze_nmd_features.R
# Optimized Version (Vectorized)

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
  library(optparse)
  library(ggpubr)
})

# --- Argument Parsing ---
option_list <- list(
  make_option(c("--gtf"), type="character", help="Path to GTF file"),
  make_option(c("--events"), type="character", help="Path to UFM1 events RDS file"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))
if(!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

message("Loading Data...")
txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")
events <- readRDS(opt$events)

# Filter UFM1 Groups
target_events <- events # Generalize to all types
# Classify Dependent vs Independent using pre-calculated Group column
dep_events <- target_events[target_events$Group == "UFM1_dependent"]
indep_events <- target_events[target_events$Group == "UFM1_independent"]

# Assign ID to mcols for tracking
if (length(dep_events) > 0) {
  mcols(dep_events)$Event_ID <- paste0("Dep_", seq_along(dep_events))
}
if (length(indep_events) > 0) {
  mcols(indep_events)$Event_ID <- paste0("Indep_", seq_along(indep_events))
}

combined_events <- c(dep_events, indep_events)

if (length(combined_events) > 0) {
  # Only keep necessary mcols for speed
  mcols(combined_events) <- mcols(combined_events)[, c("Group", "Event_ID")]
}

# --- 1. Prepare Features (Vectorized) ---
message("Preparing Transcript Features (Vectorized)...")

# Get All Introns flattened with Rank info
introns_grl <- intronsByTranscript(txdb, use.names=TRUE)
introns_flat <- unlist(introns_grl, use.names=TRUE)
# Assign metadata to flat introns
introns_flat$tx_id <- names(introns_flat)
# Calculate ranks
# unlist() preserves order of GRangesList
introns_flat$rank <- sequence(elementNROWS(introns_grl))
introns_flat$total <- rep(elementNROWS(introns_grl), elementNROWS(introns_grl))
introns_flat$is_last_intron <- (introns_flat$rank == introns_flat$total)

# Get 3'UTR Lengths for Transcripts
utr3_grl <- threeUTRsByTranscript(txdb, use.names=TRUE)
utr3_lens <- sum(width(utr3_grl))
utr3_df <- data.frame(tx_id=names(utr3_lens), spliced_utr_len=as.numeric(utr3_lens))

# --- 2. Overlap Analysis ---
message("Overlapping Events with Transcript Features...")

# Align Styles
tryCatch({seqlevelsStyle(combined_events) <- seqlevelsStyle(introns_flat)[1]}, error=function(e){})
tryCatch({seqlevelsStyle(combined_events) <- seqlevelsStyle(utr3_grl)[1]}, error=function(e){})

# A. Find corresponding Intron (EJC Rule)
# We map event to the Intron it corresponds to (conceptually equivalent for RI, for SE it might overlap)
hits_intron <- findOverlaps(combined_events, introns_flat)

results_df <- as.data.frame(hits_intron)
results_df$Event_ID <- combined_events$Event_ID[results_df$queryHits]
results_df$Group <- combined_events$Group[results_df$queryHits]
results_df$tx_id <- introns_flat$tx_id[results_df$subjectHits]
results_df$is_last_intron <- introns_flat$is_last_intron[results_df$subjectHits]

# Define NMD Status
results_df$NMD_Status <- ifelse(results_df$is_last_intron, "Escape", "Risk")

if (length(combined_events) > 0) {
    # B. Add Length Info (Long 3'UTR Rule)
    # Filter for matches where the transcript actually HAS a 3'UTR
    results_df <- results_df %>% filter(tx_id %in% names(utr3_grl))

    # Add Spliced Length
    results_df <- left_join(results_df, utr3_df, by="tx_id")

    # Add Event Length
    event_width_df <- data.frame(Event_ID = combined_events$Event_ID, Event_Len = width(combined_events))
    results_df <- left_join(results_df, event_width_df, by="Event_ID")

    # Calculate Isoform Length
    results_df$Isoform_UTR_Len <- results_df$spliced_utr_len + results_df$Event_Len

    # --- 3. Summary & Visualization ---
    message("Summarizing and Plotting...")

    # Aggregate to Unique RI (taking consensus or 1st match if multiple transcripts?)
    # RIs often overlap multiple transcripts.
    # For NMD Risk: If it escapes in ONE isoform, does it save the event? Or do we average?
    # Let's count *Transcript-RI pairs* as the unit of analysis, 
    # or pick the "Primary" transcript.
    # Let's stick to Transcript-RI pairs for robust stats.

    final_res <- results_df

    # Save Table
    write_tsv(final_res, file.path(opt$outdir, "NMD_analysis_details.tsv"))

    # Barplot
    plot_data_ejc <- final_res %>% 
      group_by(Group, NMD_Status) %>%
      summarise(Count = n()) %>%
      mutate(Percentage = Count / sum(Count) * 100)

    if (nrow(plot_data_ejc) > 0) {
        p1 <- ggplot(plot_data_ejc, aes(x=Group, y=Percentage, fill=NMD_Status)) +
          geom_bar(stat="identity", position="stack") +
          scale_fill_manual(values=c("Escape"="#66BB6A", "Risk"="#EF5350")) + 
          labs(title="NMD Risk (50-nt Rule / Downstream Introns)",
               subtitle="Escape = Event matches/overlaps Last Intron\nRisk = Downstream introns exist",
               y="Percentage of Transcripts") +
          theme_classic()

        ggsave(file.path(opt$outdir, "NMD_risk_barplot.pdf"), p1, width=6, height=5)
    }

    # Length Boxplot
    if (nrow(final_res) > 0) {
        p2 <- ggplot(final_res, aes(x=Group, y=Isoform_UTR_Len, fill=Group)) +
          geom_boxplot(outlier.shape=NA) +
          coord_cartesian(ylim=c(0, 10000)) +
          geom_hline(yintercept=2000, linetype="dashed", color="red") + 
          labs(title="Full 3'UTR Isoform Lengths",
               subtitle="Red Line: >2000bp (Long 3'UTR NMD Threshold)",
               y="Length (bp)") +
          theme_classic() +
          stat_compare_means(ref.group = "UFM1_dependent", label = "p.signif", method = "wilcox.test")

        ggsave(file.path(opt$outdir, "NMD_length_boxplot.pdf"), p2, width=5, height=6)
    }

    message("Summary of NMD Risk:")
    print(plot_data_ejc)
} else {
    message("No events to analyze.")
    # Create empty files to avoid breaking pipeline
    file.create(file.path(opt$outdir, "NMD_analysis_details.tsv"))
}
message("Done.")
