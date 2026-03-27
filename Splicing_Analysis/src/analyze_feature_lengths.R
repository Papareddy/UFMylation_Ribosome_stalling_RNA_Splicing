#!/usr/bin/env Rscript

# analyze_feature_lengths.R
# Comprehensive Length Analysis for Introns and 3'UTRs
# Categories: Constitutive, 3'UTR-Introns, UFM1-Dependent, UFM1-Independent

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

# Ensure consistent chromosome naming style (e.g. "chr1" vs "1")
if (!any(seqlevels(events) %in% seqlevels(txdb))) {
  message("Chromosome naming styles differ. Attempting to harmonize...")
  seqlevelsStyle(events) <- seqlevelsStyle(txdb)[1]
}

# Filter UFM1 Groups
target_events <- events # Generalize to all types
# Define Groups from Events
dep_events <- target_events[target_events$Group == "UFM1_dependent"]
indep_events <- target_events[target_events$Group == "UFM1_independent"]

# ==============================================================================
# PART A: INTRON LENGTHS
# ==============================================================================
message("Analyzing Intron Lengths...")

# 1. All Constitutive Introns (Background)
# We take all introns from TxDb. 
# "Constitutive" implies strictly present in all isoforms, which is hard to determine from GTF alone easily.
# We will approximate "All Introns" or "Background Introns".
# To be cleaner, we take unique introns defined in the GTF.
all_introns_grl <- intronsByTranscript(txdb, use.names=TRUE)
all_introns <- unlist(all_introns_grl)
all_introns_unique <- unique(all_introns)

# 2. Introns in 3' UTRs
utr3_grl <- threeUTRsByTranscript(txdb, use.names=TRUE)
utr3_all <- unlist(utr3_grl)

# Find overlap: Introns completely within 3'UTRs
# (Intron retention in 3'UTR is common). 
# Wait, if it's an intron *in* a 3'UTR, it's spliced out in that transcript context? 
# Or are we talking about the *region* that is an intron in some transcripts but UTR in others?
# User likely means: Introns located within the 3'UTR region of a gene model.
# Typically, standard gene models don't have "introns in 3'UTRs" unless it's a multi-exon 3'UTR.
# Yes, multi-exon 3'UTRs exist.
utr3_introns_hits <- findOverlaps(all_introns_unique, utr3_all, type="within")
utr3_introns <- all_introns_unique[unique(queryHits(utr3_introns_hits))]

# Prepare Data Frame
df_introns_list <- list()

if (length(dep_events) > 0) {
  df_introns_list[[length(df_introns_list) + 1]] <- data.frame(Length = width(dep_events), Category = "UFM1_dependent_event")
}
if (length(indep_events) > 0) {
  df_introns_list[[length(df_introns_list) + 1]] <- data.frame(Length = width(indep_events), Category = "UFM1_independent_event")
}
if (length(utr3_introns) > 0) {
  df_introns_list[[length(df_introns_list) + 1]] <- data.frame(Length = width(utr3_introns), Category = "3'UTR Introns")
}
if (length(all_introns_unique) > 0) {
  df_introns_list[[length(df_introns_list) + 1]] <- data.frame(Length = width(all_introns_unique), Category = "All Constitutive Introns")
}

df_introns <- bind_rows(df_introns_list)

# Plotting A
message("Plotting Intron Lengths...")
p_a <- ggplot(df_introns, aes(x = Category, y = Length, fill = Category)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(trans = "log2", 
                     breaks = c(100, 500, 1000, 5000, 10000, 50000),
                     labels = scales::comma) +
  theme_classic() +
  labs(title = "Event / Intron Length Distribution", y = "Length (bp, log2)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + 
  stat_compare_means(ref.group = ".all.", label = "p.signif", method = "wilcox.test")

ggsave(file.path(opt$outdir, "Length_Analysis_Introns.pdf"), p_a, width = 6, height = 6)
write_tsv(df_introns, file.path(opt$outdir, "Length_Data_Introns.tsv"))


# ==============================================================================
# PART B: 3' UTR LENGTHS
# ==============================================================================
message("Analyzing 3' UTR Lengths...")

# 1. All 3' UTRs (Aggregated Length per Transcript)
# Use transcriptLengths for accuracy
tx_lens <- transcriptLengths(txdb, with.utr3_len=TRUE)
utr3_lens_vec <- setNames(tx_lens$utr3_len, tx_lens$tx_name)
utr3_lens_vec <- utr3_lens_vec[utr3_lens_vec > 0]

# 2. 3' UTRs with Introns (Multi-exon 3'UTRs)
# We use GRangesList to check exon count (>= 2 exons means there's an intron)
utr3_grl <- threeUTRsByTranscript(txdb, use.names=TRUE)
is_multi_exon <- elementNROWS(utr3_grl) > 1
tx_multi_exon <- names(utr3_grl)[is_multi_exon]
tx_single_exon <- names(utr3_grl)[!is_multi_exon]

utr3_with_introns <- utr3_lens_vec[tx_multi_exon]
utr3_no_introns <- utr3_lens_vec[tx_single_exon]
utr3_with_introns <- na.omit(utr3_with_introns)
utr3_no_introns <- na.omit(utr3_no_introns)

# 3. 3' UTRs with UFM1-dependent/independent Introns
map_event_to_utr_len <- function(events_subset, label) {
  hits <- findOverlaps(events_subset, utr3_grl)
  
  # For each hit, get Transcript ID (subject) and Event Width (query)
  tx_ids <- names(utr3_grl)[subjectHits(hits)]
  event_widths <- width(events_subset)[queryHits(hits)]
  
  # Get spliced UTR length for these transcripts
  spliced_lens <- utr3_lens_vec[tx_ids]
  
  res <- data.frame(
    tx_id = tx_ids,
    event_width = event_widths,
    spliced_utr_len = as.numeric(spliced_lens)
  )
  res <- na.omit(res)
  
  # Effective Length with Retention = Spliced UTR + Intron Length
  res$effective_len <- res$spliced_utr_len + res$event_width
  
  if(nrow(res) == 0) return(data.frame(Length = numeric(0), Category = character(0)))
  
  return(data.frame(Length = res$effective_len, Category = label))
}

dep_utr_lens <- map_event_to_utr_len(dep_events, "3'UTRs with UFM1-dependent")
indep_utr_lens <- map_event_to_utr_len(indep_events, "3'UTRs with UFM1-independent")

# Backgrounds
bg_all_utr <- data.frame(Length = as.numeric(utr3_lens_vec), Category = "All 3'UTRs")
bg_with_int <- data.frame(Length = as.numeric(utr3_with_introns), Category = "3'UTRs with Introns")
bg_no_int <- data.frame(Length = as.numeric(utr3_no_introns), Category = "3'UTRs without Introns")

df_utr <- bind_rows(
  bg_all_utr,
  bg_with_int,
  bg_no_int,
  dep_utr_lens,
  indep_utr_lens
)

# Plotting B
message("Plotting 3'UTR Lengths...")
p_b <- ggplot(df_utr, aes(x = Category, y = Length, fill = Category)) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(trans = "log2",
                     breaks = c(100, 500, 1000, 5000, 10000, 50000),
                     labels = scales::comma) +
  theme_classic() +
  labs(title = "3'UTR Length Distribution", y = "Length (bp, log2)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") + 
  stat_compare_means(ref.group = ".all.", label = "p.signif", method = "wilcox.test")

ggsave(file.path(opt$outdir, "Length_Analysis_3UTR.pdf"), p_b, width = 6, height = 6)
write_tsv(df_utr, file.path(opt$outdir, "Length_Data_3UTR.tsv"))

message("Done.")
