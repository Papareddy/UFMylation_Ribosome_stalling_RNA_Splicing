
# ==============================================================================
# Analyze 3'UTR Properties: Length and Intron Content
# ==============================================================================

suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
    library(ggplot2)
    library(optparse)
})

# ==========================================
# 1. ARGUMENT PARSING
# ==========================================

option_list <- list(
  make_option(c("--gtf"), type="character", default=NULL, 
              help="Path to GTF file", metavar="GTF"),
  make_option(c("--events"), type="character", default=NULL, 
              help="Path to UFM1_events_rich.rds", metavar="Events_RDS"),
  make_option(c("--filter_type"), type="character", default=NULL, 
              help="Filter for specific EventType (e.g. RI, SE)", metavar="TYPE"),
  make_option(c("--outdir"), type="character", default="results/3utr_analysis", 
              help="Output directory", metavar="OUTDIR")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gtf) || is.null(opt$events)) {
  print_help(opt_parser)
  stop("Both --gtf and --events are required.", call.=FALSE)
}

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# ==========================================
# 2. EXTRACT 3'UTRs
# ==========================================

message("Loading GTF and extracting 3'UTRs...")
# Create TxDb from GTF
txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")

# Extract 3'UTRs grouped by transcript
# use.names=TRUE ensures the names of the GRangesList elements are transcript IDs
utrs_by_tx <- threeUTRsByTranscript(txdb, use.names=TRUE)

message("Found ", length(utrs_by_tx), " transcripts with 3'UTRs.")

# ==========================================
# 3. CLASSIFY 3'UTRs (Intron-containing vs Intron-less)
# ==========================================

message("Classifying 3'UTRs...")

# Determine number of exons in each 3'UTR
# A 3'UTR with > 1 exon implies it contains an intron (between those exons)
num_exons <- elementNROWS(utrs_by_tx)

# Calculate total length of 3'UTR (sum of exon widths)
utr_lengths <- sum(width(utrs_by_tx))

# Create classification DF
df <- data.frame(
    tx_id = names(utrs_by_tx),
    num_exons = num_exons,
    length = utr_lengths,
    type = ifelse(num_exons > 1, "3'UTR with Introns", "3'UTR Intron-less"),
    stringsAsFactors = FALSE
)

# ==========================================
# 4. OVERLAP WITH UFM1-DEPENDENT EVENTS
# ==========================================

message("Loading UFM1 events from: ", opt$events)
ufm_events <- readRDS(opt$events)

# Filter by Event Type if requested (e.g. "RI")
if (!is.null(opt$filter_type)) {
    message("Filtering for EventType: ", opt$filter_type)
    ufm_events <- subset(ufm_events, EventType == opt$filter_type)
    if (length(ufm_events) == 0) stop("No events found for Type: ", opt$filter_type)
}

# Split into Dependent and Independent
ufm_dep   <- subset(ufm_events, Group == "UFM1_dependent")
ufm_indep <- subset(ufm_events, Group == "UFM1_independent")

message("Checking overlaps with:")
message("  - UFM1-Dependent: ", length(ufm_dep))
message("  - UFM1-Independent: ", length(ufm_indep))

# Ensure seqlevels match
message("Aligning seqlevels style...")
tryCatch({
    seqlevelsStyle(ufm_dep) <- seqlevelsStyle(utrs_by_tx)[1]
    seqlevelsStyle(ufm_indep) <- seqlevelsStyle(utrs_by_tx)[1]
}, error = function(e) {
    message("Warning: seqlevelsStyle matching failed. Attempting manual 'chr' handling.")
})

if (any(grepl("^chr", seqlevels(ufm_dep))) && !any(grepl("^chr", seqlevels(utrs_by_tx)))) {
    seqlevels(ufm_dep) <- gsub("^chr", "", seqlevels(ufm_dep))
    seqlevels(ufm_indep) <- gsub("^chr", "", seqlevels(ufm_indep))
} else if (!any(grepl("^chr", seqlevels(ufm_dep))) && any(grepl("^chr", seqlevels(utrs_by_tx)))) {
    seqlevels(ufm_dep) <- paste0("chr", seqlevels(ufm_dep))
    seqlevels(ufm_indep) <- paste0("chr", seqlevels(ufm_indep))
}

# Find Overlaps
hits_dep   <- findOverlaps(ufm_dep, utrs_by_tx)
hits_indep <- findOverlaps(ufm_indep, utrs_by_tx)

dep_tx_ids   <- names(utrs_by_tx)[unique(subjectHits(hits_dep))]
indep_tx_ids <- names(utrs_by_tx)[unique(subjectHits(hits_indep))]

# Update Classification
# Priority: Dependent > Independent > Background Introns
df$type[df$tx_id %in% indep_tx_ids] <- "3'UTR with UFM1-Indep Introns"
df$type[df$tx_id %in% dep_tx_ids]   <- "3'UTR with UFM1-Dep Introns"

ordered_levels <- c("All 3'UTRs", "3'UTR Intron-less", "3'UTR with Introns", "3'UTR with UFM1-Indep Introns", "3'UTR with UFM1-Dep Introns")

# ==========================================
# 5. OUTPUTS
# ==========================================

# Save Table
table_file <- file.path(opt$outdir, "3UTR_classification.tsv")
write.table(df, table_file, sep="\t", quote=FALSE, row.names=FALSE)
message("Saved classification table to: ", table_file)

# Summary Stats
message("\nSummary Counts:")
print(table(df$type))

message("\nAverage Lengths:")
print(tapply(df$length, df$type, mean))

# Plotting
message("Generating Boxplot...")

# Create Plotting Data (adding "All" category)
plot_df <- rbind(
    df %>% mutate(category = type),
    df %>% mutate(category = "All 3'UTRs")
)

plot_df$category <- factor(plot_df$category, levels = ordered_levels)

# Colors
cols <- c(
    "All 3'UTRs" = "#E0E0E0",             # Light Grey
    "3'UTR Intron-less" = "#B0BEC5",      # Grey-Blue
    "3'UTR with Introns" = "#90A4AE",     # Darker Grey
    "3'UTR with UFM1-Indep Introns" = "#64B5F6", # Blue
    "3'UTR with UFM1-Dep Introns" = "#E57373"    # Red
)

p <- ggplot(plot_df, aes(x = category, y = length, fill = category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  scale_y_log10(labels = scales::comma) + 
  scale_fill_manual(values = cols) +
  labs(
    title = paste0("Length Distribution of 3'UTRs (", ifelse(is.null(opt$filter_type), "All Events", opt$filter_type), ")"),
    subtitle = "Comparison of UFM1-Associated vs Genome-Wide",
    x = "Category",
    y = "Length (bp, log10 scale)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_file <- file.path(opt$outdir, "3UTR_length_boxplot.pdf")
ggsave(plot_file, p, width = 6, height = 6)
message("Saved boxplot to: ", plot_file)
