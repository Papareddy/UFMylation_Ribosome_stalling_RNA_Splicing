
# ==============================================================================
# Analyze Intron Lengths: UFM1 RI vs Constitutive
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
  make_option(c("--outdir"), type="character", default="results/intron_analysis", 
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
# 2. LOAD DATA
# ==========================================

message("Loading UFM1 events from: ", opt$events)
ufm_events <- readRDS(opt$events)

# Filter for RI Only
ufm_ri <- subset(ufm_events, EventType == "RI")
message("Total RI Events Loaded: ", length(ufm_ri))

ufm_dep <- subset(ufm_ri, Group == "UFM1_dependent")
ufm_indep <- subset(ufm_ri, Group == "UFM1_independent")

message("  - UFM1-Dependent RI: ", length(ufm_dep))
message("  - UFM1-Independent RI: ", length(ufm_indep))


message("\nExtracting all introns from GTF (Constitutive background)...")
txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")
all_introns <- unlist(intronsByTranscript(txdb))

message("Total Introns in Genome: ", length(all_introns))


# ==========================================
# 3. OVERLAP & FILTER
# ==========================================

message("\nAligning seqlevels and filtering constitutive introns...")

# Alignment
tryCatch({
    seqlevelsStyle(ufm_dep) <- seqlevelsStyle(all_introns)[1]
    seqlevelsStyle(ufm_indep) <- seqlevelsStyle(all_introns)[1]
}, error = function(e) {
    message("Warning: seqlevelsStyle matching failed. Attempting manual 'chr' handling.")
})
# Manual fallback
if (any(grepl("^chr", seqlevels(ufm_dep))) && !any(grepl("^chr", seqlevels(all_introns)))) {
    seqlevels(ufm_dep) <- gsub("^chr", "", seqlevels(ufm_dep))
    seqlevels(ufm_indep) <- gsub("^chr", "", seqlevels(ufm_indep))
} else if (!any(grepl("^chr", seqlevels(ufm_dep))) && any(grepl("^chr", seqlevels(all_introns)))) {
    seqlevels(ufm_dep) <- paste0("chr", seqlevels(ufm_dep))
    seqlevels(ufm_indep) <- paste0("chr", seqlevels(ufm_indep))
}

# Identify Overlaps
# Any intron that overlaps with a UFM1 RI event should NOT be in the constitutive set
# We use a broad overlap check
overlaps_dep <- countOverlaps(all_introns, ufm_dep)
overlaps_indep <- countOverlaps(all_introns, ufm_indep)

# Define Constitutive: Introns that have ZERO overlap with our RI sets
const_introns <- all_introns[overlaps_dep == 0 & overlaps_indep == 0]

message("Constitutive Introns (Filtered): ", length(const_introns))

# ==========================================
# 4. PREPARE DATA
# ==========================================

df <- rbind(
    data.frame(length = width(ufm_dep), type = "UFM1-Dependent RI"),
    data.frame(length = width(ufm_indep), type = "UFM1-Independent RI"),
    data.frame(length = width(const_introns), type = "Constitutive Introns")
)

df$type <- factor(df$type, levels = c("Constitutive Introns", "UFM1-Independent RI", "UFM1-Dependent RI"))

# ==========================================
# 5. OUTPUTS
# ==========================================

# Save Table
table_file <- file.path(opt$outdir, "RI_length_stats.tsv")
stats <- df %>% group_by(type) %>% summarise(
    Count = n(),
    Mean_Length = mean(length),
    Median_Length = median(length),
    Min_Length = min(length),
    Max_Length = max(length)
)
write.table(stats, table_file, sep="\t", quote=FALSE, row.names=FALSE)
message("Saved stats table to: ", table_file)
print(stats)

# Plotting
message("Generating Boxplot...")
cols <- c(
    "Constitutive Introns" = "#B0BEC5",    # Grey
    "UFM1-Independent RI" = "#64B5F6",     # Blue
    "UFM1-Dependent RI" = "#E57373"        # Red
)

p <- ggplot(df, aes(x = type, y = length, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  scale_y_log10(labels = scales::comma) + 
  scale_fill_manual(values = cols) +
  labs(
    title = "Intron Length Comparison",
    subtitle = "UFM1-Associated RI vs Constitutive Introns",
    x = "Category",
    y = "Length (bp, log10 scale)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_file <- file.path(opt$outdir, "RI_length_comparison.pdf")
ggsave(plot_file, p, width = 6, height = 6)
message("Saved boxplot to: ", plot_file)
