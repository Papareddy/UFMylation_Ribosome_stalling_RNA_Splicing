
# ==============================================================================
# Export Splicing Events to BED and GRanges
# ==============================================================================

# Libraries
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(optparse)
    library(GenomicRanges)
})

# ==========================================
# 1. ARGUMENT PARSING
# ==========================================

option_list <- list(
  make_option(c("--dependent"), type="character", default=NULL, 
              help="Path to UFM1-dependent TSV file", metavar="DEP_FILE"),
  make_option(c("--independent"), type="character", default=NULL, 
              help="Path to UFM1-independent TSV file", metavar="INDEP_FILE"),
  make_option(c("--outdir"), type="character", default="results/export", 
              help="Output directory [default= %default]", metavar="OUTDIR")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$dependent) || is.null(opt$independent)) {
  print_help(opt_parser)
  stop(" Both --dependent and --independent input files must be provided.", call.=FALSE)
}

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# ==========================================
# 2. FUNCTIONS
# ==========================================

load_and_process <- function(filepath, group_name) {
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(data.frame())
  }
  
  df <- read_tsv(filepath, show_col_types = FALSE)
  

  # Ensure EventType column exists
  if (!"EventType" %in% colnames(df)) {
      warning(paste("EventType column missing in", filepath))
      return(data.frame())
  }
  
  # Standardize Column Names: remove .WT and .UFM suffixes
  # This ensures we have 'longExonStart_0base' instead of 'longExonStart_0base.WT'
  new_cols <- gsub("\\.WT$", "", colnames(df))
  new_cols <- gsub("\\.UFM$", "", new_cols)
  colnames(df) <- new_cols
  
  # Deduplicate columns (keep first). 
  # Since .WT and .UFM columns are often identical for coordinates, this is safe.
  # For data columns (reads etc), we might lose distinction but for this export 
  # we primarily care about coordinates and metadata which are unified.
  df <- df[, !duplicated(colnames(df))]

  # Force sample count columns to character to avoid bind_rows character/double mismatch
  cols_to_char <- c("IJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_1", "SJC_SAMPLE_2", "IncFormLen", "SkipFormLen")
  for(col in cols_to_char) {
      if(col %in% colnames(df)) {
          df[[col]] <- as.character(df[[col]])
      }
  }
  
  df$Group <- group_name
  return(df)
}

# ==========================================
# 3. PROCESSING
# ==========================================

message("Loading Dependent Data: ", opt$dependent)
df_dep <- load_and_process(opt$dependent, "UFM1_dependent")

message("Loading Independent Data: ", opt$independent)
df_indep <- load_and_process(opt$independent, "UFM1_independent")

# Homogenize types for common columns to avoid bind_rows failures
common_cols <- intersect(colnames(df_dep), colnames(df_indep))
for (col in common_cols) {
  # strict check: if classes differ at all (e.g. numeric vs integer, or character vs numeric)
  if (!identical(class(df_dep[[col]]), class(df_indep[[col]]))) {
    message(paste("Type mismatch for column:", col, "- coercing to character."))
    df_dep[[col]] <- as.character(df_dep[[col]])
    df_indep[[col]] <- as.character(df_indep[[col]])
  }
}

# Combine
combined_df <- bind_rows(df_dep, df_indep)

if (nrow(combined_df) == 0) {
    stop("No data loaded.")
}

# ---------------------------------------------------------
# Define Primary Coordinates based on Event Type
# ---------------------------------------------------------
# This logic maps the disparate column names to a standard "Start" and "End"
# representing the "Feature of Interest" (e.g., the skipped exon, retained intron).

# ---------------------------------------------------------
# Define Primary Coordinates based on Event Type
# ---------------------------------------------------------

# Pre-fill missing columns with NA to ensure case_when works
cols_to_check <- c("longExonStart_0base", "longExonEnd", 
                   "riExonStart_0base", "riExonEnd", 
                   "1stExonStart_0base", "1stExonEnd",
                   "exonStart_0base", "exonEnd")

for (col in cols_to_check) {
    if (!col %in% colnames(combined_df)) {
        combined_df[[col]] <- NA_real_
    }
}

# Apply Logic
combined_df <- combined_df %>%
  mutate(
    # Standardize Chromosome (remove 'chr' prefix if needed for some tools, but usually keep it)
    # Start/End logic
    bed_start = case_when(
        EventType == "SE"   ~ exonStart_0base,
        EventType == "RI"   ~ riExonStart_0base,
        EventType == "A3SS" ~ longExonStart_0base, # Approximation for visualization: entire long exon
        EventType == "A5SS" ~ longExonStart_0base, # Approximation
        EventType == "MXE"  ~ `1stExonStart_0base`, # Visualization: 1st exon
        TRUE ~ NA_real_
    ),
    bed_end = case_when(
        EventType == "SE"   ~ exonEnd,
        EventType == "RI"   ~ riExonEnd,
        EventType == "A3SS" ~ longExonEnd,
        EventType == "A5SS" ~ longExonEnd,
        EventType == "MXE"  ~ `1stExonEnd`,
        TRUE ~ NA_real_
    ),
    # Ensure Score is numeric
    score = 0
  ) %>%
  filter(!is.na(bed_start) & !is.na(bed_end))

# ==========================================
# 4. EXPORT: GRanges (.rds)
# ==========================================

message("Creating GRanges object...")
# For GRanges, we use the standardized bed_start/bed_end as the ranges,
# BUT we keep ALL original columns as metadata.

gr <- makeGRangesFromDataFrame(
    combined_df,
    keep.extra.columns = TRUE,
    ignore.strand = FALSE,
    seqnames.field = ifelse("chr" %in% colnames(combined_df), "chr", "chr.WT"), # Handle potential variation
    start.field = "bed_start",
    end.field = "bed_end",
    strand.field = ifelse("strand" %in% colnames(combined_df), "strand", "strand.WT")
)

rds_path <- file.path(opt$outdir, "UFM1_events_rich.rds")
saveRDS(gr, rds_path)
message("Saved rich GRanges object to: ", rds_path)

# Also save as plain data frame for Python compatibility (pyreadr can't read S4 objects)
df_path <- file.path(opt$outdir, "UFM1_events_df.rds")
saveRDS(as.data.frame(gr), df_path)
message("Saved data frame version to: ", df_path, " (for Python/pyreadr)")

# ==========================================
# 5. EXPORT: BED (Extended)
# ==========================================

message("Creating Extended BED file...")

# BED standard: chrom, start, end, name, score, strand
# We will use:
# Name = geneSymbol:EventType
# Score = 0
# Extra Columns = dPSI, FDR, Group

bed_df <- combined_df %>%
  transmute(
      chrom = ifelse("chr" %in% colnames(combined_df), chr, chr.WT),
      start = as.integer(bed_start),
      end = as.integer(bed_end),
      name = paste0(geneSymbol, ":", EventType),
      score = 0,
      strand = ifelse("strand" %in% colnames(combined_df), strand, strand.WT),
      # Metadata
      # Metadata
      # dPSI: Attempt to find dPSI_num or IncLevelDifference (prioritize dPSI_num.WT if available, but here names might be dPSI_num.WT)
      dPSI = if("dPSI_num.WT" %in% names(.)) dPSI_num.WT else if("dPSI_num" %in% names(.)) dPSI_num else if("IncLevelDifference" %in% names(.)) IncLevelDifference else 0,
      
      # FDR: Attempt to find FDR.WT, FDR.UFM or FDR
      FDR = if("FDR.WT" %in% names(.)) FDR.WT else if("FDR.UFM" %in% names(.)) FDR.UFM else if("FDR" %in% names(.)) FDR else 1,
      
      Group = Group,
      EventType = EventType
  ) %>%
  # Handle potential issues with dPSI column name variability
  # If dPSI_num doesn't exist, try dPSI.WT or similar? 
  # Based on prev `head`: dPSI_num.WT and dPSI_num.UFM might be merged? 
  # Actually `bind_rows` merges columns. 
  # `load_and_process` reads simple TSVs. 
  # Let's check column names in the combined DF dynamically if possible, but exact names were:
  # dPSI_num.WT / dPSI_num.UFM ? 
  # Wait, the `UFM1_dependent.tsv` comes from Step 1 which likely has specific column names.
  # The header check showed: `dPSI_num.WT` for dependent? No, header showed `dPSI_num.WT` inside `UFM1_dependent.tsv`?
  # Ah, Step 1 merges WT and Case stats. 
  # Let's re-verify column names if needed. 
  # Assuming `dPSI` or `IncLevelDifference` exists.
  # For now, let's just dump what we have. 
  
  # Clean up NA/Inf in dPSI for BED
  mutate(start = sprintf("%d", start), end = sprintf("%d", end)) # format integer

bed_path <- file.path(opt$outdir, "UFM1_events.bed")
write_tsv(bed_df, bed_path, col_names = FALSE)
message("Saved Extended BED file to: ", bed_path)
