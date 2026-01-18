# ==========================================
# 1. SETUP & LIBRARIES
# ==========================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# Standard Plot Settings
par(mfrow=c(4,4), las=2, tcl=-0.3, bty="n")

# -----------------------------
# ARGUMENT PARSING (Custom Implementation)
# -----------------------------
get_cli_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  arg_list <- list()
  for (arg in args) {
    if (startsWith(arg, "--")) {
      parts <- strsplit(sub("^--", "", arg), "=")[[1]]
      if (length(parts) == 2) {
        arg_list[[parts[1]]] <- parts[2]
      }
    }
  }
  return(arg_list)
}

# Defaults
args <- list(
    base_dir = "data/arabidopsis/GSE82041_RAW/",
    rmats_out_dir = "ar_results/arabidopsis/nucleus/step01_data_prep/",
    outdir = "ar_results/arabidopsis/nucleus/step_microsome_enrichment/",
    event_group = "RI"
)

# Override with CLI args
cli_args <- get_cli_args()
for (key in names(cli_args)) {
    if (key %in% names(args)) {
        args[[key]] <- cli_args[[key]]
    }
}

BASE_DIR <- args$base_dir
RMATS_OUT_DIR <- args$rmats_out_dir
OUT_DIR <- args$outdir
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


if (args$event_group == "RI") {
    EVENT_TYPE_FILTER <- "RI"
} else {
    EVENT_TYPE_FILTER <- c("SE", "MXE", "A3SS", "A5SS")
}

message(">>> Analysis Settings:")
message("   Base Dir:      ", BASE_DIR)
message("   rMATS Dir:     ", RMATS_OUT_DIR)
message("   Event Group:   ", args$event_group)
message("   Event Types:   ", paste(EVENT_TYPE_FILTER, collapse=", "))

# ==========================================
# 2. DESeq2 LOGIC (RUN ONCE, THEN LOAD)
# ==========================================
rdata_file <- file.path(RMATS_OUT_DIR, "DESeq2_results.RData") # Save in rMATS dir (or base dir if preferred, but base might be read-only)

if (file.exists(rdata_file)) {
  message(">>> Loading pre-calculated DESeq2 results from: ", rdata_file)
  load(rdata_file)
  
} else {
  message(">>> No pre-calculated results found. Starting DESeq2 analysis...")
  
  # A) Read Raw Counts
  files <- list.files(path = BASE_DIR, pattern = "*.txt", full.names = TRUE)
  if(length(files) == 0) stop("No .txt files found in BASE_DIR: ", BASE_DIR)
  
  read_raw_counts <- function(file_path) {
    df <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)
    df <- df[, 1:2] 
    colnames(df) <- c("GeneID", gsub(".txt", "", basename(file_path)))
    return(df)
  }
  
  data_list <- lapply(files, read_raw_counts)
  raw_counts_df <- Reduce(function(x, y) merge(x, y, by = "GeneID", all = TRUE), data_list)
  rownames(raw_counts_df) <- raw_counts_df$GeneID
  raw_counts_df$GeneID <- NULL
  colnames(raw_counts_df) <- sub("^.*mRNA", "mRNA", colnames(raw_counts_df))
  
  # B) Helper for DESeq2
  run_comparison <- function(counts_df, pattern_numerator, pattern_denominator, label_num, label_denom) {
    cols_num <- grep(pattern_numerator, colnames(counts_df), value=TRUE)
    cols_denom <- grep(pattern_denominator, colnames(counts_df), value=TRUE)
    
    if(length(cols_num) == 0 | length(cols_denom) == 0) return(NULL)
    
    subset_counts <- counts_df[, c(cols_denom, cols_num)]
    subset_counts <- subset_counts[rowSums(subset_counts) > 10, ] 
    
    condition <- factor(c(rep(label_denom, length(cols_denom)), rep(label_num, length(cols_num))),
                        levels = c(label_denom, label_num))
    colData <- data.frame(row.names = colnames(subset_counts), condition = condition)
    
    dds <- DESeqDataSetFromMatrix(countData = subset_counts, colData = colData, design = ~ condition)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast=c("condition", label_num, label_denom))
    return(as.data.frame(res))
  }
  
  # C) Run Comparisons
  message("   Running Localization Analysis...")
  df_micro_cyto_raw <- run_comparison(raw_counts_df, "Microsome", "Cytosol", "Microsome", "Cytosol")
  df_mbp_fp_raw     <- run_comparison(raw_counts_df, "RibosomeMBP", "RibosomeFP", "MBP", "FP")
  
  message("   Running TE Analysis...")
  df_te_cyto_raw    <- run_comparison(raw_counts_df, "RibosomeFP", "Cytosol", "FP", "Cytosol")
  df_te_micro_raw   <- run_comparison(raw_counts_df, "RibosomeMBP", "Microsome", "MBP", "Microsome")
  
  # D) Save Results
  save(df_micro_cyto_raw, df_mbp_fp_raw, df_te_cyto_raw, df_te_micro_raw, file = rdata_file)
  message(">>> Analysis complete. Results saved to ", rdata_file)
}

# ==========================================
# 3. LOAD & FILTER EVENT LISTS
# ==========================================
message(">>> Loading Event Lists from RMATS output...")

# UPDATED: Use new nomenclature
file_preserved <- file.path(RMATS_OUT_DIR, "UFM1_independent.tsv")
file_lost      <- file.path(RMATS_OUT_DIR, "UFM1_dependent.tsv")

if(!file.exists(file_preserved) | !file.exists(file_lost)) {
  stop("Could not find UFM1_independent.tsv or UFM1_dependent.tsv in: ", RMATS_OUT_DIR)
}

list_preserved <- readr::read_tsv(file_preserved, show_col_types = FALSE)
list_lost      <- readr::read_tsv(file_lost, show_col_types = FALSE)

# --- FILTERING FUNCTION (Updated for Multiple Types) ---
filter_events <- function(df, ev_types) {
  # 1. Filter by Event Type (Using %in% to allow multiple)
  # "ALL" check not strictly needed if we always pass valid types, but kept for logic
  if (!all(ev_types == "ALL")) {
    df <- df[df$EventType %in% ev_types, ]
  }
  
  # 2. Filter by dPSI (0.2 for SE, 0.1 for others)
  # Note: dPSI filtering is already done in prepare_rmats_data.R, 
  # but re-applying it ensures robustness if raw files were used.
  if (nrow(df) > 0) {
    # Assign dynamic threshold based on the specific row's type
    df$thresh <- ifelse(df$EventType == "SE", 0.2, 0.1)
    df <- df[abs(df$dPSI_num.WT) >= df$thresh, ]
  }
  return(unique(df$GeneID))
}

genes_preserved <- filter_events(list_preserved, EVENT_TYPE_FILTER)
genes_lost      <- filter_events(list_lost, EVENT_TYPE_FILTER)

# Format label for printing
type_label <- paste(EVENT_TYPE_FILTER, collapse="+")
if (args$event_group == "Others") type_label <- "Others" # Shorten for plot title

message(paste0("Event Type(s): ", type_label))
message(paste0("Genes Preserved: ", length(genes_preserved)))
message(paste0("Genes Lost:      ", length(genes_lost)))

# ==========================================
# 4. PLOTTING FUNCTION
# ==========================================
plot_comparison_boxplot <- function(full_res, genes_prev, genes_lost, title_suffix) {
  
  if (is.null(full_res)) {
       plot(1, type="n", axes=F, xlab="", ylab="", main=paste(title_suffix, "(No Data)"))
       return()
  }

  # Extract Log2FC
  vals_genome <- full_res$log2FoldChange
  vals_genome <- vals_genome[!is.na(vals_genome)]
  
  vals_prev <- full_res[rownames(full_res) %in% genes_prev, "log2FoldChange"]
  vals_prev <- vals_prev[!is.na(vals_prev)]
  
  vals_lost <- full_res[rownames(full_res) %in% genes_lost, "log2FoldChange"]
  vals_lost <- vals_lost[!is.na(vals_lost)]
  
  # --- CALCULATE P-VALUES (Wilcoxon Rank Sum Test) ---
  p_prev_gen <- NA
  if(length(vals_prev) > 0 && length(vals_genome) > 0) {
    p_prev_gen <- wilcox.test(vals_prev, vals_genome)$p.value
  }
  
  p_lost_gen <- NA
  if(length(vals_lost) > 0 && length(vals_genome) > 0) {
    p_lost_gen <- wilcox.test(vals_lost, vals_genome)$p.value
  }
  
  fmt_p <- function(p) {
    if(is.na(p)) return("NA")
    if(p < 0.001) return(formatC(p, format = "e", digits = 2))
    return(format(round(p, 3), nsmall = 3))
  }
  
  # Plot
  data_list <- list(Genome = vals_genome, UFM1_Indep = vals_prev, UFM1_Dep = vals_lost)
  cols <- c("grey80", "#2c7bb6", "#d7191c") 
  
  # Safety check for empty lists
  if (length(vals_genome) == 0) vals_genome <- 0
  if (length(vals_prev) == 0) vals_prev <- 0
  if (length(vals_lost) == 0) vals_lost <- 0

  boxplot(data_list, 
          main = title_suffix,
          ylab = "Log2 Fold Change",
          col = cols,
          outline = FALSE,
          notch = TRUE,
          las = 2,
          tcl = -0.3, 
          bty = "n")
  
  abline(h = 0, lty = 2, col = "black")
  
  # Add N counts
  text(x = 1:3, y = par("usr")[3], 
       labels = paste0("n=", c(length(vals_genome), length(vals_prev), length(vals_lost))), 
       pos = 3, cex = 0.8)
  
  # Add P-values at top
  y_max <- par("usr")[4]
  range_y <- y_max - par("usr")[3]
  
  text(x = 2, y = y_max - (range_y * 0.05), 
       labels = paste0("vs Gen: p=", fmt_p(p_prev_gen)), 
       cex = 0.8, font = 3)
  
  text(x = 3, y = y_max - (range_y * 0.05), 
       labels = paste0("vs Gen: p=", fmt_p(p_lost_gen)), 
       cex = 0.8, font = 3)
}

# ==========================================
# 5. GENERATE PLOTS
# ==========================================
outfile <- file.path(OUT_DIR, paste0("microsome_enrichment_", args$event_group, ".pdf"))
pdf(outfile, width=10, height=10)
par(mfrow=c(2,2), las=2, tcl=-0.3, bty="n")

plot_comparison_boxplot(
  df_te_cyto_raw, 
  genes_preserved, 
  genes_lost, 
  paste0("Cytosol TE (", type_label, ")")
)

plot_comparison_boxplot(
  df_te_micro_raw, 
  genes_preserved, 
  genes_lost, 
  paste0("Microsome TE (", type_label, ")")
)

# plot_comparison_boxplot(
#   df_micro_cyto_raw, 
#   genes_preserved, 
#   genes_lost, 
#   paste0("Microsome Enrichment (", type_label, ")")
# )

dev.off()
message(">>> Plotting Complete. Saved to: ", outfile)
