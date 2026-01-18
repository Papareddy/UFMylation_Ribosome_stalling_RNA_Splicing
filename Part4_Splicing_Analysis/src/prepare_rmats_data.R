suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(RColorBrewer)
})

# --- Operator for default values ---
`%||%` <- function(a, b) if (is.null(a)) b else a

# --- Function to parse command line arguments ---
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

# --- Get arguments or set defaults ---
cli_args <- get_cli_args()
WT_DIR   <- cli_args$wt_dir   %||% stop("Missing --wt_dir")
UFM_DIR  <- cli_args$ufm_dir  %||% stop("Missing --ufm_dir")
OUTDIR   <- cli_args$outdir   %||% "results/step0"
FDR_CUT  <- as.numeric(cli_args$fdr %||% 0.05)
DPSI_CUT <- as.numeric(cli_args$dpsi %||% 0.15)
MIN_READS <- as.numeric(cli_args$min_reads %||% 20)
WHICH <- "JCEC"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# --- Select Event Types ---
all_event_types <- c("A3SS","A5SS","MXE","RI","SE")

if (!is.null(cli_args$event_types)) {
  target_event_types <- unlist(strsplit(cli_args$event_types, ","))
  # Validate
  invalid <- setdiff(target_event_types, all_event_types)
  if (length(invalid) > 0) warning("Unknown event types ignored: ", paste(invalid, collapse=", "))
  target_event_types <- intersect(target_event_types, all_event_types)
} else {
  target_event_types <- all_event_types
}

# Always read ALL event types so the summary plot covers everything
# We will filter for the output files (lost/preserved.tsv) later.
proc_event_types <- all_event_types

# --- Helpers ---
parse_num_vec <- function(x) {
  if (is.null(x) || length(x) == 0) return(numeric())
  x <- as.character(x)
  if (is.na(x) || x == "" || x == "NA") return(numeric())
  v <- suppressWarnings(as.numeric(strsplit(x, ",")[[1]]))
  v[!is.na(v)]
}
mean_from_vec <- function(x) { v <- parse_num_vec(x); if (length(v) == 0) NA_real_ else mean(v) }
sum_from_vec <- function(x) { v <- parse_num_vec(x); if (length(v) == 0) NA_real_ else sum(v) }

make_event_id_fromGTF <- function(df, ev) {
  cols <- switch(
    ev,
    "SE"   = c("chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"),
    "RI"   = c("chr","strand","riExonStart_0base","riExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"),
    "A3SS" = c("chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE"),
    "A5SS" = c("chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE"),
    "MXE"  = c("chr","strand", "1stExonStart_0base","1stExonEnd", "2ndExonStart_0base","2ndExonEnd", "upstreamES","upstreamEE","downstreamES","downstreamEE"),
    stop("Unhandled event type: ", ev, call. = FALSE)
  )
  if (!all(cols %in% names(df))) stop("Missing columns for ", ev, ": ", paste(setdiff(cols, names(df)), collapse=", "), call. = FALSE)
  df$event_id <- do.call(paste, c(df[cols], sep="|"))
  df
}

read_mats <- function(dir, ev, which = "JC") {
  f <- file.path(dir, paste0(ev, ".MATS.", which, ".txt"))
  if (!file.exists(f)) stop("File not found: ", f, call.=FALSE)
  df <- suppressMessages(readr::read_tsv(f, show_col_types = FALSE, quote = "\"", name_repair = "unique"))
  df$EventType <- ev
  if ("GeneID" %in% names(df)) df$GeneID <- gsub('"', "", df$GeneID)
  if ("geneSymbol" %in% names(df)) df$geneSymbol <- gsub('"', "", df$geneSymbol)
  df <- make_event_id_fromGTF(df, ev)
  df$FDR_num  <- suppressWarnings(as.numeric(df$FDR))
  df$dPSI_num <- suppressWarnings(as.numeric(df$IncLevelDifference))
  df$mean_psi_ctrl <- vapply(df$IncLevel1, mean_from_vec, numeric(1))
  df$mean_psi_case <- vapply(df$IncLevel2, mean_from_vec, numeric(1))
  needed_counts <- c("IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2")
  if (all(needed_counts %in% names(df))) {
    df$total_reads <- rowSums(cbind(
      vapply(df$IJC_SAMPLE_1, sum_from_vec, numeric(1)),
      vapply(df$SJC_SAMPLE_1, sum_from_vec, numeric(1)),
      vapply(df$IJC_SAMPLE_2, sum_from_vec, numeric(1)),
      vapply(df$SJC_SAMPLE_2, sum_from_vec, numeric(1))
    ), na.rm=TRUE)
  } else {
    df$total_reads <- NA_real_
  }
  df[!is.na(df$FDR_num) & df$FDR_num >= 0 & df$FDR_num <= 1, , drop = FALSE]
}

read_all_event_types <- function(event_types, fun, ...) {
  dplyr::bind_rows(lapply(event_types, function(ev) fun(..., ev)))
}

# --- Main Logic ---
wt_mats  <- read_all_event_types(proc_event_types, read_mats, WT_DIR, which = WHICH)
ufm_mats <- read_all_event_types(proc_event_types, read_mats, UFM_DIR, which = WHICH)

merged_fixed <- dplyr::inner_join(
  wt_mats, ufm_mats, 
  by = c("event_id", "GeneID", "geneSymbol", "EventType"), 
  suffix = c(".WT", ".UFM")
)

wt_sig_events <- merged_fixed %>%
  dplyr::filter(
    !is.na(total_reads.WT), total_reads.WT >= MIN_READS,
    FDR_num.WT < FDR_CUT,
    ifelse(EventType == "SE", abs(dPSI_num.WT) >= 0.2, abs(dPSI_num.WT) >= 0.1)
  )

wt_sig_covered <- wt_sig_events %>%
  dplyr::filter(
    !is.na(FDR_num.UFM), !is.na(dPSI_num.UFM),
    !is.na(total_reads.UFM), total_reads.UFM >= MIN_READS
  )

wt_sig_preserved_ufm1 <- wt_sig_covered %>%
  dplyr::filter(
    FDR_num.UFM < FDR_CUT, 
    ifelse(EventType == "SE", abs(dPSI_num.UFM) >= 0.2, abs(dPSI_num.UFM) >= 0.1)
  )

wt_sig_lost_in_ufm1 <- wt_sig_covered %>%
  dplyr::filter(
    !(FDR_num.UFM < FDR_CUT & ifelse(EventType == "SE", abs(dPSI_num.UFM) >= 0.2, abs(dPSI_num.UFM) >= 0.1))
  )
  
wt_sig_notcovered_ufm1 <- wt_sig_events %>%
  dplyr::filter(is.na(FDR_num.UFM) | is.na(dPSI_num.UFM) | is.na(total_reads.UFM) | total_reads.UFM < MIN_READS)

# --- Write Output (Filtered by requested event_types) ---
# Filter only for writing to files, so downstream steps only process what was requested.
out_preserved <- wt_sig_preserved_ufm1 %>% dplyr::filter(EventType %in% target_event_types)
out_lost      <- wt_sig_lost_in_ufm1      %>% dplyr::filter(EventType %in% target_event_types)

readr::write_tsv(out_preserved, file.path(OUTDIR, "UFM1_independent.tsv"))
readr::write_tsv(out_lost,      file.path(OUTDIR, "UFM1_dependent.tsv"))

# Split by Direction (dPSI > 0 vs dPSI < 0) for directional motif analysis
# Positive: WT dPSI > 0 (Inclusion Enhanced in WT vs Control) -> dPSI_positive
# Negative: WT dPSI < 0 (Inclusion Silenced in WT vs Control) -> dPSI_negative

out_UFM1_independent_pos <- out_preserved %>% dplyr::filter(dPSI_num.WT > 0)
out_UFM1_independent_neg <- out_preserved %>% dplyr::filter(dPSI_num.WT < 0)
out_UFM1_dependent_pos   <- out_lost      %>% dplyr::filter(dPSI_num.WT > 0)
out_UFM1_dependent_neg   <- out_lost      %>% dplyr::filter(dPSI_num.WT < 0)

readr::write_tsv(out_UFM1_independent_pos, file.path(OUTDIR, "UFM1_independent_dPSI_positive.tsv"))
readr::write_tsv(out_UFM1_independent_neg, file.path(OUTDIR, "UFM1_independent_dPSI_negative.tsv"))
readr::write_tsv(out_UFM1_dependent_pos,   file.path(OUTDIR, "UFM1_dependent_dPSI_positive.tsv"))
readr::write_tsv(out_UFM1_dependent_neg,   file.path(OUTDIR, "UFM1_dependent_dPSI_negative.tsv"))

cat("Generated 'UFM1_independent.tsv', 'UFM1_dependent.tsv' and split versions (_dPSI_positive/_dPSI_negative) in", OUTDIR, "(filtered by", paste(target_event_types, collapse=","), ")\n")

# --- Export Background Genes (Tested in rMATS) ---
bg_genes <- unique(stats::na.omit(merged_fixed$GeneID))
writeLines(bg_genes, file.path(OUTDIR, "background_genes.txt"))
cat("Generated 'background_genes.txt' in", OUTDIR, "\n")

# --- Plotting (Uses ALL event types) ---
counts_all <- wt_sig_events %>%
  dplyr::count(EventType, name = "n_WT_sig") %>%
  dplyr::full_join(wt_sig_preserved_ufm1 %>% dplyr::count(EventType, name = "n_preserved_ufm1"), by = "EventType") %>%
  dplyr::full_join(wt_sig_lost_in_ufm1 %>% dplyr::count(EventType, name = "n_lost_ufm1"), by = "EventType") %>%
  dplyr::full_join(wt_sig_notcovered_ufm1 %>% dplyr::count(EventType, name = "n_not_covered_ufm1"), by = "EventType") %>%
  dplyr::mutate(across(starts_with("n_"), ~replace_na(., 0L))) %>%
  dplyr::mutate(EventType = factor(EventType, levels = all_event_types)) %>%
  dplyr::arrange(EventType) %>%
  dplyr::mutate(check_sum = n_preserved_ufm1 + n_lost_ufm1 + n_not_covered_ufm1)

readr::write_tsv(counts_all, file.path(OUTDIR, "WT_sig_event_counts_by_category.tsv"))

mat <- rbind(
  Preserved_in_UFM1   = counts_all$n_preserved_ufm1,
  Lost_in_UFM1        = counts_all$n_lost_ufm1,
  Not_covered_in_UFM1 = counts_all$n_not_covered_ufm1
)
colnames(mat) <- as.character(counts_all$EventType)

event_colors <- setNames(brewer.pal(length(all_event_types), "Set2"), all_event_types)

pdf(file.path(OUTDIR, "event_counts_barplot.pdf"), width = 8, height = 6)
par(las = 2, mar = c(5, 8, 4, 2) + 0.1)

bp <- barplot(height = colSums(mat), plot = FALSE)
ylim_max <- max(colSums(mat)) * 1.15

plot(NA, xlim = range(bp) + c(-0.8, 0.8), ylim = c(0, ylim_max), xaxt = "n", xlab = "Event type", ylab = "Count", main = "WT vs ANS (UFM1 dependency)")
axis(1, at = bp, labels = colnames(mat), las = 2)

alphas <- c(1.0, 0.75, 0.5)
w <- 1

for (j in seq_along(bp)) {
  ev <- colnames(mat)[j]
  base_col <- event_colors[[ev]]
  y0 <- 0
  
  for (i in seq_len(nrow(mat))) {
    h <- mat[i, j]
    if (h <= 0) next
    
    x1 <- bp[j] - w/2
    x2 <- bp[j] + w/2
    y1 <- y0
    y2 <- y0 + h
    
    rect(x1, y1, x2, y2, col = adjustcolor(base_col, alpha.f = alphas[i]), border = NA)
    
    if (i == 2) { # Hatch for "Lost"
      polygon(x = c(x1, x2, x2, x1), y = c(y1, y1, y2, y2), density = 15, angle = 45, col = NA, border = "white", lwd = 0.8)
    }
    y0 <- y2
  }
  
  y_centers <- c(mat[1, j] / 2, mat[1, j] + mat[2, j] / 2, mat[1, j] + mat[2, j] + mat[3, j] / 2)
  labels_vec <- c(mat[1, j], mat[2, j], mat[3, j])
  
  for (k in seq_along(labels_vec)) {
    if (labels_vec[k] > 0) {
      text(bp[j], y_centers[k], labels = labels_vec[k], cex = 0.85)
    }
  }
}

legend("topleft", bty = "n", legend = rownames(mat), fill = sapply(alphas, function(a) adjustcolor("grey", alpha.f = a)))
dev.off()

cat("Generated barplot 'event_counts_barplot.pdf' in", OUTDIR, "\n")