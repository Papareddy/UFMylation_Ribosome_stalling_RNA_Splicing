#!/usr/bin/env Rscript
# =============================================================================
# Deep Feature Hunter: ML Pipeline for UFM1-Dependent vs Independent Introns
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(Biostrings)
  library(rtracklayer)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(randomForest)
  library(pROC)
  library(caret)
})

# --- CLI Arguments ---
option_list <- list(
  make_option("--dependent", type = "character", help = "Path to UFM1_dependent.tsv"),
  make_option("--independent", type = "character", help = "Path to UFM1_independent.tsv"),
  make_option("--genome", type = "character", help = "Path to genome FASTA"),
  make_option("--meme", type = "character", default = NULL, help = "Path to .meme motif file"),
  make_option("--outdir", type = "character", default = "machine_learning/outputs", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
cat("[INFO] Deep Feature Hunter Pipeline Starting...\n")

# =============================================================================
# STEP 0: Load Data
# =============================================================================
cat("[INFO] Loading intron data...\n")

load_introns <- function(tsv_file, label) {
  df <- read_tsv(tsv_file, show_col_types = FALSE)
  
  # Extract RI coordinates from event_id
  coords <- strsplit(df$event_id, "\\|")
  
  gr <- GRanges(
    seqnames = sapply(coords, `[`, 1),
    strand = sapply(coords, `[`, 2),
    ranges = IRanges(
      start = as.integer(sapply(coords, `[`, 3)),
      end = as.integer(sapply(coords, `[`, 4))
    )
  )
  
  mcols(gr)$gene_id <- df$GeneID
  mcols(gr)$gene_symbol <- df$geneSymbol
  mcols(gr)$dpsi <- df$dPSI_num.WT
  mcols(gr)$class <- label
  
  return(gr)
}

dep_gr <- load_introns(opt$dependent, "Dependent")
indep_gr <- load_introns(opt$independent, "Independent")

cat(sprintf("[INFO] Loaded %d Dependent and %d Independent introns.\n", length(dep_gr), length(indep_gr)))

# Combine
all_gr <- c(dep_gr, indep_gr)
all_gr$class_label <- ifelse(all_gr$class == "Dependent", 1, 0)

# =============================================================================
# STEP 1A: Splice Site Sequence Extraction
# =============================================================================
cat("[INFO] Loading genome...\n")
genome <- readDNAStringSet(opt$genome)
names(genome) <- gsub(" .*", "", names(genome))

# Harmonize chromosome names (remove 'chr' prefix from intron data if genome uses numbers)
genome_has_chr <- any(grepl("^chr", names(genome)))
cat(sprintf("[INFO] Genome has 'chr' prefix: %s\n", genome_has_chr))

# Fix chromosome names in GRanges to match genome
harmonize_chr <- function(gr, genome_names) {
  gr_chr <- as.character(seqnames(gr))
  gr_has_chr <- any(grepl("^chr", gr_chr))
  genome_has_chr <- any(grepl("^chr", genome_names))
  
  if (gr_has_chr && !genome_has_chr) {
    # Remove 'chr' from GRanges
    new_chr <- gsub("^chr", "", gr_chr)
    seqlevels(gr) <- unique(new_chr)
    gr <- GRanges(seqnames = new_chr, ranges = ranges(gr), strand = strand(gr), mcols(gr))
  } else if (!gr_has_chr && genome_has_chr) {
    # Add 'chr' to GRanges
    new_chr <- paste0("chr", gr_chr)
    seqlevels(gr) <- unique(new_chr)
    gr <- GRanges(seqnames = new_chr, ranges = ranges(gr), strand = strand(gr), mcols(gr))
  }
  return(gr)
}

all_gr <- harmonize_chr(all_gr, names(genome))
cat(sprintf("[INFO] Chromosome name example: %s\n", as.character(seqnames(all_gr)[1])))

extract_ss_sequences <- function(gr, genome) {
  # 5' SS: -3 to +6 (9-mer) relative to intron start
  # 3' SS: -20 to +3 (23-mer) relative to intron end
  
  ss5_seqs <- character(length(gr))
  ss3_seqs <- character(length(gr))
  
  for (i in seq_along(gr)) {
    chr <- as.character(seqnames(gr[i]))
    strand_i <- as.character(strand(gr[i]))
    
    if (!chr %in% names(genome)) next
    chr_seq <- genome[[chr]]
    
    if (strand_i == "+") {
      # 5' SS: exon(-3) + intron(+6)
      ss5_start <- start(gr[i]) - 3
      ss5_end <- start(gr[i]) + 5
      # 3' SS: intron(-20) + exon(+3)
      ss3_start <- end(gr[i]) - 19
      ss3_end <- end(gr[i]) + 3
    } else {
      # Reverse complement for minus strand
      ss5_start <- end(gr[i]) - 5
      ss5_end <- end(gr[i]) + 3
      ss3_start <- start(gr[i]) - 3
      ss3_end <- start(gr[i]) + 19
    }
    
    # Bounds check
    if (ss5_start < 1 || ss5_end > length(chr_seq)) next
    if (ss3_start < 1 || ss3_end > length(chr_seq)) next
    
    ss5_seq <- subseq(chr_seq, ss5_start, ss5_end)
    ss3_seq <- subseq(chr_seq, ss3_start, ss3_end)
    
    if (strand_i == "-") {
      ss5_seq <- reverseComplement(ss5_seq)
      ss3_seq <- reverseComplement(ss3_seq)
    }
    
    ss5_seqs[i] <- as.character(ss5_seq)
    ss3_seqs[i] <- as.character(ss3_seq)
  }
  
  return(list(ss5 = ss5_seqs, ss3 = ss3_seqs))
}

cat("[INFO] Extracting splice site sequences...\n")
ss_seqs <- extract_ss_sequences(all_gr, genome)

# =============================================================================
# STEP 1B: Sequence Features
# =============================================================================
cat("[INFO] Calculating sequence features...\n")

calc_gc_content <- function(seq) {
  if (nchar(seq) == 0) return(NA)
  bases <- strsplit(toupper(seq), "")[[1]]
  sum(bases %in% c("G", "C")) / length(bases)
}

calc_ppt_density <- function(seq) {
  # Pyrimidine (C/T) density in last 15bp before 3' SS
  if (nchar(seq) < 15) return(NA)
  region <- substr(seq, nchar(seq) - 14, nchar(seq))
  bases <- strsplit(toupper(region), "")[[1]]
  sum(bases %in% c("C", "T")) / length(bases)
}

calc_t_block <- function(seq) {
  # Longest run of Ts
  if (nchar(seq) == 0) return(0)
  runs <- rle(strsplit(toupper(seq), "")[[1]])
  t_runs <- runs$lengths[runs$values == "T"]
  if (length(t_runs) == 0) return(0)
  max(t_runs)
}

extract_intron_sequence <- function(gr, genome) {
  seqs <- character(length(gr))
  for (i in seq_along(gr)) {
    chr <- as.character(seqnames(gr[i]))
    if (!chr %in% names(genome)) next
    chr_seq <- genome[[chr]]
    s <- max(1, start(gr[i]))
    e <- min(length(chr_seq), end(gr[i]))
    if (s >= e) next
    intron_seq <- subseq(chr_seq, s, e)
    if (as.character(strand(gr[i])) == "-") {
      intron_seq <- reverseComplement(intron_seq)
    }
    seqs[i] <- as.character(intron_seq)
  }
  return(seqs)
}

intron_seqs <- extract_intron_sequence(all_gr, genome)

# Build feature matrix
features <- data.frame(
  class = all_gr$class_label,
  intron_length = width(all_gr),
  GC_intron = sapply(intron_seqs, calc_gc_content),
  PPT_density = sapply(ss_seqs$ss3, calc_ppt_density),
  T_block_length = sapply(intron_seqs, calc_t_block)
)

# =============================================================================
# STEP 1C: Motif Features (SRSF3/PCBP2)
# =============================================================================
cat("[INFO] Scanning for SRSF3 (WCWWC) motifs...\n")

# SRSF3 canonical: WCWWC (W = A/T)
srsf3_pattern <- "([AT]C[AT][AT]C)"

count_srsf3 <- function(seq) {
  if (is.na(seq) || nchar(seq) == 0) return(0)
  length(gregexpr(srsf3_pattern, toupper(seq), perl = TRUE)[[1]])
}

features$SRSF3_count <- sapply(intron_seqs, count_srsf3)
features$SRSF3_density <- features$SRSF3_count / (features$intron_length / 1000)

# PCBP2 pattern: C-rich (CCUYCCC or UWCCC simplified as CCC+)
pcbp2_pattern <- "CCC+"

count_pcbp2 <- function(seq) {
  if (is.na(seq) || nchar(seq) == 0) return(0)
  matches <- gregexpr(pcbp2_pattern, toupper(seq), perl = TRUE)[[1]]
  if (matches[1] == -1) return(0)
  length(matches)
}

features$PCBP2_count <- sapply(intron_seqs, count_pcbp2)
features$PCBP2_density <- features$PCBP2_count / (features$intron_length / 1000)

# =============================================================================
# STEP 1D: Additional Features
# =============================================================================

# CpG count
calc_cpg <- function(seq) {
  if (is.na(seq) || nchar(seq) == 0) return(0)
  length(gregexpr("CG", toupper(seq))[[1]])
}

features$CpG_count <- sapply(intron_seqs, calc_cpg)
features$CpG_density <- features$CpG_count / (features$intron_length / 1000)

# Remove rows with NAs
features <- features[complete.cases(features), ]
cat(sprintf("[INFO] Feature matrix: %d introns x %d features\n", nrow(features), ncol(features) - 1))

# Save feature matrix
write_csv(features, file.path(opt$outdir, "feature_matrix.csv"))

# =============================================================================
# STEP 2: Machine Learning (Random Forest)
# =============================================================================
cat("[INFO] Training Random Forest classifier...\n")

features$class <- as.factor(features$class)

# 5-fold CV
set.seed(42)
train_control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

# Rename classes for caret
levels(features$class) <- c("Independent", "Dependent")

rf_model <- train(
  class ~ . - class,
  data = features,
  method = "rf",
  trControl = train_control,
  metric = "ROC",
  ntree = 500
)

cat(sprintf("[INFO] Cross-validated AUC: %.3f\n", max(rf_model$results$ROC)))

# =============================================================================
# STEP 3A: Variable Importance Plot
# =============================================================================
cat("[INFO] Generating Variable Importance Plot...\n")

importance_df <- varImp(rf_model)$importance
importance_df$Feature <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$Overall), ]

p_imp <- ggplot(importance_df, aes(x = reorder(Feature, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Variable Importance: UFM1-Dependent vs Independent",
    x = "Feature",
    y = "Importance (Mean Decrease Gini)"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(opt$outdir, "Variable_Importance.pdf"), p_imp, width = 8, height = 6)

# =============================================================================
# STEP 3B: Feature Boxplots
# =============================================================================
cat("[INFO] Generating Feature Effect Boxplots...\n")

top_features <- head(importance_df$Feature, 4)
features_long <- features %>%
  select(class, all_of(top_features)) %>%
  tidyr::pivot_longer(cols = -class, names_to = "Feature", values_to = "Value")

p_box <- ggplot(features_long, aes(x = class, y = Value, fill = class)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~Feature, scales = "free_y") +
  scale_fill_manual(values = c("Dependent" = "red", "Independent" = "blue")) +
  labs(title = "Top 4 Discriminative Features", x = "", y = "Value") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(opt$outdir, "Feature_Boxplots.pdf"), p_box, width = 10, height = 8)

# =============================================================================
# STEP 3C: ROC Curve
# =============================================================================
cat("[INFO] Generating ROC Curve...\n")

pred_probs <- rf_model$pred
roc_obj <- roc(pred_probs$obs, pred_probs$Dependent)

pdf(file.path(opt$outdir, "ROC_Curve.pdf"), width = 6, height = 6)
plot(roc_obj, main = sprintf("ROC Curve (AUC = %.3f)", auc(roc_obj)), col = "darkblue", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "grey")
dev.off()

# =============================================================================
# STEP 3D: Positional Metagene (SRSF3 density)
# =============================================================================
cat("[INFO] Generating Positional Metagene Profile...\n")

# This would require more complex positional scanning
# Placeholder for now

cat("[DONE] Deep Feature Hunter Pipeline Completed!\n")
cat(sprintf("[DONE] Results saved to: %s\n", opt$outdir))
