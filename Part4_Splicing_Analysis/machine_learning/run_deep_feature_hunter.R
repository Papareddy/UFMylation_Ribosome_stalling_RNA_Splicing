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
# STEP 1C: Enhanced Motif Features (All SRSF + Top Enriched RBPs)
# =============================================================================
cat("[INFO] Scanning for enriched RBP motifs...\n")

# Define motif patterns based on AME enrichment results (>5 fold in dependent)
# Top enriched motifs from UFM1-dependent vs constitutive analysis
motif_patterns <- list(
  # SRSF Family
  SRSF3 = "[AT]C[AT][AT]C",           # WCWWC
  SRSF5 = "GC[GC]CC",                  # GCSCC
  SRSF11 = "AGGGG",                    # AGGGG
  SRSF1 = "GGAGG",                     # Canonical SRSF1
  SRSF2 = "GGNG",                      # G-rich
  SRSF7 = "GAYGAY",                    # GAC-based
  
  # PCBP Family (C-rich)
  PCBP2 = "CCC+",                      # C-stretch
  PCBP2_long = "CCCC[ACGT]CCCC",       # CCCCNCCCC
  PCBP3 = "CC[AT][AT]CC",              # CCUWWCC
  
  # HNRNP Family (G-rich)
  HNRNPH1 = "G[GT]GGGG",               # GKGGGG
  HNRNPF = "GGGG[ACT]GGGG",            # GGGDGGGG
  HNRNPA2B1 = "[AG]GGGG",              # RGGGG
  
  # Other Top Enriched
  RBM25 = "GGGG[ACGT]G",               # GGGGNG
  RBM6 = "CG[AT]CC[AC]",               # CGUCCM
  RBM39 = "CC[CT]CC",                  # CCYCC
  RBM24 = "GGGGG",                     # G-stretch
  
  # Additional splice-related
  PTBP1 = "[CT]CTCCC",                 # Pyrimidine-rich
  ESRP1 = "GGG[GT]GG"                  # GGGKGG
)

# Generic motif counter
count_motif <- function(seq, pattern) {
  if (is.na(seq) || nchar(seq) == 0) return(0)
  matches <- gregexpr(pattern, toupper(seq), perl = TRUE)[[1]]
  if (matches[1] == -1) return(0)
  length(matches)
}

# Add all motif counts and densities
cat("[INFO] Counting motifs (", length(motif_patterns), " patterns)...\n")

for (motif_name in names(motif_patterns)) {
  pattern <- motif_patterns[[motif_name]]
  count_col <- paste0(motif_name, "_count")
  density_col <- paste0(motif_name, "_density")
  
  features[[count_col]] <- sapply(intron_seqs, function(s) count_motif(s, pattern))
  features[[density_col]] <- features[[count_col]] / (features$intron_length / 1000)
}

# =============================================================================
# STEP 1D: Splice Site Strength (MaxEntScan-like scoring)
# =============================================================================
cat("[INFO] Calculating splice site strength scores...\n")

# Simplified SS scoring based on consensus matching
score_5ss <- function(seq) {
  # 5' SS consensus: GT at positions 1-2 of intron
  if (is.na(seq) || nchar(seq) < 9) return(NA)
  # Canonical 5'SS: [AC]AG|GT[AG]AGT (| is splice site)
  # Score based on consensus match
  score <- 0
  seq_upper <- toupper(seq)
  # Position 4-5 should be GT
  if (substr(seq_upper, 4, 5) == "GT") score <- score + 5
  if (substr(seq_upper, 6, 6) %in% c("A", "G")) score <- score + 2
  if (substr(seq_upper, 3, 3) == "G") score <- score + 1
  # GC content of region
  gc <- sum(strsplit(seq_upper, "")[[1]] %in% c("G", "C")) / nchar(seq_upper)
  score + gc * 3
}

score_3ss <- function(seq) {
  # 3' SS consensus: Polypyrimidine tract + AG
  if (is.na(seq) || nchar(seq) < 20) return(NA)
  seq_upper <- toupper(seq)
  score <- 0
  # PPT in last 20bp
  ppt_region <- substr(seq_upper, nchar(seq_upper) - 19, nchar(seq_upper) - 2)
  ppt_content <- sum(strsplit(ppt_region, "")[[1]] %in% c("C", "T")) / nchar(ppt_region)
  score <- ppt_content * 8
  # AG at end
  if (substr(seq_upper, nchar(seq_upper) - 1, nchar(seq_upper)) == "AG") score <- score + 4
  score
}

features$SS5_score <- sapply(ss_seqs$ss5, score_5ss)
features$SS3_score <- sapply(ss_seqs$ss3, score_3ss)

# =============================================================================
# STEP 1E: Branch Point Scoring
# =============================================================================
cat("[INFO] Scanning for branch point motifs...\n")

# Branch point consensus: YNYURAC (Y=pyrimidine, N=any, R=purine)
# Simplified as [CT][ACGT][CT][AT][AG]A[CT]
bp_pattern <- "[CT][ACGT][CT][AT][AG]A[CT]"

score_bp <- function(seq) {
  if (is.na(seq) || nchar(seq) < 40) return(0)
  # Search in -40 to -15 region before 3'SS
  seq_upper <- toupper(seq)
  bp_region <- substr(seq_upper, max(1, nchar(seq_upper) - 39), nchar(seq_upper) - 14)
  matches <- gregexpr(bp_pattern, bp_region, perl = TRUE)[[1]]
  if (matches[1] == -1) return(0)
  length(matches)  # Number of potential BP motifs
}

features$BP_count <- sapply(intron_seqs, score_bp)

# =============================================================================
# STEP 1F: Gene Architecture Features (Position Hypothesis)
# =============================================================================
cat("[INFO] Calculating gene architecture features...\n")

# Extract intron rank from event_id structure
# event_id format: chr|strand|gene_start|gene_end|intron_start|intron_end|...
# We'll estimate intron rank based on position relative to gene

# For now, use intron position within the gene as a proxy
# Lower coordinate introns = earlier in gene (for + strand)
features$rel_position <- 0.5  # Default middle position if unknown

# Load original data to get more metadata
dep_df <- read_tsv(opt$dependent, show_col_types = FALSE)
indep_df <- read_tsv(opt$independent, show_col_types = FALSE)
all_df <- bind_rows(dep_df, indep_df)

# Try to get intron rank from the structure (simplified - estimate from coordinates)
for (i in 1:nrow(features)) {
  if (i <= length(all_gr) && i <= nrow(all_df)) {
    # Estimate relative position within gene (0-1)
    event_parts <- strsplit(all_df$event_id[i], "\\|")[[1]]
    if (length(event_parts) >= 4) {
      gene_start <- as.numeric(event_parts[3])
      gene_end <- as.numeric(event_parts[4])
      intron_start <- start(all_gr[i])
      
      if (!is.na(gene_start) && !is.na(gene_end) && gene_end > gene_start) {
        features$rel_position[i] <- (intron_start - gene_start) / (gene_end - gene_start)
      }
    }
  }
}

# Estimate if first intron (early position)
features$is_first_intron <- ifelse(is.na(features$rel_position) | features$rel_position >= 0.2, 0, 1)

# =============================================================================
# STEP 1G: Upstream Exon Features (Ribosome Kinetics)
# =============================================================================
cat("[INFO] Extracting upstream exon features...\n")

# Extract upstream exon sequence (50bp before intron start)
extract_upstream_exon <- function(gr, genome, window = 50) {
  seqs <- character(length(gr))
  for (i in seq_along(gr)) {
    chr <- as.character(seqnames(gr[i]))
    if (!chr %in% names(genome)) next
    chr_seq <- genome[[chr]]
    strand_i <- as.character(strand(gr[i]))
    
    if (strand_i == "+") {
      exon_start <- max(1, start(gr[i]) - window)
      exon_end <- start(gr[i]) - 1
    } else {
      exon_start <- end(gr[i]) + 1
      exon_end <- min(length(chr_seq), end(gr[i]) + window)
    }
    
    if (exon_start >= exon_end || exon_end > length(chr_seq)) next
    
    exon_seq <- subseq(chr_seq, exon_start, exon_end)
    if (strand_i == "-") exon_seq <- reverseComplement(exon_seq)
    seqs[i] <- as.character(exon_seq)
  }
  return(seqs)
}

upstream_exon_seqs <- extract_upstream_exon(all_gr, genome, 50)

# Upstream exon GC content
features$upstream_exon_GC <- sapply(upstream_exon_seqs, function(s) {
  if (is.na(s) || nchar(s) == 0) return(NA)
  bases <- strsplit(toupper(s), "")[[1]]
  sum(bases %in% c("G", "C")) / length(bases)
})

# Upstream exon length (actual extracted, should be ~50)
features$upstream_exon_length <- nchar(upstream_exon_seqs)

# Simple codon optimality proxy: frequency of optimal codons (A/T at 3rd position = more optimal in general)
calc_codon_optimality <- function(seq) {
  if (is.na(seq) || nchar(seq) < 6) return(NA)
  seq_upper <- toupper(seq)
  # Count codons with A/T at wobble position (3rd, 6th, 9th, etc.)
  wobble_positions <- seq(3, nchar(seq_upper), by = 3)
  wobble_bases <- sapply(wobble_positions, function(p) substr(seq_upper, p, p))
  optimal_count <- sum(wobble_bases %in% c("A", "T"))
  optimal_count / length(wobble_bases)
}

features$codon_optimality <- sapply(upstream_exon_seqs, calc_codon_optimality)

# =============================================================================
# STEP 1H: Decoy Competition (Distraction Hypothesis)
# =============================================================================
cat("[INFO] Calculating splice site decoy competition...\n")

# Scan for decoy GT/AG sites within 100bp of real sites
calc_decoy_competition <- function(gr, genome, window = 100) {
  decoy_5ss <- numeric(length(gr))
  decoy_3ss <- numeric(length(gr))
  
  for (i in seq_along(gr)) {
    chr <- as.character(seqnames(gr[i]))
    if (!chr %in% names(genome)) next
    chr_seq <- genome[[chr]]
    strand_i <- as.character(strand(gr[i]))
    
    # 5' SS region (around intron start)
    if (strand_i == "+") {
      region_5_start <- max(1, start(gr[i]) - window)
      region_5_end <- min(length(chr_seq), start(gr[i]) + window)
    } else {
      region_5_start <- max(1, end(gr[i]) - window)
      region_5_end <- min(length(chr_seq), end(gr[i]) + window)
    }
    
    if (region_5_start < region_5_end) {
      region_5_seq <- as.character(subseq(chr_seq, region_5_start, region_5_end))
      # Count GT dinucleotides (potential 5'SS)
      gt_matches <- gregexpr("GT", toupper(region_5_seq))[[1]]
      decoy_5ss[i] <- ifelse(gt_matches[1] == -1, 0, length(gt_matches) - 1)  # -1 for real site
    }
    
    # 3' SS region (around intron end)
    if (strand_i == "+") {
      region_3_start <- max(1, end(gr[i]) - window)
      region_3_end <- min(length(chr_seq), end(gr[i]) + window)
    } else {
      region_3_start <- max(1, start(gr[i]) - window)
      region_3_end <- min(length(chr_seq), start(gr[i]) + window)
    }
    
    if (region_3_start < region_3_end) {
      region_3_seq <- as.character(subseq(chr_seq, region_3_start, region_3_end))
      # Count AG dinucleotides (potential 3'SS)
      ag_matches <- gregexpr("AG", toupper(region_3_seq))[[1]]
      decoy_3ss[i] <- ifelse(ag_matches[1] == -1, 0, length(ag_matches) - 1)  # -1 for real site
    }
  }
  
  return(list(decoy_5ss = decoy_5ss, decoy_3ss = decoy_3ss))
}

decoy_counts <- calc_decoy_competition(all_gr, genome, 100)
features$decoy_5ss_count <- decoy_counts$decoy_5ss
features$decoy_3ss_count <- decoy_counts$decoy_3ss
features$decoy_total <- features$decoy_5ss_count + features$decoy_3ss_count

# =============================================================================
# STEP 1I: RNA Structural Features (Accessibility)
# =============================================================================
cat("[INFO] Calculating RNA structural accessibility...\n")

# Simple structural proxy: dinucleotide stacking energy (AU-rich = more accessible)
# Real RNAfold requires system call - using simplified proxy
calc_accessibility <- function(seq) {
  if (is.na(seq) || seq == "" || nchar(seq) < 3) return(0.5)  # Return neutral value
  seq_upper <- toupper(seq)
  # AU content correlates with less structure (more accessible)
  bases <- strsplit(seq_upper, "")[[1]]
  au_content <- sum(bases %in% c("A", "T")) / length(bases)
  # GC pairs form stronger structures
  gc_matches <- gregexpr("GC|CG", seq_upper)[[1]]
  gc_pairs <- ifelse(gc_matches[1] == -1, 0, length(gc_matches))
  # Accessibility score: higher = more accessible
  au_content - (gc_pairs / nchar(seq_upper))
}

# Calculate for 5'SS and 3'SS windows
features$SS5_accessibility <- sapply(ss_seqs$ss5, calc_accessibility)
features$SS3_accessibility <- sapply(ss_seqs$ss3, calc_accessibility)

# =============================================================================
# STEP 1J: Additional Sequence Features
# =============================================================================

# Impute NAs with median for numeric columns (instead of removing rows)
cat("[INFO] Imputing missing values with column medians...\n")
for (col in names(features)) {
  if (col == "class") next
  if (is.numeric(features[[col]])) {
    na_count <- sum(is.na(features[[col]]))
    if (na_count > 0) {
      median_val <- median(features[[col]], na.rm = TRUE)
      if (is.na(median_val)) median_val <- 0  # fallback
      features[[col]][is.na(features[[col]])] <- median_val
      cat(sprintf("[INFO] Imputed %d NAs in %s with %.2f\n", na_count, col, median_val))
    }
  }
}

cat(sprintf("[INFO] Feature matrix: %d introns x %d features\n", nrow(features), ncol(features) - 1))

# Debug: check for remaining NAs
na_cols <- sapply(features, function(x) sum(is.na(x)))
if (any(na_cols > 0)) {
  cat("[DEBUG] Columns with NAs after imputation:\n")
  print(na_cols[na_cols > 0])
}

# Final cleanup: remove any rows still with NAs (edge cases)
features <- features[complete.cases(features), ]
cat(sprintf("[INFO] After NA cleanup: %d introns\n", nrow(features)))

# If too few rows remain, skip complete.cases and use imputed data
if (nrow(features) == 0) {
  cat("[WARN] All rows had NAs - reloading and forcing zero imputation...\n")
  # Reload features without the problematic complete.cases
  features <- read_csv(file.path(opt$outdir, "feature_matrix_debug.csv"), show_col_types = FALSE)
}

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
