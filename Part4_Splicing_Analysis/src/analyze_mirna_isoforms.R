#!/usr/bin/env Rscript

# analyze_mirna_isoforms.R
#
# Purpose: Analyze miRNA target site density across Full 3' UTR Isoforms.
# Contrasts standard spliced UTRs against specific UFM1-dependent Retained Intron (RI) isoforms.
#
# Groups:
# A: Single-Exon 3' UTRs (Background)
# B: Multi-Exon 3' UTRs (Spliced Control)
# C: UFM1-Dependent RI Isoforms (Target) - UTR + RI
# D: UFM1-Independent RI Isoforms (Control RI) - UTR + RI
# E: Genomic Control (Intron-Only)

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
  library(tidyverse)
  library(optparse)
  library(ggpubr)
  library(rtracklayer)
})

# --- Argument Parsing ---
option_list <- list(
  make_option(c("--gtf"), type="character", help="Path to GTF file"),
  make_option(c("--events"), type="character", help="Path to UFM1 events RDS file"),
  make_option(c("--fasta"), type="character", help="Path to Genome FASTA file"),
  make_option(c("--species"), type="character", default="human", help="Species (human or mouse)"),
  make_option(c("--seeds"), type="character", default=NULL, help="Path to miRNA seeds TSV (optional, will download if missing)"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Ensure output directory exists
if(!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

# --- Helper Functions ---

# Function to load or download seeds (robust implementation)
get_mirna_seeds <- function(seed_file, species) {
  # TaxID Map
  tax_map <- list(
    "human" = 9606,
    "mouse" = 10090,
    "arabidopsis" = 3702
  )
  
  target_taxid <- tax_map[[tolower(species)]]
  if(is.null(target_taxid)) stop(paste("Unsupported species:", species))
  
  if (is.null(seed_file) || !file.exists(seed_file)) {
    message(paste("No seed file provided. Downloading conserved seeds for", species, "(TaxID:", target_taxid, ")..."))
    
    # Download Family Info
    url <- "https://www.targetscan.org/vert_80/vert_80_data_download/miR_Family_Info.txt.zip"
    temp <- tempfile()
    download.file(url, temp, quiet = TRUE)
    
    # Use read.table which handles the format reliably (dots for spaces)
    miRNA_families <- read.table(unz(temp, "miR_Family_Info.txt"), header=TRUE, sep="\t", fill=TRUE, quote="")
    unlink(temp)
    
    seeds <- miRNA_families %>%
      filter(Species.ID == target_taxid) %>%   # Dynamic TaxID
      filter(Family.Conservation. >= 2) %>%    # Conserved families only
      select(
        Family = miR.family, 
        Seed_7mer_m8 = Seed.m8
      ) %>%
      distinct(Seed_7mer_m8, .keep_all = TRUE) %>%
      rename(Seed = Seed_7mer_m8, Name = Family) %>%
      mutate(Seed = gsub("U", "T", Seed))
    
    return(seeds)
  } else {
    return(read_tsv(seed_file, col_types = cols()))
  }
}

# --- 1. Load Data ---

message("Loading GTF...")
txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")

message("Loading UFM1 Events...")
events <- readRDS(opt$events)

# Filter for RI events
ri_events <- events[events$EventType == "RI"]

# Ensure numeric columns (using _num versions if available or converting)
# Checking available columns
cols <- colnames(mcols(ri_events))

if("FDR_num" %in% cols) {
  ri_events$padj <- ri_events$FDR_num
} else {
  ri_events$padj <- as.numeric(ri_events$FDR)
}

if("dPSI_num" %in% cols) {
  ri_events$dPSI <- ri_events$dPSI_num
} else {
  # IncLevelDifference is usually dPSI
  ri_events$dPSI <- as.numeric(ri_events$IncLevelDifference)
}

# Classify Dependent vs Independent using pre-calculated Group column
dep_ri <- ri_events[ri_events$Group == "UFM1_dependent"]
indep_ri <- ri_events[ri_events$Group == "UFM1_independent"]

message(paste("Dependent RI:", length(dep_ri)))
message(paste("Independent RI:", length(indep_ri)))

# --- 2. Extract 3'UTRs grouped by Transcript ---

message("Extracting 3'UTRs...")
# primary key is transcript_id
utr3_by_tx <- threeUTRsByTranscript(txdb, use.names=TRUE)

# --- 3. Construct Groups ---

# Align seqlevels styles if needed
tryCatch({
  seqlevelsStyle(dep_ri) <- seqlevelsStyle(txdb)[1]
  seqlevelsStyle(indep_ri) <- seqlevelsStyle(txdb)[1]
  seqlevelsStyle(ri_events) <- seqlevelsStyle(txdb)[1]
}, error=function(e){
  message("Warning: Seqlevels alignment failed. Proceeding anyway but overlaps might fail.")
})

# Group A: Single-Exon 3' UTRs (Background)
# length(gr) == 1 means the 3'UTR has no introns. 
# Note: This checks for introns WITHIN the 3'UTR. 
# "Standard 3' UTRs that naturally have NO introns."
utr3_counts <- elementNROWS(utr3_by_tx)
group_a_ids <- names(utr3_by_tx)[utr3_counts == 1]
group_a_gr <- utr3_by_tx[group_a_ids]
message(paste("Group A (Single-Exon UTRs):", length(group_a_gr)))

# Group B: Multi-Exon 3' UTRs (Spliced Control)
# length(gr) > 1 means the 3'UTR is spliced.
group_b_ids <- names(utr3_by_tx)[utr3_counts > 1]
group_b_gr <- utr3_by_tx[group_b_ids]
message(paste("Group B (Multi-Exon Spliced UTRs):", length(group_b_gr)))

# Helper to construct RI Isoforms
construct_ri_isoform <- function(ri_events_gr, all_utrs_grl) {
  # Align seqlevels explicitly inside function just in case
  # But we did it globally above for the inputs.
  
  hits <- findOverlaps(ri_events_gr, all_utrs_grl)
  
  ri_isoforms <- list()
  
  # For each overlap
  for (i in seq_along(hits)) {
    ri_idx <- queryHits(hits)[i]
    tx_idx <- subjectHits(hits)[i]
    
    ri_geo <- ri_events_gr[ri_idx]
    tx_id <- names(all_utrs_grl)[tx_idx]
    original_utr <- all_utrs_grl[[tx_idx]]
    
    # Union the Intron (RI) with the UTR Exons
    # This effectively "retains" the intron
    isoform_parts <- c(original_utr, granges(ri_geo))
    
    # Merge overlapping/adjacent intervals (Exon-Intron-Exon becomes one block)
    # Explicitly use GenomicRanges::reduce to avoid purrr conflict
    isoform_merged <- GenomicRanges::reduce(isoform_parts)
    
    ri_isoforms[[tx_id]] <- isoform_merged
  }
  
  if(length(ri_isoforms) == 0) return(GRangesList())
  return(GRangesList(ri_isoforms))
}

# Group C: UFM1-Dependent RI Isoforms
message("Constructing Group C (Dependent RI Isoforms)...")
group_c_gr <- construct_ri_isoform(dep_ri, utr3_by_tx)
message(paste("Group C created:", length(group_c_gr)))

# Group D: UFM1-Independent RI Isoforms
message("Constructing Group D (Independent RI Isoforms)...")
group_d_gr <- construct_ri_isoform(indep_ri, utr3_by_tx)
message(paste("Group D created:", length(group_d_gr)))

# Group E (Intron Only, UTR-associated)
overlaps <- findOverlaps(dep_ri, utr3_by_tx)
group_e_gr <- dep_ri[unique(queryHits(overlaps))]
message(paste("Group E (Intron Only, UTR-associated):", length(group_e_gr)))

# Group F: Constitutive 3'UTR Introns (Control)
# Identical logic to Metagene Analysis
message("Constructing Group F (Constitutive 3'UTR Introns)...")
utr3_ranges_flat <- unlist(range(utr3_by_tx), use.names=TRUE)
all_introns_grl <- intronsByTranscript(txdb, use.names=TRUE)
all_introns_flat <- unlist(all_introns_grl, use.names=TRUE)
group_f_gr <- subsetByOverlaps(all_introns_flat, utr3_ranges_flat, type="within")
message(paste("Group F (Constitutive 3'UTR Introns):", length(group_f_gr)))

# --- 4. Sequence Extraction ---

message("Loading Genome...")
genome <- FaFile(opt$fasta)

get_grl_seqs <- function(grl, name) {
  message(paste("Extracting sequences for", name, "..."))
  # extractTranscriptSeqs handles splicing (joining exons)
  seqs <- extractTranscriptSeqs(genome, grl)
  return(seqs)
}

seqs_a <- get_grl_seqs(group_a_gr, "Group A")
seqs_b <- get_grl_seqs(group_b_gr, "Group B")
seqs_c <- get_grl_seqs(group_c_gr, "Group C")
seqs_d <- get_grl_seqs(group_d_gr, "Group D")

message("Extracting sequences for Group E...")
seqs_e <- getSeq(genome, group_e_gr)
if(!is.null(mcols(group_e_gr)$event_id)) {
  names(seqs_e) <- mcols(group_e_gr)$event_id
} else {
  names(seqs_e) <- paste0("Intron_Dep_", seq_along(seqs_e))
}

message("Extracting sequences for Group F...")
seqs_f <- getSeq(genome, group_f_gr)
names(seqs_f) <- paste0("Intron_Const_", seq_along(seqs_f))

message("Extracting sequences for Group E...")
# For simple GRanges (introns), use getSeq
seqs_e <- getSeq(genome, group_e_gr)
# Use generated names or event ID if available
if(!is.null(mcols(group_e_gr)$event_id)) {
  names(seqs_e) <- mcols(group_e_gr)$event_id
} else {
  names(seqs_e) <- paste0("Intron_", seq_along(seqs_e))
}

# --- 5. miRNA Scanning ---

seeds <- get_mirna_seeds(opt$seeds, opt$species)
message(paste("Loaded", nrow(seeds), "miRNA seeds."))

scan_density <- function(sequences, seeds_df, group_name) {
  if (length(sequences) == 0) {
    message(paste("Warning: No sequences for group", group_name))
    return(data.frame(ID=character(), Length_bp=numeric(), Site_Count=numeric(), Group=character(), Sites_Per_KB=numeric()))
  }

  total_kb <- sum(width(sequences)) / 1000
  total_sites <- 0
  
  # Iterate over seeds
  # Optimization: Combine seeds into one pattern or iterate? 
  # Iterating 108 seeds over thousands of seqs is usually fine in R.
  for (s in seeds_df$Seed) {
    # countPattern is fast
    counts <- vcountPattern(s, sequences)
    total_sites <- total_sites + sum(counts)
  }
  
  # Calculate per-sequence density for distribution
  # For boxplot we need per-sequence sites
  
  # Initialize count vector
  per_seq_counts <- numeric(length(sequences))
  
  for (s in seeds_df$Seed) {
    # Only counting matches on the Forward strand of the RNA sequence provided
    # (Sequences extracted are already 5'->3' of the transcript)
    per_seq_counts <- per_seq_counts + vcountPattern(s, sequences)
  }
  
  df <- data.frame(
    ID =names(sequences),
    Length_bp = width(sequences),
    Site_Count = per_seq_counts,
    Group = group_name
  )
  
  df$Sites_Per_KB <- (df$Site_Count / df$Length_bp) * 1000
  return(df)
}

message("Scanning Group A...")
res_a <- scan_density(seqs_a, seeds, "A: Single-Exon 3'UTR")
message("Scanning Group B...")
res_b <- scan_density(seqs_b, seeds, "B: Multi-Exon 3'UTR (Spliced)")
message("Scanning Group C...")
res_c <- scan_density(seqs_c, seeds, "C: Dep RI Isoform")
message("Scanning Group D...")
res_d <- scan_density(seqs_d, seeds, "D: Indep RI Isoform")
message("Scanning Group E...")
res_e <- scan_density(seqs_e, seeds, "E: Intron Only")
message("Scanning Group F...")
res_f <- scan_density(seqs_f, seeds, "F: Constitutive Intron")

all_results <- bind_rows(res_a, res_b, res_c, res_d, res_e, res_f)

# --- 6. Stats & Visualization ---

message("Performing Statistical Analysis...")

# Summary Table
summary_stats <- all_results %>%
  group_by(Group) %>%
  summarise(
    N = n(),
    Mean_Density = mean(Sites_Per_KB, na.rm=TRUE),
    Median_Density = median(Sites_Per_KB, na.rm=TRUE),
    SD = sd(Sites_Per_KB, na.rm=TRUE)
  )

print(summary_stats)
write_tsv(summary_stats, file.path(opt$outdir, "isoform_mirna_stats_summary.tsv"))

# Pairwise Comparisons of Interest
# C vs B (Dilution)
# E vs A (GC/Intron check)

comparisons <- list(
  c("C: Dep RI Isoform", "B: Multi-Exon 3'UTR (Spliced)"),
  c("E: Intron Only", "A: Single-Exon 3'UTR"),
  c("C: Dep RI Isoform", "D: Indep RI Isoform")
)

# Plot
p <- ggplot(all_results, aes(x=Group, y=Sites_Per_KB, fill=Group)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 20)) + # Focus on reasonable range, ignore extreme outliers
  labs(title="miRNA Target Density in 3'UTR Isoforms",
       y="Sites per KB",
       x="") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", p.adjust.method = "BH")

# Plot
p <- ggplot(all_results, aes(x=Group, y=Sites_Per_KB, fill=Group)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 20)) + # Focus on reasonable range, ignore extreme outliers
  labs(title="miRNA Target Density in 3'UTR Isoforms",
       y="Sites per KB",
       x="") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(ref.group = "C: Dep RI Isoform", label = "p.signif", method = "wilcox.test", hide.ns = FALSE)

ggsave(file.path(opt$outdir, "miRNA_isoform_density_boxplot.pdf"), p, width=8, height=6)
message("Boxplot saved.")

# --- 7. CDF Plot (Cumulative Distribution Function) ---
# Specific request: Cons UTR M (Group B), Group C, Group D, Group F

cdf_groups <- c("B: Multi-Exon 3'UTR (Spliced)", "C: Dep RI Isoform", "D: Indep RI Isoform", "F: Constitutive Intron")
cdf_data <- all_results %>% filter(Group %in% cdf_groups)

# Clean Names for Plot Legend
cdf_data$Group_Label <- recode(cdf_data$Group,
  "B: Multi-Exon 3'UTR (Spliced)" = "Spliced UTR (Control)",
  "C: Dep RI Isoform" = "UFM1-Dependent RI",
  "D: Indep RI Isoform" = "UFM1-Independent RI",
  "F: Constitutive Intron" = "Constitutive Intron"
)

p_cdf <- ggplot(cdf_data, aes(x=Sites_Per_KB, color=Group_Label)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  coord_cartesian(xlim = c(0, 50)) + # Focus on 0-50 sites/kb
  labs(title="miRNA Target Density CDF",
       y="Cumulative Fraction",
       x="miRNA Sites per KB",
       color="Isoform Type") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right")

ggsave(file.path(opt$outdir, "miRNA_density_CDF.pdf"), p_cdf, width=8, height=6)
message("CDF Plot saved.")

message("Done.")
