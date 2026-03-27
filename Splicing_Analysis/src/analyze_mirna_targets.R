
# ==============================================================================
# Analyze miRNA Targeting: UFM1 RIs vs Controls
# ==============================================================================

suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(GenomicRanges)
    library(rtracklayer)
    library(Biostrings)
    library(Rsamtools) # Use Rsamtools for FASTA access
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
  make_option(c("--fasta"), type="character", default=NULL, 
              help="Path to Genomic FASTA", metavar="FASTA"),
  make_option(c("--events"), type="character", default=NULL, 
              help="Path to UFM1_events_rich.rds", metavar="RDS"),
  make_option(c("--seeds"), type="character", default=NULL, 
              help="Path to miRNA seeds TSV", metavar="SEEDS"),
  make_option(c("--outdir"), type="character", default="results/mirna_analysis", 
              help="Output directory", metavar="OUTDIR")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gtf) || is.null(opt$events) || is.null(opt$fasta)) {
  print_help(opt_parser)
  stop("Missing required arguments (GTF, Events, FASTA).", call.=FALSE)
}

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# ==========================================
# 2. LOAD DATA & SEQUENCES
# ==========================================

message("Loading Data...")
ufm_events <- readRDS(opt$events)

if (!is.null(opt$seeds)) {
    message("Loading provided seed file: ", opt$seeds)
    seeds <- read.table(opt$seeds, header=TRUE, sep="\t", stringsAsFactors=FALSE)
} else {
    message("No seed file provided. Downloading robust list from TargetScan 8.0...")
    
    # 1. Download Family Info
    url <- "https://www.targetscan.org/vert_80/vert_80_data_download/miR_Family_Info.txt.zip"
    temp <- tempfile()
    download.file(url, temp)
    miRNA_families <- read.table(unz(temp, "miR_Family_Info.txt"), header=TRUE, sep="\t", fill=TRUE, quote="")
    unlink(temp)
    
    # Check column names if they differ
    # TargetScan headers can vary.
    # Expected: "Species.ID", "Family.Conservation.", "miR.family", "Seed.m8", "MiRBase.ID"
    # User's snippet used "Species.ID", "Family.Conservation?", "miR.family", "Seed.m8"
    
    # 2. Filter for Human & Conserved
    seeds <- miRNA_families %>%
      filter(Species.ID == 9606) %>%           # Human Only
      filter(Family.Conservation. >= 2) %>%    # Conserved families only
      select(
        Family = miR.family, 
        Seed_7mer_m8 = Seed.m8, 
        MiRBase_ID = MiRBase.ID
      ) %>%
      distinct(Seed_7mer_m8, .keep_all = TRUE)   # Remove duplicates
      
    message("Generated Robust Seed List: ", nrow(seeds), " unique conserved seeds.")
}

message("Loading Genomic Fasta: ", opt$fasta)
genome_fa <- FaFile(opt$fasta)

# Filter RIs overlapping 3'UTRs
# We need to replicate the 3'UTR overlap logic to be precise
txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")
utrs_by_tx <- threeUTRsByTranscript(txdb, use.names=TRUE)

get_3utr_overlapping_ris <- function(events) {
    events <- subset(events, EventType == "RI")
    # Align styles
    tryCatch({seqlevelsStyle(events) <- seqlevelsStyle(utrs_by_tx)[1]}, error=function(e){})
    
    hits <- findOverlaps(events, utrs_by_tx)
    subset(events, unique(queryHits(hits))) # Return the EVENTS that overlap any 3'UTR
}

# Define Groups
ufm_dep_ri   <- get_3utr_overlapping_ris(subset(ufm_events, Group == "UFM1_dependent"))
ufm_indep_ri <- get_3utr_overlapping_ris(subset(ufm_events, Group == "UFM1_independent"))

message("Group Sizes (RI in 3'UTR):")
message("  - Dependent: ", length(ufm_dep_ri))
message("  - Independent: ", length(ufm_indep_ri))

# Controls
# 1. Random Constitutive Introns
message("Generating Controls...")
all_introns <- unique(unlist(intronsByTranscript(txdb)))
# Remove overlaps with our RIs
overlaps <- countOverlaps(all_introns, c(ufm_dep_ri, ufm_indep_ri))
const_introns <- all_introns[overlaps == 0]
# Sample 5000
set.seed(123)
ctrl_introns <- sample(const_introns, min(5000, length(const_introns)))

# 2. Random Constitutive 3'UTRs
all_utrs <- unlist(utrs_by_tx)
ctrl_utrs <- sample(all_utrs, min(5000, length(all_utrs)))

# Extract Sequences
message("Extracting Sequences...")
get_seqs <- function(gr) {
    # Ensure chromosomes present in genome
    # FaFile requires a check or try-catch usually, but getSeq handles it if seqlevels match
    # Seqlevels in FaFile are usually "1", "2" or "chr1", "chr2"
    fa_levels <- seqlevels(seqinfo(genome_fa))
    
    # Align seqlevels
    tryCatch({seqlevelsStyle(gr) <- seqlevelsStyle(fa_levels)[1]}, error=function(e){})
    
    # Filter for valid chroms
    gr <- gr[seqnames(gr) %in% fa_levels]
    
    # Extract
    getSeq(genome_fa, gr)
}

seqs_list <- list(
    "UFM1-Dep RI" = get_seqs(ufm_dep_ri),
    "UFM1-Indep RI" = get_seqs(ufm_indep_ri),
    "Const Introns" = get_seqs(ctrl_introns),
    "Const 3'UTRs" = get_seqs(ctrl_utrs)
)

# ==========================================
# 3. SCAN MIRNA SEEDS
# ==========================================

scan_seeds <- function(sequences, seeds_df) {
    # Simple scan: count occurrences of each seed in each sequence
    # Returns density (sites/kb)
    
    total_sites <- numeric(length(sequences))
    
    for (s in seeds_df$Seed_7mer_m8) {
        # DNA-RNA conversion: Seed is usually given as RNA (U), Genome is DNA (T)
        dna_seed <- gsub("U", "T", s)
        counts <- vcountPattern(dna_seed, sequences)
        total_sites <- total_sites + counts
    }
    
    # Calculate Density (sites per kb)
    # Avoid division by zero
    lens <- width(sequences)
    norm_density <- (total_sites / lens) * 1000
    
    return(norm_density)
}

message("Scanning for ", nrow(seeds), " miRNA seeds...")
results <- data.frame()

for (grp in names(seqs_list)) {
    dens <- scan_seeds(seqs_list[[grp]], seeds)
    tmp <- data.frame(
        Group = grp,
        Density = dens
    )
    results <- rbind(results, tmp)
}

# ==========================================
# 4. PLOTTING & STATS
# ==========================================

# Order
results$Group <- factor(results$Group, levels = c("Const Introns", "UFM1-Indep RI", "UFM1-Dep RI", "Const 3'UTRs"))

# Stats
message("\nDensity Stats (Sites/kb):")
stats <- results %>% group_by(Group) %>% summarise(
    Mean = mean(Density),
    Median = median(Density)
)
print(stats)

# Wilcox Test
pval <- wilcox.test(
    results$Density[results$Group == "UFM1-Dep RI"],
    results$Density[results$Group == "UFM1-Indep RI"]
)$p.value
message("\nWilcoxon P-value (Dep vs Indep): ", pval)

# Plot CDF
p <- ggplot(results, aes(x = Density, color = Group)) +
  stat_ecdf(geom = "step", size = 1.2) +
  scale_color_manual(values = c(
      "Const Introns" = "grey70",
      "UFM1-Indep RI" = "#64B5F6", # Blue
      "UFM1-Dep RI" = "#E57373",   # Red
      "Const 3'UTRs" = "black"
  )) +
  coord_cartesian(xlim = c(0, 20)) + # Zoom in meaningful range
  labs(
    title = "miRNA Target Site Density (CDF)",
    subtitle = paste0("Conserved Seeds (Top ", nrow(seeds), ")"),
    x = "Site Density (Sites per kb)",
    y = "Cumulative Fraction"
  ) +
  theme_bw(base_size = 14) +
  annotate("text", x = 15, y = 0.2, label = paste0("P(Dep vs Indep) = ", format.pval(pval, digits=3)))

out_plot <- file.path(opt$outdir, "miRNA_density_plot.pdf")
ggsave(out_plot, p, width = 7, height = 5)
message("Saved plot to: ", out_plot)
