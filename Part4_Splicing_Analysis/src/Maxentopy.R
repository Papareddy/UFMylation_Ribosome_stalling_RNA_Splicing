# ==============================================================================
# SCRIPT: Asymmetric Splice Site Sequence Logos
# FEATURES: Separate Exon/Intron window sizes, Auto-Strand, Chromosome Fix
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
  library(ggseqlogo)
  library(gridExtra)
  library(tidyverse)
})

# ==============================================================================
# 1. CONFIGURATION (Adjust variables here)
# ==============================================================================
# Paths
FASTA_PATH <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/pcg_gencode.v45.annotation.gtf"
RDS_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/results/human/nucleus/step01_data_prep/UFM1_events_rich.rds"

# --- ASYMMETRIC VARIABLES ---
WIN_INTRON <- 20  # Bases to show INTO the Intron (e.g., for PPT)
WIN_EXON   <- 20   # Bases to show INTO the Exon (e.g., for consensus)

# ==============================================================================
# 2. DATA LOADING & PREP
# ==============================================================================
message("--- Loading Data ---")
txdb <- makeTxDbFromGFF(GTF_PATH, format="gtf")
events <- readRDS(RDS_PATH)
genome_fa <- FaFile(FASTA_PATH)

# Define Groups
ri_events <- events[events$EventType == "RI", ]
dep_ri_raw    <- ri_events[ri_events$Group == "UFM1_dependent", ]
indep_ri_raw  <- ri_events[ri_events$Group == "UFM1_independent", ]

# Define Control
all_introns_grl  <- intronsByTranscript(txdb, use.names=TRUE)
all_introns_flat <- unlist(all_introns_grl, use.names=TRUE)
utr3_grl     <- threeUTRsByTranscript(txdb, use.names=TRUE)
utr3_ranges  <- unlist(range(utr3_grl), use.names=TRUE)
utr3_introns <- subsetByOverlaps(all_introns_flat, utr3_ranges, type="within")

targets <- c(dep_ri_raw, indep_ri_raw)
hits <- findOverlaps(utr3_introns, targets)
control_introns_raw <- if(length(hits) > 0) utr3_introns[-queryHits(hits)] else utr3_introns

# Convert to Unique Loci
message("--- Collapsing to Unique Loci ---")
dep_ri_uniq   <- unique(granges(dep_ri_raw))
indep_ri_uniq <- unique(granges(indep_ri_raw))
control_uniq  <- unique(granges(control_introns_raw))

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

# --- A. Chromosome Name Fixer ---
ensure_seqlevels_match <- function(gr, fasta_file) {
  fasta_levels <- seqlevels(seqinfo(fasta_file))
  has_chr_fasta <- any(grepl("^chr", fasta_levels))
  has_chr_gr    <- any(grepl("^chr", seqlevels(gr)))
  
  if (!has_chr_fasta && has_chr_gr) {
    seqlevelsStyle(gr) <- "NCBI"
  } else if (has_chr_fasta && !has_chr_gr) {
    seqlevelsStyle(gr) <- "UCSC"
  }
  return(gr)
}

# --- B. Strand Fixer ---
fix_strand_based_on_sequence <- function(gr, fasta_file) {
  gr <- ensure_seqlevels_match(gr, fasta_file)
  gr_check <- gr
  strand(gr_check) <- "+"
  
  seq_start <- getSeq(fasta_file, GRanges(seqnames(gr_check), IRanges(start(gr_check), start(gr_check)+1)))
  seq_end   <- getSeq(fasta_file, GRanges(seqnames(gr_check), IRanges(end(gr_check)-1, end(gr_check))))
  s_start <- as.character(seq_start)
  s_end   <- as.character(seq_end)
  
  new_strand <- rep("*", length(gr))
  new_strand[s_start == "GT" & s_end == "AG"] <- "+"
  new_strand[s_start == "CT" & s_end == "AC"] <- "-"
  
  strand(gr) <- new_strand
  gr_fixed <- gr[strand(gr) != "*"]
  return(gr_fixed)
}

# --- C. Sequence Extractor (ASYMMETRIC LOGIC) ---
get_site_sequences <- function(gr, fasta_file, type="5p", w_exon=3, w_intron=15) {
  
  gr <- ensure_seqlevels_match(gr, fasta_file)
  is_plus <- as.character(strand(gr)) == "+"
  
  if(type == "5p") {
    # 5' Donor (Start of Intron)
    # Structure: [Exon] | [Intron]
    center <- ifelse(is_plus, start(gr), end(gr))
    
    # Plus: [Center - w_exon, Center + w_intron - 1]
    # Note: start(gr) is the 1st base of Intron. 
    # So center-1 is last base of Exon.
    gr_win <- GRanges(seqnames(gr), IRanges(center - w_exon, center + w_intron - 1), strand=strand(gr))
    
  } else {
    # 3' Acceptor (End of Intron)
    # Structure: [Intron] | [Exon]
    center <- ifelse(is_plus, end(gr), start(gr))
    
    # Plus: [Center - w_intron + 1, Center + w_exon]
    # Note: end(gr) is last base of Intron.
    # So center+1 is 1st base of Exon.
    gr_win <- GRanges(seqnames(gr), IRanges(center - w_intron + 1, center + w_exon), strand=strand(gr))
  }
  
  gr_win <- trim(gr_win)
  gr_win <- gr_win[width(gr_win) == (w_exon + w_intron)]
  if(length(gr_win) == 0) return(NULL)
  
  return(as.character(getSeq(fasta_file, gr_win)))
}

# ==============================================================================
# 4. EXECUTION
# ==============================================================================
message("--- Auto-Correcting Strands ---")
dep_fixed   <- fix_strand_based_on_sequence(dep_ri_uniq, genome_fa)
indep_fixed <- fix_strand_based_on_sequence(indep_ri_uniq, genome_fa)
ctrl_fixed  <- fix_strand_based_on_sequence(control_uniq, genome_fa)

message(paste0("Extracting: Exon=", WIN_EXON, "bp / Intron=", WIN_INTRON, "bp..."))

# 5' Donor
s5_dep   <- get_site_sequences(dep_fixed, genome_fa, "5p", WIN_EXON, WIN_INTRON)
s5_indep <- get_site_sequences(indep_fixed, genome_fa, "5p", WIN_EXON, WIN_INTRON)
s5_ctrl  <- get_site_sequences(ctrl_fixed, genome_fa, "5p", WIN_EXON, WIN_INTRON)

# 3' Acceptor
s3_dep   <- get_site_sequences(dep_fixed, genome_fa, "3p", WIN_EXON, WIN_INTRON)
s3_indep <- get_site_sequences(indep_fixed, genome_fa, "3p", WIN_EXON, WIN_INTRON)
s3_ctrl  <- get_site_sequences(ctrl_fixed, genome_fa, "3p", WIN_EXON, WIN_INTRON)

# ==============================================================================
# 5. PLOTTING (ASYMMETRIC)
# ==============================================================================
message("--- Generating Logos ---")

# --- 5' PLOT (Exon | Intron) ---
# Line position is after the Exon bases
line_pos_5p <- WIN_EXON + 0.5 

p5_list <- list("Control"=s5_ctrl, "UFM1 Indep"=s5_indep, "UFM1 Dep"=s5_dep)
p_5p <- ggseqlogo(p5_list, ncol=1) + 
  theme_classic() +
  ggtitle(paste0("5' Donor Site (Exon ", WIN_EXON, "bp | Intron ", WIN_INTRON, "bp)")) +
  geom_vline(xintercept = line_pos_5p, linetype="dashed", color="black", size=0.8) +
  annotate("text", x=WIN_EXON/2 + 0.5, y=2, label="Exon", color="black", fontface="italic", size=3) +
  annotate("text", x=WIN_EXON + WIN_INTRON/2, y=2, label="Intron", color="black", fontface="italic", size=3) +
  theme(axis.text.x = element_blank()) 

# --- 3' PLOT (Intron | Exon) ---
# Line position is after the Intron bases
line_pos_3p <- WIN_INTRON + 0.5 

p3_list <- list("Control"=s3_ctrl, "UFM1 Indep"=s3_indep, "UFM1 Dep"=s3_dep)
p_3p <- ggseqlogo(p3_list, ncol=1) + 
  theme_classic() +
  ggtitle(paste0("3' Acceptor Site (Intron ", WIN_INTRON, "bp | Exon ", WIN_EXON, "bp)")) +
  geom_vline(xintercept = line_pos_3p, linetype="dashed", color="black", size=0.8) +
  annotate("text", x=WIN_INTRON/2, y=2, label="Intron", color="black", fontface="italic", size=3) +
  annotate("text", x=WIN_INTRON + WIN_EXON/2 + 0.5, y=2, label="Exon", color="black", fontface="italic", size=3) +
  theme(axis.text.x = element_blank())

# COMBINE
grid.arrange(p_5p, p_3p, ncol=2)

message("Done.")

