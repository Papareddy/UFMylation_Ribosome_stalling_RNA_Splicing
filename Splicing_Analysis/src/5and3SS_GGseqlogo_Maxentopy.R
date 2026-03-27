# ==============================================================================
# SCRIPT: Asymmetric Splice Site Sequence Logos (Two Controls + Method Toggle)
# GROUPS: 1. CDS Introns, 2. 3'UTR Introns, 3. UFM1 Dep, 4. UFM1 Indep
# FEATURES: Auto-Strand, Custom Windows, Bits/Probability Toggle
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
# 1. CONFIGURATION
# ==============================================================================
FASTA_PATH <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/pcg_gencode.v45.annotation.gtf"
RDS_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/results/human/nucleus/step01_data_prep/UFM1_events_rich.rds"


FASTA_PATH <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Mus_musculus.GRCm39.112.gtf"
RDS_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/results/mouse/total/step01_data_prep/UFM1_events_rich.rds"

FASTA_PATH <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/arabidopsis/Arabidopsis_thaliana_TAIR10.dna.primary_assembly.fa"
GTF_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/no_plastid_no_rRNA.Arabidopsis_thaliana.TAIR10.56.gtf"
RDS_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/results/arabidopsis/nucleus/step01_data_prep/UFM1_events_rich.rds"



# --- VISUALIZATION SETTINGS ---
# Choose: "bits" (Information Content) OR "prob" (Probability/Frequency)
LOGO_METHOD <- "probability"  

# --- 5' DONOR WINDOWS ---
D_EXON   <- 20  # Upstream Exon
D_INTRON <- 6 # Downstream Intron

# --- 3' ACCEPTOR WINDOWS ---
A_INTRON <- 20  # Upstream Intron (PPT)
A_EXON   <- 6  # Downstream Exon


# ==============================================================================
# 2. DATA LOADING & GROUP DEFINITION
# ==============================================================================
message("--- Loading Data ---")
txdb <- makeTxDbFromGFF(GTF_PATH, format="gtf")
events <- readRDS(RDS_PATH)
genome_fa <- FaFile(FASTA_PATH)

# --- A. Define UFM1 Groups ---
ri_events <- events[events$EventType == "RI", ]
dep_ri_raw    <- ri_events[ri_events$Group == "UFM1_dependent", ]
indep_ri_raw  <- ri_events[ri_events$Group == "UFM1_independent", ]

# --- B. Define Controls (CDS vs 3'UTR) ---
all_introns_grl  <- intronsByTranscript(txdb, use.names=TRUE)
all_introns_flat <- unlist(all_introns_grl, use.names=TRUE)

# 1. Get Region Ranges
cds_grl <- cdsBy(txdb, by="tx", use.names=TRUE)
cds_ranges <- unlist(range(cds_grl), use.names=TRUE) 

utr3_grl <- threeUTRsByTranscript(txdb, use.names=TRUE)
utr3_ranges <- unlist(range(utr3_grl), use.names=TRUE)

# 2. Subset Introns
cds_introns_raw <- subsetByOverlaps(all_introns_flat, cds_ranges, type="within")
utr3_introns_raw <- subsetByOverlaps(all_introns_flat, utr3_ranges, type="within")

# 3. Clean Controls (Remove overlaps with Target RIs)
targets <- c(dep_ri_raw, indep_ri_raw)

hits_cds <- findOverlaps(cds_introns_raw, targets)
ctrl_cds_clean <- if(length(hits_cds) > 0) cds_introns_raw[-queryHits(hits_cds)] else cds_introns_raw

hits_utr <- findOverlaps(utr3_introns_raw, targets)
ctrl_utr_clean <- if(length(hits_utr) > 0) utr3_introns_raw[-queryHits(hits_utr)] else utr3_introns_raw

# --- C. Unique Loci ---
message("--- Collapsing to Unique Loci ---")
dep_ri_uniq   <- unique(granges(dep_ri_raw))
indep_ri_uniq <- unique(granges(indep_ri_raw))
ctrl_cds_uniq <- unique(granges(ctrl_cds_clean))
ctrl_utr_uniq <- unique(granges(ctrl_utr_clean))

message(paste("N CDS Control:", length(ctrl_cds_uniq)))
message(paste("N 3'UTR Control:", length(ctrl_utr_uniq)))

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

# --- A. Chromosome Name Fixer ---
ensure_seqlevels_match <- function(gr, fasta_file) {
  fasta_levels <- seqlevels(seqinfo(fasta_file))
  has_chr_fasta <- any(grepl("^chr", fasta_levels))
  has_chr_gr    <- any(grepl("^chr", seqlevels(gr)))
  if (!has_chr_fasta && has_chr_gr) { seqlevelsStyle(gr) <- "NCBI" } 
  else if (has_chr_fasta && !has_chr_gr) { seqlevelsStyle(gr) <- "UCSC" }
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

# --- C. Robust Sequence Extractor ---
get_site_sequences <- function(gr, fasta_file, type="5p", w_exon, w_intron) {
  gr <- ensure_seqlevels_match(gr, fasta_file)
  target_width <- w_exon + w_intron
  final_seqs <- character()
  
  idx_plus <- which(as.character(strand(gr)) == "+")
  if(length(idx_plus) > 0) {
    gr_p <- gr[idx_plus]
    if(type == "5p") {
      gr_win <- GRanges(seqnames(gr_p), IRanges(start(gr_p) - w_exon, start(gr_p) + w_intron - 1), strand="+")
    } else {
      gr_win <- GRanges(seqnames(gr_p), IRanges(end(gr_p) - w_intron + 1, end(gr_p) + w_exon), strand="+")
    }
    gr_win <- trim(gr_win)
    gr_win <- gr_win[width(gr_win) == target_width]
    if(length(gr_win) > 0) final_seqs <- c(final_seqs, as.character(getSeq(fasta_file, gr_win)))
  }
  
  idx_minus <- which(as.character(strand(gr)) == "-")
  if(length(idx_minus) > 0) {
    gr_m <- gr[idx_minus]
    if(type == "5p") {
      gr_win <- GRanges(seqnames(gr_m), IRanges(end(gr_m) - w_intron, end(gr_m) + w_exon - 1), strand="-")
    } else {
      gr_win <- GRanges(seqnames(gr_m), IRanges(start(gr_m) - w_exon, start(gr_m) + w_intron), strand="-")
    }
    gr_win <- trim(gr_win)
    gr_win <- gr_win[width(gr_win) == target_width]
    if(length(gr_win) > 0) final_seqs <- c(final_seqs, as.character(getSeq(fasta_file, gr_win)))
  }
  return(final_seqs)
}

# ==============================================================================
# 4. EXECUTION
# ==============================================================================
message("--- Auto-Correcting Strands ---")
dep_fixed   <- fix_strand_based_on_sequence(dep_ri_uniq, genome_fa)
indep_fixed <- fix_strand_based_on_sequence(indep_ri_uniq, genome_fa)
cds_fixed   <- fix_strand_based_on_sequence(ctrl_cds_uniq, genome_fa)
utr_fixed   <- fix_strand_based_on_sequence(ctrl_utr_uniq, genome_fa)

message("--- Extracting Sequences ---")

# 5' Donor
s5_dep   <- get_site_sequences(dep_fixed, genome_fa, "5p", w_exon=D_EXON, w_intron=D_INTRON)
s5_indep <- get_site_sequences(indep_fixed, genome_fa, "5p", w_exon=D_EXON, w_intron=D_INTRON)
s5_cds   <- get_site_sequences(cds_fixed, genome_fa, "5p", w_exon=D_EXON, w_intron=D_INTRON)
s5_utr   <- get_site_sequences(utr_fixed, genome_fa, "5p", w_exon=D_EXON, w_intron=D_INTRON)

# 3' Acceptor
s3_dep   <- get_site_sequences(dep_fixed, genome_fa, "3p", w_exon=A_EXON, w_intron=A_INTRON)
s3_indep <- get_site_sequences(indep_fixed, genome_fa, "3p", w_exon=A_EXON, w_intron=A_INTRON)
s3_cds   <- get_site_sequences(cds_fixed, genome_fa, "3p", w_exon=A_EXON, w_intron=A_INTRON)
s3_utr   <- get_site_sequences(utr_fixed, genome_fa, "3p", w_exon=A_EXON, w_intron=A_INTRON)

# ==============================================================================
# 5. PLOTTING
# ==============================================================================
message(paste("--- Generating Logos (Method:", LOGO_METHOD, ") ---"))

# --- 5' PLOT ---
line_pos_5p <- D_EXON + 0.5 
p5_list <- list(
  
  "Control (3'UTR)"    = s5_utr, 
  "UFM1 Dependent"     = s5_dep, 
  "UFM1 Independent"   = s5_indep
)

p_5p <- ggseqlogo(p5_list, ncol=1, seq_type="dna",method=LOGO_METHOD) + 
  theme_classic() +
  ggtitle(paste0("5' Donor Site (Exon ", D_EXON, "bp | Intron ", D_INTRON, "bp)")) +
  geom_vline(xintercept = line_pos_5p, linetype="dashed", color="black", size=0.8) +
  # Only draw bits threshold line if using 'bits'
  {if(LOGO_METHOD=="bits") geom_hline(yintercept = 1, linetype="dashed", color="black", size=0.8)} +
  annotate("text", x=D_EXON/2 + 0.5, y=ifelse(LOGO_METHOD=="bits", 2, 1), label="Exon", color="black", fontface="italic", size=3) +
  annotate("text", x=D_EXON + D_INTRON/2, y=ifelse(LOGO_METHOD=="bits", 2, 1), label="Intron", color="black", fontface="italic", size=3) +
  theme(axis.text.x = element_blank())+my_plot_theme()

# --- 3' PLOT ---
line_pos_3p <- A_INTRON + 0.5 
p3_list <- list(
  
  "Control (3'UTR)"    = s3_utr, 
  "UFM1 Dependent"     = s3_dep, 
  "UFM1 Independent"   = s3_indep
)

p_3p <- ggseqlogo(p3_list, ncol=1, seq_type="dna",method=LOGO_METHOD) + 
  theme_classic() +
  ggtitle(paste0("3' Acceptor Site (Intron ", A_INTRON, "bp | Exon ", A_EXON, "bp)")) +
  geom_vline(xintercept = line_pos_3p, linetype="dashed", color="black", size=0.8) +
  {if(LOGO_METHOD=="bits") geom_hline(yintercept = 1, linetype="dashed", color="black", size=0.8)} +
  annotate("text", x=A_INTRON/2, y=ifelse(LOGO_METHOD=="bits", 2, 1), label="Intron", color="black", fontface="italic", size=3) +
  annotate("text", x=A_INTRON + A_EXON/2 + 0.5, y=ifelse(LOGO_METHOD=="bits", 2, 1), label="Exon", color="black", fontface="italic", size=3) +
  theme(axis.text.x = element_blank())+my_plot_theme()

# COMBINE
grid.arrange(p_5p, p_3p, ncol=2)

message("Done.")

