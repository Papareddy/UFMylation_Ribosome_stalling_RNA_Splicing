# ==============================================================================
# SCRIPT: UFM1 Intron Translation Analysis (Full CDS Stitching)
# PURPOSE: Determine if retained introns introduce Stop Codons (PTCs)
# METHOD:  Stitches Intron into Full CDS -> Translates from authentic ATG
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(Biostrings)
  library(Rsamtools)
  library(tidyverse)
  library(BSgenome.Hsapiens.UCSC.hg38) # Standard Human Genome
  # NOTE: If using Mouse/Arabidopsis, load the corresponding BSgenome package
  # e.g., library(BSgenome.Mmusculus.UCSC.mm10)
})

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================
# Use the paths you provided
FASTA_PATH <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/pcg_gencode.v45.annotation.gtf"
RDS_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/results/human/nucleus/step01_data_prep/UFM1_events_rich.rds"

# Output settings
OUTPUT_FILE <- "UFM1_Dependent_Intron_Translation_Results.csv"

# ==============================================================================
# 2. DATA LOADING
# ==============================================================================
message("--- Loading Data ---")

# A. Load Annotations
txdb <- makeTxDbFromGFF(GTF_PATH, format="gtf")
# Pre-load CDS grouped by transcript (for speed)
cds_by_tx <- cdsBy(txdb, by="tx", use.names=TRUE)

# B. Load Events
events <- readRDS(RDS_PATH)

# C. Filter for UFM1 Dependent Introns (RI)
# We strictly check for 'RI' (Retained Intron) events
ri_events <- events[events$EventType == "RI", ]
dep_ri    <- ri_events[ri_events$Group == "UFM1_dependent", ]

message(paste("Found", length(dep_ri), "UFM1-Dependent Retained Introns to analyze."))

# ==============================================================================
# 3. ANALYSIS FUNCTION: STITCH & TRANSLATE
# ==============================================================================

translate_retention <- function(intron_gr, cds_list, genome_seq) {
  
  # 1. Find the Parent Transcript
  # We look for transcripts that overlap this intron
  hits <- findOverlaps(intron_gr, cds_list)
  
  if(length(hits) == 0) return(NULL) # No coding transcript found for this intron location
  
  # Pick the first valid transcript (canonical/principal is usually 1st in list)
  tx_name <- names(cds_list)[subjectHits(hits)[1]]
  parent_cds <- cds_list[[tx_name]]
  
  # 2. Create "Mutant" CDS (Simulate Retention)
  # We simply add the intron range to the existing CDS ranges
  # Note: This handles strands automatically because GRanges objects are strand-aware
  mutant_gr <- c(parent_cds, intron_gr)
  
  # 3. Merge and Sort
  # reduce() fuses the intron into the exons, creating one continuous block
  # sort() ensures they are in correct genomic order
  mutant_merged <- reduce(sort(mutant_gr))
  
  # 4. Extract DNA Sequence
  # This gets the sequence 5' -> 3' respecting the strand
  dna <- getSeq(genome_seq, mutant_merged)
  tx_seq <- unlist(dna) # Collapse into single string
  
  # 5. Translate
  # if.fuzzy.codon="solve" handles end-of-seq fragments gracefully
  aa_seq <- translate(tx_seq, if.fuzzy.codon="solve")
  protein_str <- as.character(aa_seq)
  
  # 6. Analyze Stop Codons
  # Find the first "*" (Stop)
  first_stop <- regexpr("\\*", protein_str)[1]
  
  # Calculate where the intron actually started in protein space
  # (To see if the stop is INSIDE the intron or AFTER)
  # Strategy: Measure length of upstream exons
  
  # Get upstream CDS only
  if (as.character(strand(intron_gr)) == "+") {
    upstream_cds <- parent_cds[end(parent_cds) < start(intron_gr)]
  } else {
    upstream_cds <- parent_cds[start(parent_cds) > end(intron_gr)]
  }
  
  upstream_dna_len <- sum(width(upstream_cds))
  upstream_aa_len  <- floor(upstream_dna_len / 3)
  
  # Return detailed stats
  return(data.frame(
    Intron_UID = paste0(seqnames(intron_gr), ":", start(intron_gr), "-", end(intron_gr)),
    Gene_ID = tx_name,
    Strand = as.character(strand(intron_gr)),
    Protein_Length_WT = sum(width(parent_cds))/3,
    Protein_Length_Mutant = nchar(protein_str),
    Intron_Start_AA_Index = upstream_aa_len,
    First_Stop_AA_Index = first_stop,
    # Status Logic:
    # If Stop is -1 -> No Stop (Readthrough)
    # If Stop <= Intron_Start -> Stop was already there (Annotated PTC? Rare)
    # If Stop > Intron_Start -> Stop caused by Retention
    Result = case_when(
      first_stop == -1 ~ "No_Stop_Codon",
      first_stop <= upstream_aa_len ~ "Pre_Existing_Stop", 
      first_stop > upstream_aa_len & first_stop < (upstream_aa_len + 30) ~ "Immediate_PTC",
      TRUE ~ "Downstream_PTC"
    ),
    Snippet_At_Junction = substr(protein_str, max(1, upstream_aa_len - 5), min(nchar(protein_str), upstream_aa_len + 10))
  ))
}

# ==============================================================================
# 4. EXECUTION LOOP
# ==============================================================================
message("--- Running Translation Analysis (This may take a moment) ---")

# Load the genome FASTA into memory-efficient object
genome_fa <- FaFile(FASTA_PATH)

results_list <- list()

# Use progress bar if available, otherwise simple loop
pb <- txtProgressBar(min = 0, max = length(dep_ri), style = 3)

for(i in seq_along(dep_ri)) {
  res <- translate_retention(dep_ri[i], cds_by_tx, genome_fa)
  if(!is.null(res)) {
    results_list[[i]] <- res
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# Combine Results
final_df <- do.call(rbind, results_list)

# ==============================================================================
# 5. SUMMARY & EXPORT
# ==============================================================================

if (!is.null(final_df)) {
  message("\n--- Analysis Complete ---")
  print(table(final_df$Result))
  
  # Show a few examples of "Immediate_PTC"
  message("\n--- Examples of Immediate Stops (The 'TAG' Trap) ---")
  print(head(final_df[final_df$Result == "Immediate_PTC", c("Gene_ID", "Intron_Start_AA_Index", "First_Stop_AA_Index", "Snippet_At_Junction")], 5))
  
  # Save
  # write.csv(final_df, OUTPUT_FILE, row.names = FALSE)
  # message(paste("Results saved to:", OUTPUT_FILE))
  
} else {
  message("No valid coding transcripts found for overlapping introns.")
}