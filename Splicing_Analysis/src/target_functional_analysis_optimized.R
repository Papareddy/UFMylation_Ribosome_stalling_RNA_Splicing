#!/usr/bin/env Rscript
# ==============================================================================
# UNIVERSAL UFM1 ANALYSIS: HUMAN / MOUSE / ARABIDOPSIS
# OPTIM

IZED VERSION - Improved performance and caching
# ==============================================================================

# --- 1. SPECIES CONFIGURATION (CHANGE THIS!) ---
# Options: "Human", "Mouse", "Arabidopsis"
CURRENT_SPECIES <- "Arabidopsis" 

# --- PATHS (UPDATE THESE FOR YOUR CHOSEN SPECIES) ---
if(CURRENT_SPECIES == "Mouse") {
  FASTA_PATH <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa"
  GTF_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Mus_musculus.GRCm39.112.gtf"
  RDS_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/results/mouse/total/step01_data_prep/UFM1_events_rich.rds"
  
  ORG_DB    <- "org.Mm.eg.db"
  MART_HOST <- "https://www.ensembl.org"
  MART_NAME <- "ensembl"
  MART_DATA <- "mmusculus_gene_ensembl"
  CAI_SPECIES <- "mm"
  
} else if(CURRENT_SPECIES == "Arabidopsis") {
  FASTA_PATH <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/arabidopsis/Arabidopsis_thaliana_TAIR10.dna.primary_assembly.fa"
  GTF_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/arabidopsis/no_plastid_no_rRNA.Arabidopsis_thaliana.TAIR10.56.gtf"
  RDS_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/results/arabidopsis/nucleus/step01_data_prep/UFM1_events_rich.rds"
  
  ORG_DB    <- "org.At.tair.db"
  MART_HOST <- "https://plants.ensembl.org"
  MART_NAME <- "plants_mart"
  MART_DATA <- "athaliana_eg_gene"
  CAI_SPECIES <- "at" 
  
} else { # Default to Human
  FASTA_PATH <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
  GTF_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/pcg_gencode.v45.annotation.gtf"
  RDS_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/results/human/nucleus/step01_data_prep/UFM1_events_rich.rds"
  
  ORG_DB    <- "org.Hs.eg.db"
  MART_HOST <- "https://www.ensembl.org"
  MART_NAME <- "ensembl"
  MART_DATA <- "hsapiens_gene_ensembl"
  CAI_SPECIES <- "hs"
}

# --- LIBRARIES ---
suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(Biostrings)
  library(tidyverse)
  library(ggpubr)
  library(biomaRt)
  library(clusterProfiler)
  library(seqinr)
  
  if(!require(ORG_DB, character.only=TRUE)) {
    BiocManager::install(ORG_DB)
    library(ORG_DB, character.only=TRUE)
  }
})

# ==============================================================================
# 2. DATA LOADING & EXTRACTION (OPTIMIZED)
# ==============================================================================
message(paste("\n>>> RUNNING ANALYSIS FOR:", CURRENT_SPECIES, "<<<"))
message("Loading data...")

# Load all data upfront
txdb <- makeTxDbFromGFF(GTF_PATH, format="gtf")
events <- readRDS(RDS_PATH)
genome_fa <- FaFile(FASTA_PATH)

# OPTIMIZATION: Cache genes and CDS once
message("Caching genome annotations...")
genes_gr_cache <- genes(txdb)
cds_by_tx_cache <- cdsBy(txdb, by="tx", use.names=TRUE)
tx_gene_map_cache <- AnnotationDbi::select(txdb, keys=names(cds_by_tx_cache), 
                                           columns="GENEID", keytype="TXNAME")

# Define Groups
ri_events <- events[events$EventType == "RI", ]
dep_ri    <- ri_events[ri_events$Group == "UFM1_dependent", ]
indep_ri  <- ri_events[ri_events$Group == "UFM1_independent", ]

# Remove Overlaps
overlaps <- findOverlaps(indep_ri, dep_ri)
if (length(overlaps) > 0) indep_ri <- indep_ri[-unique(queryHits(overlaps))]

message(sprintf("Groups: %d dependent, %d independent", length(dep_ri), length(indep_ri)))

# OPTIMIZED: Clean IDs (use cached genes)
get_clean_genes <- function(gr, species, genes_cache) {
  seqlevelsStyle(gr) <- seqlevelsStyle(genes_cache)[1]
  hits <- findOverlaps(gr, genes_cache)
  ids <- unique(genes_cache$gene_id[subjectHits(hits)])
  if(species != "Arabidopsis") {
    return(sub("\\..*", "", ids))
  } else {
    return(ids)
  }
}

genes_dep   <- get_clean_genes(dep_ri, CURRENT_SPECIES, genes_gr_cache)
genes_indep <- get_clean_genes(indep_ri, CURRENT_SPECIES, genes_gr_cache)
universe_genes <- unique(get_clean_genes(genes_gr_cache, CURRENT_SPECIES, genes_gr_cache))
background_genes <- setdiff(universe_genes, c(genes_dep, genes_indep))

# OPTIMIZED: Get C-Terminus (use cached CDS)
get_c_term <- function(gr, genes_cache, cds_cache, tx_map) {
  seqlevelsStyle(gr) <- seqlevelsStyle(genes_cache)[1]
  hits <- findOverlaps(gr, genes_cache)
  if(length(hits)==0) return(NULL)
  
  targets <- unique(genes_cache$gene_id[subjectHits(hits)])
  map <- tx_map[tx_map$GENEID %in% targets & !is.na(tx_map$GENEID),]
  target_cds <- cds_cache[map$TXNAME]
  
  # Select best transcript per gene (longest CDS)
  sel <- data.frame(
    tx = names(target_cds), 
    gene = map$GENEID[match(names(target_cds), map$TXNAME)], 
    len = sapply(target_cds, function(x) sum(width(x)))
  )
  best <- sel %>% group_by(gene) %>% slice_max(len, n=1, with_ties=FALSE)
  
  seqs <- extractTranscriptSeqs(genome_fa, target_cds[best$tx])
  valid <- width(seqs) >= 90
  return(list(dna=subseq(seqs[valid], start=width(seqs[valid])-89, end=width(seqs[valid]))))
}

message("Extracting C-terminal sequences...")
data_dep_c <- get_c_term(dep_ri, genes_gr_cache, cds_by_tx_cache, tx_gene_map_cache)
data_indep_c <- get_c_term(indep_ri, genes_gr_cache, cds_by_tx_cache, tx_gene_map_cache)

# ==============================================================================
# 3. NEGATIVE CONTROLS (OPTIMIZED)
# ==============================================================================
message("\n>>> STEP 3: Intrinsic Features <<<")

# OPTIMIZED: Vectorized peptide charge calculation
calc_pep <- function(dna, grp) {
  if(is.null(dna)) return(NULL)
  aa <- as.character(Biostrings::translate(dna$dna, if.fuzzy.codon="solve"))
  # Vectorized string counting
  basic_counts <- str_count(aa, "[KR]")
  aa_lengths <- nchar(aa)
  data.frame(Group=grp, Basic_Density=basic_counts/aa_lengths)
}

df_pep <- rbind(calc_pep(data_dep_c, "Dependent"), calc_pep(data_indep_c, "Independent"))
p_pep <- ggboxplot(df_pep, x="Group", y="Basic_Density", palette="npg", 
                   add="jitter", title="A. Peptide Charge") + stat_compare_means()

# B. Codon Usage
data(caitab)
w_species <- caitab[[CAI_SPECIES]]

calc_cai <- function(dna, grp) {
  if(is.null(dna)) return(NULL)
  vals <- sapply(as.character(dna$dna), function(x) {
    s <- s2c(x)
    if(length(s)%%3!=0 || any(s=="N")) return(NA)
    tryCatch(cai(s, w=w_species), error=function(e) NA)
  })
  data.frame(Group=grp, CAI=vals)
}

df_cai <- rbind(calc_cai(data_dep_c, "Dependent"), calc_cai(data_indep_c, "Independent"))
p_cai <- ggboxplot(df_cai, x="Group", y="CAI", palette="npg", 
                   add="jitter", title="B. Codon Usage") + stat_compare_means()

# ==============================================================================
# 4. LOCALIZATION (OPTIMIZED - BATCH QUERIES)
# ==============================================================================
message("\n>>> STEP 4: Localization (BioMart) <<<")

p_loc <- tryCatch({
  message("Connecting to BioMart...")
  mart <- useMart(MART_NAME, dataset = MART_DATA, host = MART_HOST)
  
  filter_type <- if(CURRENT_SPECIES == "Arabidopsis") "ensembl_gene_id" else "ensembl_gene_id"
  
  # OPTIMIZATION: Single BioMart query instead of multiple
  all_query_genes <- c(genes_dep, genes_indep, background_genes)
  message(sprintf("Querying %d genes from BioMart...", length(all_query_genes)))
  
  go_data <- getBM(
    attributes=c('ensembl_gene_id','go_id'), 
    filters=filter_type, 
    values=all_query_genes, 
    mart=mart
  )
  
  # GO Terms
  loc_terms <- list(
    "Nucleus"="GO:0005634", 
    "Cytoplasm"="GO:0005737", 
    "Mitochondria"="GO:0005739",
    "ER_Network"="GO:0005783", 
    "Secreted"=c("GO:0005576", "GO:0005615"), 
    "Membrane"="GO:0016020"
  )
  if(CURRENT_SPECIES == "Arabidopsis") loc_terms[["Chloroplast"]] <- "GO:0009507"
  
  # OPTIMIZATION: Vectorized enrichment calculation
  calc_enrich <- function(target, bg, lbl, term) {
    hits_t <- sum(target %in% go_data$ensembl_gene_id[go_data$go_id %in% term])
    hits_b <- sum(bg %in% go_data$ensembl_gene_id[go_data$go_id %in% term])
    
    mat <- matrix(c(hits_t, length(target)-hits_t, hits_b, length(bg)-hits_b), nrow=2)
    ft <- fisher.test(mat)
    data.frame(Group=lbl, OR=ft$estimate, Lower=ft$conf.int[1], 
               Upper=ft$conf.int[2], P=ft$p.value)
  }
  
  # Calculate all enrichments
  results <- lapply(names(loc_terms), function(loc) {
    term <- loc_terms[[loc]]
    rbind(
      cbind(Category=loc, calc_enrich(genes_dep, background_genes, "Dependent", term)),
      cbind(Category=loc, calc_enrich(genes_indep, background_genes, "Independent", term))
    )
  })
  
  df_loc <- do.call(rbind, results)
  df_loc$Sig <- ifelse(df_loc$P < 0.05, "*", "")
  
  ggplot(df_loc, aes(x=Category, y=OR, fill=Group)) +
    geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper), position=position_dodge(width=0.8), width=0.25) +
    geom_text(aes(y=Upper*1.1, label=Sig), position=position_dodge(width=0.8), size=5) +
    scale_y_log10() + coord_flip() + theme_classic() + 
    scale_fill_manual(values=c("#E57373", "#64B5F6")) +
    labs(title=paste(CURRENT_SPECIES, "Localization"), y="Odds Ratio", x="")
  
}, error = function(e) {
  message(sprintf("BioMart failed: %s", e$message))
  ggplot() + theme_void() + labs(title="BioMart Error")
})

# ==============================================================================
# 5. FUNCTIONAL ENRICHMENT (OPTIMIZED)
# ==============================================================================
message("\n>>> STEP 5: Functional Enrichment <<<")

key_type <- if(CURRENT_SPECIES == "Arabidopsis") "TAIR" else "ENSEMBL"

p_enrich <- tryCatch({
  message("Converting gene IDs...")
  # OPTIMIZATION: Batch ID conversion
  all_conv <- bitr(c(genes_dep, genes_indep), fromType=key_type, 
                   toType="ENTREZID", OrgDb=get(ORG_DB))
  
  entrez_dep <- all_conv$ENTREZID[all_conv[[1]] %in% genes_dep]
  entrez_ind <- all_conv$ENTREZID[all_conv[[1]] %in% genes_indep]
  
  message(sprintf("Running enrichment: %d dep, %d indep", 
                  length(entrez_dep), length(entrez_ind)))
  
  ck <- compareCluster(
    geneCluster = list(Dependent=entrez_dep, Independent=entrez_ind),
    fun = "enrichGO", 
    OrgDb = get(ORG_DB), 
    ont = "BP", 
    pvalueCutoff = 0.05
  )
  
  # Save results
  ego_df <- as.data.frame(ck) %>% filter(Cluster == "Dependent")
  output_file <- paste0("UFM1_", CURRENT_SPECIES, "_Pathway_Drivers.csv")
  write.csv(ego_df, output_file, row.names=FALSE)
  message(sprintf("Saved pathway drivers to: %s", output_file))
  
  dotplot(ck, showCategory=10, title=paste(CURRENT_SPECIES, "Pathways")) + 
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
}, error = function(e) {
  message(sprintf("Enrichment failed: %s", e$message))
  ggplot() + theme_void() + labs(title="Enrichment Error")
})

# ==============================================================================
# 6. FINAL SUMMARY
# ==============================================================================
message("\n>>> Generating final figure <<<")

final_fig <- ggarrange(
  ggarrange(p_pep, p_cai, p_loc, ncol=3), 
  p_enrich,
  nrow=2, 
  heights=c(1, 1.5)
)

print(final_fig)
message(paste("\n✓ Analysis for", CURRENT_SPECIES, "complete!"))
