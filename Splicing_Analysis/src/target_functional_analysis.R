#!/usr/bin/env Rscript
# ==============================================================================
# UNIVERSAL UFM1 ANALYSIS: HUMAN / MOUSE / ARABIDOPSIS
# Original version - See target_functional_analysis_optimized.R for improvements
# ==============================================================================

# --- 1. SPECIES CONFIGURATION (CHANGE THIS!) ---
# Options: "Human", "Mouse", "Arabidopsis"
CURRENT_SPECIES <- "Arabidopsis" 

# --- PATHS (UPDATE THESE FOR YOUR CHOSEN SPECIES) ---
if(CURRENT_SPECIES == "Mouse") {
  FASTA_PATH <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa"
  GTF_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Mus_musculus.GRCm39.112.gtf"
  RDS_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/results/mouse/total/step01_data_prep/UFM1_events_rich.rds"
  
  
  # Database Settings
  ORG_DB    <- "org.Mm.eg.db"
  MART_HOST <- "https://www.ensembl.org"
  MART_NAME <- "ensembl"
  MART_DATA <- "mmusculus_gene_ensembl"
  CAI_SPECIES <- "mm" # for seqinr caitab
  
} else if(CURRENT_SPECIES == "Arabidopsis") {
  FASTA_PATH <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/arabidopsis/Arabidopsis_thaliana_TAIR10.dna.primary_assembly.fa"
  GTF_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/arabidopsis/no_plastid_no_rRNA.Arabidopsis_thaliana.TAIR10.56.gtf"
  RDS_PATH   <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/results/arabidopsis/nucleus/step01_data_prep/UFM1_events_rich.rds"
  
  
  # Database Settings
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
  
  # Load Organism DB dynamically
  if(!require(ORG_DB, character.only=TRUE)) {
    BiocManager::install(ORG_DB)
    library(ORG_DB, character.only=TRUE)
  }
})

# ==============================================================================
# 2. DATA LOADING & EXTRACTION
# ==============================================================================
message(paste("\n>>> RUNNING ANALYSIS FOR:", CURRENT_SPECIES, "<<<"))

txdb <- makeTxDbFromGFF(GTF_PATH, format="gtf")
events <- readRDS(RDS_PATH)
genome_fa <- FaFile(FASTA_PATH)

# Define Groups
ri_events <- events[events$EventType == "RI", ]
dep_ri    <- ri_events[ri_events$Group == "UFM1_dependent", ]
indep_ri  <- ri_events[ri_events$Group == "UFM1_independent", ]

# Remove Overlaps
overlaps <- findOverlaps(indep_ri, dep_ri)
if (length(overlaps) > 0) indep_ri <- indep_ri[-unique(queryHits(overlaps))]

# Helper: Clean IDs (Handle Ensembl versions vs TAIR IDs)
get_clean_genes <- function(gr, species) {
  seqlevelsStyle(gr) <- seqlevelsStyle(txdb)[1]
  genes_gr <- genes(txdb); hits <- findOverlaps(gr, genes_gr)
  ids <- unique(genes_gr$gene_id[subjectHits(hits)])
  if(species != "Arabidopsis") {
    return(sub("\\..*", "", ids)) # Strip version for animals
  } else {
    return(ids) # Keep TAIR IDs (AT1G01010) as is
  }
}

genes_dep   <- get_clean_genes(dep_ri, CURRENT_SPECIES)
genes_indep <- get_clean_genes(indep_ri, CURRENT_SPECIES)
all_genes_gr <- genes(txdb)
universe_genes <- unique(get_clean_genes(all_genes_gr, CURRENT_SPECIES))
background_genes <- setdiff(universe_genes, c(genes_dep, genes_indep))

# Helper: Get C-Terminus (90bp)
get_c_term <- function(gr) {
  seqlevelsStyle(gr) <- seqlevelsStyle(txdb)[1]
  genes_gr <- genes(txdb); hits <- findOverlaps(gr, genes_gr)
  if(length(hits)==0) return(NULL)
  targets <- unique(genes_gr$gene_id[subjectHits(hits)])
  cds_by_tx <- cdsBy(txdb, by="tx", use.names=TRUE)
  map <- AnnotationDbi::select(txdb, keys=names(cds_by_tx), columns="GENEID", keytype="TXNAME")
  map <- map[map$GENEID %in% targets & !is.na(map$GENEID),]
  target_cds <- cds_by_tx[map$TXNAME]
  sel <- data.frame(tx=names(target_cds), gene=map$GENEID[match(names(target_cds), map$TXNAME)], len=sum(width(target_cds)))
  best <- sel %>% group_by(gene) %>% slice_max(len, n=1, with_ties=FALSE)
  seqs <- extractTranscriptSeqs(genome_fa, target_cds[best$tx])
  valid <- width(seqs) >= 90
  return(list(dna=subseq(seqs[valid], start=width(seqs[valid])-89, end=width(seqs[valid]))))
}

data_dep_c <- get_c_term(dep_ri)
data_indep_c <- get_c_term(indep_ri)

# ==============================================================================
# 3. NEGATIVE CONTROLS
# ==============================================================================
message("\n>>> STEP 3: Intrinsic Features <<<")

# A. Peptide Charge
calc_pep <- function(dna, grp) {
  if(is.null(dna)) return(NULL)
  aa <- as.character(Biostrings::translate(dna$dna, if.fuzzy.codon="solve"))
  data.frame(Group=grp, Basic_Density=str_count(aa, "[KR]")/nchar(aa))
}
df_pep <- rbind(calc_pep(data_dep_c, "Dependent"), calc_pep(data_indep_c, "Independent"))
p_pep <- ggboxplot(df_pep, x="Group", y="Basic_Density", palette="npg", add="jitter", title="A. Peptide Charge") + stat_compare_means()

# B. Codon Usage (Dynamic Weights)
data(caitab)
w_species <- caitab[[CAI_SPECIES]] # Loads 'hs', 'mm', or 'at' weights

calc_cai <- function(dna, grp) {
  if(is.null(dna)) return(NULL)
  vals <- sapply(as.character(dna$dna), function(x) {
    s <- s2c(x); if(length(s)%%3!=0 || any(s=="N")) return(NA); tryCatch(cai(s, w=w_species), error=function(e) NA)
  })
  data.frame(Group=grp, CAI=vals)
}
df_cai <- rbind(calc_cai(data_dep_c, "Dependent"), calc_cai(data_indep_c, "Independent"))
p_cai <- ggboxplot(df_cai, x="Group", y="CAI", palette="npg", add="jitter", title="B. Codon Usage") + stat_compare_means()

# ==============================================================================
# 4. LOCALIZATION (SPECIES SPECIFIC BIOMART)
# ==============================================================================
message("\n>>> STEP 4: Localization (BioMart) <<<")

tryCatch({
  mart <- useMart(MART_NAME, dataset = MART_DATA, host = MART_HOST)
  
  # Note: Arabidopsis often uses TAIR IDs for filter
  filter_type <- if(CURRENT_SPECIES == "Arabidopsis") "ensembl_gene_id" else "ensembl_gene_id"
  
  go_data <- getBM(attributes=c('ensembl_gene_id','go_id'), filters=filter_type, values=c(genes_dep, genes_indep, background_genes), mart=mart)
  
  # Standard GO Terms apply to all eukaryotes
  loc_terms <- list("Nucleus"="GO:0005634", "Cytoplasm"="GO:0005737", "Mitochondria"="GO:0005739",
                    "ER_Network"="GO:0005783", "Secreted"=c("GO:0005576", "GO:0005615"), "Membrane"="GO:0016020")
  
  # Add Chloroplast for Plants
  if(CURRENT_SPECIES == "Arabidopsis") loc_terms[["Chloroplast"]] <- "GO:0009507"
  
  results_list <- list()
  for(loc in names(loc_terms)) {
    term <- loc_terms[[loc]]
    
    calc_enrich <- function(target, bg, lbl) {
      hits_t <- sum(target %in% go_data$ensembl_gene_id[go_data$go_id %in% term])
      hits_b <- sum(bg %in% go_data$ensembl_gene_id[go_data$go_id %in% term])
      ft <- fisher.test(matrix(c(hits_t, length(target)-hits_t, hits_b, length(bg)-hits_b), nrow=2))
      data.frame(Category=loc, Group=lbl, OR=ft$estimate, Lower=ft$conf.int[1], Upper=ft$conf.int[2], P=ft$p.value)
    }
    results_list[[paste0(loc,"_dep")]] <- calc_enrich(genes_dep, background_genes, "Dependent")
    results_list[[paste0(loc,"_ind")]] <- calc_enrich(genes_indep, background_genes, "Independent")
  }
  
  df_loc <- do.call(rbind, results_list)
  df_loc$Sig <- ifelse(df_loc$P < 0.05, "*", "")
  
  p_loc <- ggplot(df_loc, aes(x=Category, y=OR, fill=Group)) +
    geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper), position=position_dodge(width=0.8), width=0.25) +
    geom_text(aes(y=Upper*1.1, label=Sig), position=position_dodge(width=0.8), size=5) +
    scale_y_log10() + coord_flip() + theme_classic() + scale_fill_manual(values=c("#E57373", "#64B5F6")) +
    labs(title=paste(CURRENT_SPECIES, "Localization"), y="Odds Ratio", x="")
  
}, error = function(e) {
  message("BioMart failed (Connection/Timeout). Skipping Localization Plot.")
  p_loc <- ggplot() + theme_void() + labs(title="BioMart Error")
})

# ==============================================================================
# 5. FUNCTIONAL ENRICHMENT (SPECIES SPECIFIC)
# ==============================================================================
message("\n>>> STEP 5: Functional Enrichment <<<")

# ID Mapping (Ensembl -> Entrez is safest for clusterProfiler)
# Arabidopsis mapping keytype is typically 'TAIR' or 'ENSEMBL' depending on DB version
key_type <- if(CURRENT_SPECIES == "Arabidopsis") "TAIR" else "ENSEMBL"

tryCatch({
  entrez_dep <- bitr(genes_dep, fromType=key_type, toType="ENTREZID", OrgDb=get(ORG_DB))$ENTREZID
  entrez_ind <- bitr(genes_indep, fromType=key_type, toType="ENTREZID", OrgDb=get(ORG_DB))$ENTREZID
  
  ck <- compareCluster(geneCluster = list(Dependent=entrez_dep, Independent=entrez_ind),
                       fun = "enrichGO", OrgDb = get(ORG_DB), ont = "BP", pvalueCutoff = 0.05)
  
  p_enrich <- dotplot(ck, showCategory=10, title=paste(CURRENT_SPECIES, "Pathways")) + 
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  # Save Drivers
  ego_df <- as.data.frame(ck) %>% filter(Cluster == "Dependent")
  write.csv(ego_df, paste0("UFM1_", CURRENT_SPECIES, "_Pathway_Drivers.csv"))
  
}, error = function(e) {
  message("Enrichment failed (ID Mapping?). Skipping Dotplot.")
  p_enrich <- ggplot() + theme_void() + labs(title="Enrichment Error")
})

# ==============================================================================
# 6. FINAL SUMMARY
# ==============================================================================
final_fig <- ggarrange(ggarrange(p_pep, p_cai, p_loc, ncol=3), p_enrich, nrow=2, heights=c(1, 1.5))
print(final_fig)
message(paste("\nDone! Analysis for", CURRENT_SPECIES, "complete."))
