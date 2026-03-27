
library(dplyr)
library(tidyr)
library(readr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.At.tair.db)
library(AnnotationDbi)

# --- CONFIG ---
HUMAN_RESULTS <- "GO_Final_test/human/nucleus/step14_go_enrichment/Comparative_GO_Results.tsv"
MOUSE_RESULTS <- "GO_Final_test/mouse/total/step14_go_enrichment/Comparative_GO_Results.tsv"
ARAB_RESULTS  <- "GO_Final_test/arabidopsis/nucleus/step14_go_enrichment/Comparative_GO_Results.tsv"
OUTPUT_DIR <- "GO_Final_test/cross_species_visualization"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- LOAD DATA ---
cat("Loading results...\n")
h_df <- readr::read_tsv(HUMAN_RESULTS, show_col_types = FALSE)
m_df <- readr::read_tsv(MOUSE_RESULTS, show_col_types = FALSE)
a_df <- readr::read_tsv(ARAB_RESULTS, show_col_types = FALSE)

# --- HELPER: MAP METADATA ---
get_metadata <- function(ids, species) {
  if(species == "Human") {
    orgdb <- org.Hs.eg.db; kt <- "ENTREZID"
  } else if(species == "Mouse") {
    orgdb <- org.Mm.eg.db; kt <- "ENTREZID"
  } else {
    orgdb <- org.At.tair.db; kt <- "ENTREZID"
  }
  
  target_cols <- intersect(columns(orgdb), c("SYMBOL", "GENENAME"))
  if(length(target_cols) == 0) return(data.frame(GeneID=ids, SYMBOL=ids, GENENAME=ids))
  
  meta <- tryCatch({
    AnnotationDbi::select(orgdb, keys=as.character(ids), 
                          columns=target_cols, 
                          keytype=kt)
  }, error = function(e) return(data.frame(GeneID=ids, SYMBOL=ids, GENENAME=ids)))
  
  # Normalize column names
  colnames(meta)[colnames(meta) == kt] <- "GeneID"
  return(meta)
}

# --- UNIFIED GENE EXTRACTION ---
process_cluster <- function(cluster_name, arab_keywords, p_mam=0.01, p_at=0.05) {
  cat("\n>>> Processing Eukaryotic", cluster_name, "gene lists...\n")
  
  # 1. Mammalian Shared IDs
  h_ids <- h_df$ID[h_df$Cluster == cluster_name & h_df$pvalue < p_mam]
  m_ids <- m_df$ID[m_df$Cluster == cluster_name & m_df$pvalue < p_mam]
  mam_shared <- intersect(h_ids, m_ids)
  
  # 2. Arabidopsis Thematic Terms
  at_themes <- a_df %>% 
    filter(Cluster == cluster_name & pvalue < p_at) %>%
    filter(grepl(arab_keywords, Description, ignore.case=TRUE)) %>%
    filter(!grepl("rRNA", Description, ignore.case=TRUE)) %>%
    pull(ID)
  
  # Filter Mammalian Shared as well
  mam_shared_filtered <- h_df %>% 
    filter(ID %in% mam_shared) %>% 
    filter(!grepl("rRNA", Description, ignore.case=TRUE)) %>% 
    pull(ID) %>% unique()

  all_ids <- unique(c(mam_shared_filtered, at_themes))
  
  if(length(all_ids) == 0) {
    cat("No shared terms found for", cluster_name, "\n")
    return(NULL)
  }
  
  # Create directory
  cluster_dir <- file.path(OUTPUT_DIR, "Genes_By_Conserved_Term", paste0("Eukaryotic_", cluster_name))
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
  
  cluster_expanded <- data.frame()
  
  for (go_id in all_ids) {
    # Column IDs
    h_col <- colnames(h_df)[grep("geneid", colnames(h_df), ignore.case=TRUE)[1]]
    m_col <- colnames(m_df)[grep("geneid", colnames(m_df), ignore.case=TRUE)[1]]
    a_col <- colnames(a_df)[grep("geneid", colnames(a_df), ignore.case=TRUE)[1]]
    
    # Human
    h_sub <- h_df[h_df$ID == go_id & h_df$Cluster == cluster_name, c("ID", "Description", h_col)]
    if(nrow(h_sub)>0) {
      colnames(h_sub)[3] <- "GeneID"
      h_exp <- h_sub %>% mutate(GeneID = strsplit(as.character(GeneID), "/")) %>% unnest(GeneID) %>% mutate(Species = "Human")
      h_exp <- h_exp %>% left_join(get_metadata(unique(h_exp$GeneID), "Human"), by = "GeneID")
    } else { h_exp <- data.frame() }
    
    # Mouse
    m_sub <- m_df[m_df$ID == go_id & m_df$Cluster == cluster_name, c("ID", "Description", m_col)]
    if(nrow(m_sub)>0) {
      colnames(m_sub)[3] <- "GeneID"
      m_exp <- m_sub %>% mutate(GeneID = strsplit(as.character(GeneID), "/")) %>% unnest(GeneID) %>% mutate(Species = "Mouse")
      m_exp <- m_exp %>% left_join(get_metadata(unique(m_exp$GeneID), "Mouse"), by = "GeneID")
    } else { m_exp <- data.frame() }
    
    # Arabidopsis
    a_sub <- a_df[a_df$ID == go_id & a_df$Cluster == cluster_name, c("ID", "Description", a_col)]
    if(nrow(a_sub)>0) {
      colnames(a_sub)[3] <- "GeneID"
      a_exp <- a_sub %>% mutate(GeneID = strsplit(as.character(GeneID), "/")) %>% unnest(GeneID) %>% mutate(Species = "Arabidopsis")
      a_exp <- a_exp %>% left_join(get_metadata(unique(a_exp$GeneID), "Arabidopsis"), by = "GeneID")
    } else { a_exp <- data.frame() }
    
    term_meta <- rbind(h_exp, m_exp, a_exp)
    if(nrow(term_meta) > 0) {
      term_meta <- term_meta %>% dplyr::select(ID, Description, Species, GeneID, SYMBOL, GENENAME)
      safe_name <- gsub(":", "_", go_id)
      write.table(term_meta, file.path(cluster_dir, paste0(safe_name, ".tsv")), 
                  sep="\t", quote=FALSE, row.names=FALSE)
      cluster_expanded <- rbind(cluster_expanded, term_meta)
    }
  }
  
  # Save Master
  out_path <- file.path(OUTPUT_DIR, paste0("Eukaryotic_Conserved_", cluster_name, "_GO_Genes.tsv"))
  write.table(cluster_expanded, out_path, sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saved results to:", out_path, "\n")
}

# --- EXECUTION ---
process_cluster("Dependent", "lipid|endoplasmic|ER organization")
process_cluster("Independent", "splicing|RNA metabolic|RNA processing")
