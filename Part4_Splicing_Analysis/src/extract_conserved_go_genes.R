
library(dplyr)
library(tidyr)

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)

library(readr)

# --- CONFIG ---
HUMAN_RESULTS <- "GO_Final_test/human/nucleus/step14_go_enrichment/Comparative_GO_Results.tsv"
MOUSE_RESULTS <- "GO_Final_test/mouse/total/step14_go_enrichment/Comparative_GO_Results.tsv"
OUTPUT_DIR <- "GO_Final_test/cross_species_visualization"
OUT_DEP <- file.path(OUTPUT_DIR, "Conserved_Dependent_GO_Genes.tsv")
OUT_INDEP <- file.path(OUTPUT_DIR, "Conserved_Independent_GO_Genes.tsv")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- LOAD DATA ---
cat("Loading results...\n")
h_df <- readr::read_tsv(HUMAN_RESULTS, show_col_types = FALSE)
m_df <- readr::read_tsv(MOUSE_RESULTS, show_col_types = FALSE)

cat("Human columns:", paste(colnames(h_df), collapse=", "), "\n")
cat("Mouse columns:", paste(colnames(m_df), collapse=", "), "\n")

# --- HELPER: MAP METADATA ---
get_metadata <- function(ids, species) {
  orgdb <- if(species == "Human") org.Hs.eg.db else org.Mm.eg.db
  meta <- AnnotationDbi::select(orgdb, keys=as.character(ids), 
                                columns=c("SYMBOL", "GENENAME"), 
                                keytype="ENTREZID")
  return(meta)
}

# --- HELPER FUNCTION: EXTRACT AND FORMAT ---
extract_conserved_genes <- function(cluster_name, p_thresh = 0.01) {
  # 1. Identify shared GO IDs
  h_sig <- h_df$ID[h_df$Cluster == cluster_name & h_df$pvalue < p_thresh]
  m_sig <- m_df$ID[m_df$Cluster == cluster_name & m_df$pvalue < p_thresh]
  shared_ids <- intersect(h_sig, m_sig)
  
  if(length(shared_ids) == 0) return(data.frame())
  
  # Create cluster-specific directory
  cluster_dir <- file.path(OUTPUT_DIR, "Genes_By_Conserved_Term", cluster_name)
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Robust column identification
  h_gi_col <- colnames(h_df)[grep("geneid", colnames(h_df), ignore.case=TRUE)[1]]
  m_gi_col <- colnames(m_df)[grep("geneid", colnames(m_df), ignore.case=TRUE)[1]]
  
  all_expanded <- data.frame()
  
  for (go_id in shared_ids) {
    # Extract Human
    h_sub <- h_df[h_df$ID == go_id & h_df$Cluster == cluster_name, c("ID", "Description", h_gi_col)]
    colnames(h_sub)[3] <- "EntrezID"
    h_expanded <- h_sub %>%
      mutate(EntrezID = strsplit(as.character(EntrezID), "/")) %>%
      unnest(EntrezID) %>%
      mutate(Species = "Human")
    
    # Extract Mouse
    m_sub <- m_df[m_df$ID == go_id & m_df$Cluster == cluster_name, c("ID", "Description", m_gi_col)]
    colnames(m_sub)[3] <- "EntrezID"
    m_expanded <- m_sub %>%
      mutate(EntrezID = strsplit(as.character(EntrezID), "/")) %>%
      unnest(EntrezID) %>%
      mutate(Species = "Mouse")
    
    # Add Metadata
    if(nrow(h_expanded) > 0) {
      h_meta <- get_metadata(unique(h_expanded$EntrezID), "Human")
      h_expanded <- h_expanded %>% left_join(h_meta, by = c("EntrezID" = "ENTREZID"))
    }
    if(nrow(m_expanded) > 0) {
      m_meta <- get_metadata(unique(m_expanded$EntrezID), "Mouse")
      m_expanded <- m_expanded %>% left_join(m_meta, by = c("EntrezID" = "ENTREZID"))
    }
    
    combined_meta <- rbind(h_expanded, m_expanded) %>%
      dplyr::select(ID, Description, Species, EntrezID, SYMBOL, GENENAME)
    
    # Save individual table
    safe_name <- gsub(":", "_", go_id)
    write.table(combined_meta, file.path(cluster_dir, paste0(safe_name, ".tsv")), 
                sep="\t", quote=FALSE, row.names=FALSE)
    
    all_expanded <- rbind(all_expanded, combined_meta)
  }
    
  return(all_expanded)
}

# --- PROCESS DEPENDENT ---
cat("Processing Dependent cluster...\n")
dep_genes <- extract_conserved_genes("Dependent")
if(nrow(dep_genes) > 0) {
  write.table(dep_genes, OUT_DEP, sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saved Dependent genes to:", OUT_DEP, "\n")
}

# --- PROCESS INDEPENDENT ---
cat("Processing Independent cluster...\n")
indep_genes <- extract_conserved_genes("Independent")
if(nrow(indep_genes) > 0) {
  write.table(indep_genes, OUT_INDEP, sep="\t", quote=FALSE, row.names=FALSE)
  cat("Saved Independent genes to:", OUT_INDEP, "\n")
}
