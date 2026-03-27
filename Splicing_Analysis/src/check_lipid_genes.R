library(tidyverse)
library(biomaRt)

# --- PATHS ---
FASTA_PATH <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/data/mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF_PATH   <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Mus_musculus.GRCm39.112.gtf"
RDS_PATH   <- "/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/results/mouse/total/step01_data_prep/UFM1_events_rich.rds"

library(GenomicFeatures)
library(GenomicRanges)

# Load Mouse events
txdb <- makeTxDbFromGFF(GTF_PATH, format="gtf")
events <- readRDS(RDS_PATH)

# Helper: Clean IDs
get_clean_genes <- function(gr) {
  seqlevelsStyle(gr) <- seqlevelsStyle(txdb)[1]
  genes_gr <- genes(txdb); hits <- findOverlaps(gr, genes_gr)
  ids <- unique(genes_gr$gene_id[subjectHits(hits)])
  return(sub("\\..*", "", ids))
}

# Get Dependent Genes (RI)
ri_events <- events[events$EventType == "RI", ]
dep_ri <- ri_events[ri_events$Group == "UFM1_dependent", ]
genes_dep <- get_clean_genes(dep_ri)

# Get Independent Genes (RI)
indep_ri <- ri_events[ri_events$Group == "UFM1_independent", ]
genes_indep <- get_clean_genes(indep_ri)

# Query BioMart for all RI genes
all_ri_genes <- unique(c(genes_dep, genes_indep))
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://www.ensembl.org")
go_data <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'go_id', 'name_1006'), 
                 filters = 'ensembl_gene_id', 
                 values = all_ri_genes, 
                 mart = mart)

# Filter for phospholipid/glycerophospholipid terms
lipid_genes <- go_data %>% 
  filter(str_detect(name_1006, fixed("phospholipid", ignore_case=TRUE)) | 
         str_detect(name_1006, fixed("glycerophospholipid", ignore_case=TRUE))) %>%
  dplyr::select(ensembl_gene_id, external_gene_name, go_id, name_1006) %>%
  distinct()

# Assign Groups
lipid_genes$Group <- case_when(
  lipid_genes$ensembl_gene_id %in% genes_dep & lipid_genes$ensembl_gene_id %in% genes_indep ~ "Both",
  lipid_genes$ensembl_gene_id %in% genes_dep ~ "UFM1_dependent",
  lipid_genes$ensembl_gene_id %in% genes_indep ~ "UFM1_independent",
  TRUE ~ "Unknown"
)

# Summarize lipid genes
result_summary <- lipid_genes %>%
  group_by(external_gene_name, Group) %>%
  summarize(Terms = paste(unique(name_1006), collapse = "; "), .groups = "drop") %>%
  arrange(Group, external_gene_name)

# Query BioMart for ER involvement of these specific lipid genes
lipid_gene_ids <- result_summary$external_gene_name
go_cc_data <- getBM(attributes = c('external_gene_name', 'name_1006'), 
                    filters = 'external_gene_name', 
                    values = lipid_gene_ids, 
                    mart = mart)

# Filter for ER terms
er_involvement <- go_cc_data %>%
  filter(str_detect(name_1006, fixed("endoplasmic reticulum", ignore_case=TRUE))) %>%
  group_by(external_gene_name) %>%
  summarize(ER_Terms = paste(unique(name_1006), collapse = "; "), .groups = "drop")

# Merge with previous summary
final_summary <- result_summary %>%
  left_join(er_involvement, by = "external_gene_name") %>%
  mutate(Involved_in_ER = ifelse(!is.na(ER_Terms), "Yes", "No"))

print("--- Phospholipid Genes and ER Involvement ---")
print(as.data.frame(final_summary %>% dplyr::select(external_gene_name, Group, Involved_in_ER, ER_Terms)))
