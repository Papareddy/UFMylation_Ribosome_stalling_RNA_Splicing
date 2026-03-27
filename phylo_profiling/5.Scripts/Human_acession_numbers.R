library(dplyr)
library(readr)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)

par(mfrow=c(4,4),las=2, tcl=-0.3, bty="n")

setwd("~/Downloads/UFM1_profiling/UFM1_correlation/")

###############################################################
## 1. Load clusters and mark significance
###############################################################
clusters <- read_tsv("UFM1.correlated_clusters.n100000.tab", col_types = cols(.default = "c")) %>%
  mutate(p_value = as.numeric(pvalue), pbc = as.numeric(pbc))

sig_clusters <- clusters %>% filter(p_value <= 0.05, pbc >= 0.25)
message("Significant clusters: ", nrow(sig_clusters))

###############################################################
## 2. Load orthogroup assignments and extract HUMAN entries
###############################################################
assignments <- read_tsv("UFM1.orthogroup_assignment.tab", col_names = FALSE, col_types = cols(.default = "c"))
colnames(assignments) <- c("species_uniprot", "cluster")

human_hits <- assignments %>%
  mutate(
    species = str_extract(species_uniprot, "^[0-9]+"),
    uniprot = str_extract(species_uniprot, "(?<=\\.)[A-Za-z0-9]+")
  ) %>%
  filter(cluster %in% sig_clusters$cluster, species == "9606") %>%
  select(cluster, uniprot) %>%
  distinct()

message("Human UniProt entries extracted: ", nrow(human_hits))

###############################################################
## 3. Fast Local Mapping (Replaces UniProt API and BioMart)
###############################################################
message("Mapping annotations locally...")

uniprot_keys <- unique(human_hits$uniprot)

# Map UniProt to Gene Symbol and Ensembl ID
gene_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = uniprot_keys,
  columns = c("SYMBOL", "ENSEMBL", "GENENAME"),
  keytype = "UNIPROT"
))

# Map UniProt to GO IDs and Categories (Ontology)
go_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys = uniprot_keys,
  columns = c("GO", "ONTOLOGY"),
  keytype = "UNIPROT"
)) %>% filter(!is.na(GO))

# Get GO Term Names from GO.db
go_terms <- suppressMessages(AnnotationDbi::select(
  GO.db,
  keys = unique(go_map$GO),
  columns = "TERM",
  keytype = "GOID"
))

# Combine GO IDs with their names and collapse per UniProt ID
go_full <- go_map %>%
  left_join(go_terms, by = c("GO" = "GOID")) %>%
  group_by(UNIPROT) %>%
  summarise(
    GO_terms = paste(unique(GO), collapse = "; "),
    GO_names = paste(unique(TERM), collapse = "; "),
    GO_category = paste(unique(ONTOLOGY), collapse = "; ")
  )

###############################################################
## 4. Merge ALL annotations together + cluster metrics
###############################################################
final_annotated <- human_hits %>%
  left_join(gene_map, by = c("uniprot" = "UNIPROT")) %>%
  left_join(go_full, by = c("uniprot" = "UNIPROT")) %>%
  left_join(clusters, by = "cluster") %>%
  distinct()

###############################################################
## 5. Save output
###############################################################
final_annotated_reorder=as.data.frame(final_annotated[c(3,1:2,4:7,17,18,11:16)])%>% 
  arrange(SYMBOL)
write.table(
  final_annotated_reorder,
  "~/Downloads/UFM1_profiling/UFM1_correlation/UFMylation_CoevolvedGenes_human_orthologs_annotated_with_GO.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Saved: human_significant_clusters_annotated_with_GO.tsv")




library(dplyr)
library(ggplot2)
library(ggrepel)

par(mfrow=c(4,4),las=2, tcl=-0.3, bty="n")

###############################################################
## 1. Prepare Background Cluster Data
###############################################################
clusters <- clusters %>%
  mutate(
    pbc = as.numeric(pbc),
    p_value = as.numeric(p_value),
    neglog10p = -log10(p_value),
    ycap = pmin(neglog10p, 5)
  )

###############################################################
## 2. Categorize and Prepare Overlay Data (FORPLOTING)
###############################################################
df <- FORPLOTING %>%
  select(gene_name = SYMBOL, pbc, p_value) %>%
  na.omit() %>%
  mutate(
    pbc = as.numeric(pbc),
    p_value = as.numeric(p_value),
    neglog10p = -log10(p_value),
    ycap = pmin(neglog10p, 5),
    category = case_when(
      grepl("^FANC", gene_name, ignore.case = TRUE) | 
        gene_name %in% c("FANCI", "FANCL", "FANCD2", "INTS7", "INTS3", "NABP2", "NABP1", "NHEJ1", "TREX1", "TREX2", "CCAR2") ~ "DNA damage",
      
      grepl("^PRDM|^PCGF", gene_name, ignore.case = TRUE) | 
        gene_name == "SETD7" ~ "Chromatin modifiers",
      
      grepl("^RBM|^ERI", gene_name, ignore.case = TRUE) | 
        gene_name %in% c("RNPC3", "RBM41", "NUP153", "THOC6", "PNN", "ZNF830", "TTF2", "GPATCH8", "CSTF2", "CSTF2T", "ZCCHC13", "SMG9", "ZCCHC7") ~ "RNA processing",
      
      gene_name %in% c("UBA5", "UFC1", "UFL1", "UFM1", "UFSP1", "UFSP2", "ODR4", "DDRGK1", "CDK5RAP3") ~ "UFMylation",
      
      TRUE ~ "Other"
    )
  ) %>%
  filter(category != "Other")





df <- FORPLOTING %>%
  select(gene_name = SYMBOL, pbc, p_value) %>%
  na.omit() %>%
  mutate(
    pbc = as.numeric(pbc),
    p_value = as.numeric(p_value),
    neglog10p = -log10(p_value),
    ycap = pmin(neglog10p, 5),
    category = case_when(
   
      
      grepl("^RTN|^TRAP", gene_name, ignore.case = TRUE) | 
        gene_name == "SETD7" ~ "Chromatin modifiers",
      
      
      
      TRUE ~ "Other"
    )
  ) %>%
  filter(category != "Other")



###############################################################
## 3. Color Palette
###############################################################
category_colors <- c(
  "RNA processing"      = "#2a568c",
  "DNA damage"          = "#4db7b1",
  "Chromatin modifiers" = "#d69545",
  "UFMylation"          = "#d96d83",
  "RNA degradation"     = "#aaaaaa",
  "Other"               = "grey60"
)

###############################################################
## 4. GGPlot Generation
###############################################################
p <- ggplot() +
  geom_point(data = clusters, aes(x = pbc, y = ycap), 
             shape = 21, fill = "grey80", color = "white", size = 2.8) +
  geom_vline(xintercept = 0.3, color = "grey20", linewidth = 0.8, linetype = "dashed") +
  geom_hline(yintercept = 1.30103, color = "grey20", linewidth = 0.8, linetype = "dashed") +
  geom_point(data = df, aes(x = pbc, y = ycap, fill = category), 
             shape = 21, color = "black", size = 4) +
  geom_text_repel(data = df, aes(x = pbc, y = ycap, label = gene_name, color = category), 
                  size = 3.5, max.overlaps = Inf, show.legend = FALSE) +
  scale_fill_manual(values = category_colors) +
  scale_color_manual(values = category_colors) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 5)) +
  labs(
    x = "PBC",
    y = "-log10(p-value)",
    title = "PBC vs -log10(p-value) with Category Overlay",
    fill = "Gene Category"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid = element_line(color = "grey90")
  )

print(p)
