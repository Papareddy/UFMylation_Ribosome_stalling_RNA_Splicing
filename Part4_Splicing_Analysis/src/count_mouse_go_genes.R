
library(dplyr)
library(tidyr)

# --- CONFIG ---
HUMAN_DATA <- "mammalian_RI_dpsi01_fdr05/human/nucleus/step14_go_enrichment/Differential_GO_Barplot_Data.tsv"
MOUSE_RESULTS <- "mammalian_RI_dpsi01_fdr05/mouse/total/step14_go_enrichment/Comparative_GO_Results.tsv"
OUTPUT_FILE <- "mammalian_RI_dpsi01_fdr05/mouse/total/step14_go_enrichment/Mouse_Counts_For_Human_Sig_GO.tsv"

# --- LOAD DATA ---
cat("Loading Human data from:", HUMAN_DATA, "\n")
human_df <- read.table(HUMAN_DATA, header=TRUE, sep="\t", quote="", comment.char="")

# Column 3 is the ID in the Human file (Description, Cluster, ID, ...)
# But let's check columns names of human_df
# head output suggested: Description Cluster ID GeneRatio ...
# Wait, let's look at the head output I got again:
# Description Cluster ID GeneRatio BgRatio ...
# So human_df$ID should be the GO IDs.

cat("Loading Mouse data from:", MOUSE_RESULTS, "\n")
mouse_df <- read.table(MOUSE_RESULTS, header=TRUE, sep="\t", quote="", comment.char="")

# Mouse columns: Cluster ID Description GeneRatio ...
# So mouse_df$ID should be the GO IDs.

# --- PROCESS ---

# 1. Get significant GO IDs from Human (Only from Dependent cluster)
human_sig_ids <- human_df %>%
  filter(Cluster == "Dependent") %>%
  pull(ID)

cat("Found", length(human_sig_ids), "significant Human GO terms.\n")

# 2. Extract Mouse counts for these IDs (Only from Dependent cluster)
mouse_counts <- mouse_df %>%
  filter(Cluster == "Dependent" & ID %in% human_sig_ids) %>%
  select(ID, Description, Count)

# 3. Join with Human descriptions to ensure we have the names
# Use human_df to get descriptions for the IDs we are looking for
human_names <- human_df %>%
  filter(Cluster == "Dependent") %>%
  select(ID, Description) %>%
  distinct()

result_df <- human_names %>%
  left_join(mouse_counts, by = "ID", suffix = c(".human", ".mouse")) %>%
  mutate(Mouse_Count = ifelse(is.na(Count), 0, Count)) %>%
  select(ID, Description.human, Mouse_Count)

# Rename columns for clarity
colnames(result_df) <- c("GO_ID", "Description", "Mouse_UFM1_Dependent_Gene_Count")

# --- OUTPUT ---
cat("Saving results to:", OUTPUT_FILE, "\n")
write.table(result_df, OUTPUT_FILE, sep="\t", quote=FALSE, row.names=FALSE)

cat("\nTop 10 results:\n")
print(head(result_df, 10))
