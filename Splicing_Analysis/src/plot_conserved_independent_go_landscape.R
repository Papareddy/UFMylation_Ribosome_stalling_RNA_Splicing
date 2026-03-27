
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# --- CONFIG ---
HUMAN_RESULTS <- "GO_Final_test/human/nucleus/step14_go_enrichment/Comparative_GO_Results.tsv"
MOUSE_RESULTS <- "GO_Final_test/mouse/total/step14_go_enrichment/Comparative_GO_Results.tsv"
ARAB_RESULTS  <- "GO_Final_test/arabidopsis/nucleus/step14_go_enrichment/Comparative_GO_Results.tsv"
OUTPUT_DIR <- "GO_Final_test/cross_species_visualization"
OUT_PREFIX <- file.path(OUTPUT_DIR, "Eukaryotic_Conserved_Independent_Landscape")

# --- LOAD DATA ---
cat("Loading results...\n")
h_full <- readr::read_tsv(HUMAN_RESULTS, show_col_types = FALSE)
m_full <- readr::read_tsv(MOUSE_RESULTS, show_col_types = FALSE)
a_full <- readr::read_tsv(ARAB_RESULTS, show_col_types = FALSE)

# 1. Mammalian Shared IDs (Intersection p < 0.01)
h_ind_sig <- h_full$ID[h_full$Cluster == "Independent" & h_full$pvalue < 0.01]
m_ind_sig <- m_full$ID[m_full$Cluster == "Independent" & m_full$pvalue < 0.01]
mam_shared_ids <- intersect(h_ind_sig, m_ind_sig)

# 2. Arabidopsis Thematic Terms (Splicing/RNA metabolic p < 0.05)
at_themes_ids <- a_full %>% 
  filter(Cluster == "Independent" & pvalue < 0.05) %>%
  filter(grepl("splicing|RNA metabolic|RNA processing", Description, ignore.case=TRUE)) %>%
  filter(!grepl("rRNA", Description, ignore.case=TRUE)) %>%
  pull(ID)

combined_landscape_ids <- unique(c(mam_shared_ids, at_themes_ids))

cat("Found", length(mam_shared_ids), "mammalian shared independent terms and", length(at_themes_ids), "Arabidopsis thematic terms.\n")

# --- PREPARE PLOT DATA ---
h_df <- h_full %>% filter(Cluster == "Independent" & ID %in% combined_landscape_ids) %>% select(ID, Description, pvalue, Count) %>% mutate(Species = "Human")
m_df <- m_full %>% filter(Cluster == "Independent" & ID %in% combined_landscape_ids) %>% select(ID, Description, pvalue, Count) %>% mutate(Species = "Mouse")
a_df <- a_full %>% filter(Cluster == "Independent" & ID %in% combined_landscape_ids) %>% select(ID, Description, pvalue, Count) %>% mutate(Species = "Arabidopsis")

plot_df <- rbind(h_df, m_df, a_df)
plot_df$logP <- -log10(plot_df$pvalue)
plot_df$Species <- factor(plot_df$Species, levels=c("Human", "Mouse", "Arabidopsis"))

# Create Label with GO ID
plot_df$Label <- paste0(plot_df$Description, " (", plot_df$ID, ")")

# Sort by average logP across species
avg_sig <- plot_df %>% 
  group_by(Label) %>% 
  summarise(avg = mean(logP, na.rm=TRUE)) %>% 
  arrange(avg)

plot_df$Label <- factor(plot_df$Label, levels = avg_sig$Label)

# --- PLOT ---
cat("Generating Bubble Plot...\n")
p <- ggplot(plot_df, aes(x = Species, y = Label)) +
  geom_point(aes(size = Count, color = logP)) +
  scale_color_viridis_c(option = "viridis", name = "-log10(p-value)", na.value="transparent", end = 0.9) +
  scale_size_continuous(name = "Gene Count") +
  theme_bw() +
  labs(title = "Eukaryotic Conserved UFM1-Independent Functional Landscape",
       subtitle = "Conserved themes (RNA Splicing) across Human, Mouse, and Arabidopsis",
       x = "", y = "") +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 12, face="bold"),
        panel.grid.major = element_line(color = "gray95", linetype = "dashed"))

# Save
ggsave(paste0(OUT_PREFIX, ".pdf"), p, width=10, height=8)
ggsave(paste0(OUT_PREFIX, ".png"), p, width=10, height=8)

# Save Data
write.table(plot_df, paste0(OUT_PREFIX, "_Data.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cat("Saved results and data to:", OUT_PREFIX, "\n")
