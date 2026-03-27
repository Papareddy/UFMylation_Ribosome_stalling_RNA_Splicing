# ========================================================
# TurboID: Optimized Annotation, Categorization, and Plot
# ========================================================

# Load required libraries
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db") # Run once if not installed
if (!require("beeswarm", quietly = TRUE)) install.packages("beeswarm")

library(dplyr)
library(stringr)
library(readr)
library(org.Hs.eg.db)
library(beeswarm)

# -------------------------
# 1) Load Data
# -------------------------
df <- read_tsv("/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/SuplimentalTables/Table-S2_hsUFM_Coevolved_proteins_AF2_MM.tsv")

# -------------------------
# 2) Local Annotation (UniProt -> Gene Symbol)
# -------------------------
ids <- unique(df$prey)

gene_map <- mapIds(org.Hs.eg.db, 
                   keys = ids, 
                   column = "SYMBOL", 
                   keytype = "UNIPROT", 
                   multiVals = "first")

bm <- data.frame(prey = names(gene_map), 
                 gene_name = as.character(gene_map), 
                 stringsAsFactors = FALSE)

df_annot <- df %>% left_join(bm, by = "prey")

# -------------------------
# 3) Categorize and Filter
# -------------------------
thr <- 0.75

genes_rna_degradation <- c("SMG9", "ZCCHC7", "MEX3A")
genes_dna_damage      <- c("FANCI","FANCL","FANCD2","INTS7","INTS3","NABP2","NABP1","NHEJ1","TREX1","TREX2","CCAR2")
genes_chromatin       <- c("SETD7")
genes_rna_processing  <- c("RNPC3","RBM41","NUP153","THOC6","PNN","ZNF830","TTF2","GPATCH8","CSTF2","CSTF2T","ZCCHC13")
genes_ufmylation      <- c("UBA5","UFC1","UFL1","UFM1","UFSP1","UFSP2","ODR4","DDRGK1","CDK5RAP3")

df_cat <- df_annot %>%
  mutate(
    category = case_when(
      str_detect(gene_name, "^ERI") | gene_name %in% genes_rna_degradation ~ "RNA degradation",
      str_detect(gene_name, "^FANC") | gene_name %in% genes_dna_damage ~ "DNA damage",
      str_detect(gene_name, "^(PRDM|PCGF)") | gene_name %in% genes_chromatin ~ "Chromatin modifiers",
      str_detect(gene_name, "^RBM") | gene_name %in% genes_rna_processing ~ "RNA processing",
      gene_name %in% genes_ufmylation ~ "UFMylation",
      TRUE ~ NA_character_
    ),
    plot_group = factor(ifelse(mean_top3_scaled_PEAK > thr & !is.na(category), category, "Neutral"),
                        levels = c("Neutral", "UFMylation", "DNA damage", "Chromatin modifiers", "RNA degradation", "RNA processing"))
  ) %>%
  arrange(plot_group)

# -------------------------
# 4) Plotting (Base R)
# -------------------------
par(mfrow=c(4,4),las=2, tcl=-0.3, bty="n")

my_colors <- c(
  "UFMylation"          = "#7991c9",
  "DNA damage"          = "#4cb6b0",
  "Chromatin modifiers" = "#d59543",
  "RNA degradation"     = "#d96d83",
  "RNA processing"      = "#d96d83",
  "Neutral"             = "#d1d3d3"
)

baits <- unique(df_cat$bait)

for (b in baits) {
  sub_df <- df_cat[df_cat$bait == b, ]
  
  bg_colors <- my_colors[as.character(sub_df$plot_group)]
  
  bs_coords <- beeswarm(sub_df$mean_top3_scaled_PEAK,
                        pwbg = bg_colors,
                        pch = 21,
                        col = "white",
                        ylim = c(0, 1.0),
                        main = b,
                        ylab = "Mean of top-k scaled_PEAK",
                        xlab = "",
                        labels = "All proteins",
                        method = "compactswarm")
  
  abline(h = thr, lty = 2, col = "grey50")
  
  label_idx <- which(sub_df$plot_group != "Neutral")
  
  if (length(label_idx) > 0) {
    text(x = bs_coords$x[label_idx], 
         y = bs_coords$y[label_idx], 
         labels = sub_df$gene_name[label_idx], 
         pos = 4, 
         cex = 0.8, 
         offset = 0.5)
  }
}




# -------------------------
# 4) Plotting
# -------------------------
par(mfrow=c(4,4),las=2, tcl=-0.3, bty="n")

# Set the order of the facets by converting 'bait' to a factor
df_cat$bait <- factor(df_cat$bait, levels = c("HsUFM1", "hsUFL1", "hsDDRGK1", "hsUBA5", "hsUFC1", "hsCDK5RAP3"))

# Sort data to render "Neutral" points on the bottom layer
df_cat <- df_cat %>% arrange(plot_group != "Neutral")

p <- ggplot(df_cat, aes(x = "All proteins", y = mean_top3_scaled_PEAK)) +
  geom_hline(yintercept = thr, linetype = "dashed", color = "grey50") +
  # Map size conditionally (3 * 1.25 = 3.75)
  geom_beeswarm(
    aes(
      fill = plot_group, 
      size = ifelse(plot_group == "Neutral", "Background", "Target")
    ), 
    shape = 21, 
    color = "white", 
    stroke = 0.25,
    method = "compactswarm", 
    cex = 2 
  ) +
  geom_text_repel(
    data = subset(df_cat, plot_group != "Neutral"),
    aes(label = gene_name),
    size = 3.5, 
    max.overlaps = 50, 
    box.padding = 0.5, 
    point.padding = 0.3,
    segment.color = "grey50",
    min.segment.length = 0
  ) +
  scale_fill_manual(values = my_colors) +
  # Apply the specific sizes and hide the size legend
  scale_size_manual(values = c("Background" = 5, "Target" = 5.75), guide = "none") +
  scale_y_continuous(limits = c(0, 1.0)) +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Mean of top-k scaled_PEAK") +
  facet_wrap(~bait) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.length.y = unit(-0.15, "cm"), 
    axis.text.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), color = "black"), 
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.5)
  )

print(p)
