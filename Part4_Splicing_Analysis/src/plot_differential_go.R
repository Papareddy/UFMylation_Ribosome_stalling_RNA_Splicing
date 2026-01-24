
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# --- CONFIG ---
RESULTS_FILE <- "mammalian_RI_dpsi01_fdr05/human/nucleus/step14_go_enrichment/Comparative_GO_Results.tsv"
OUTPUT_DIR <- "mammalian_RI_dpsi01_fdr05/human/nucleus/step14_go_enrichment"
OUTPUT_PLOT <- file.path(OUTPUT_DIR, "Differential_GO_Barplot")

# --- LOAD DATA ---
cat("Loading results from:", RESULTS_FILE, "\n")
df <- read.table(RESULTS_FILE, header=TRUE, sep="\t", quote="", comment.char="")

# Check Columns
if(!"Cluster" %in% colnames(df)) stop("Column 'Cluster' not found!")
if(!"p.adjust" %in% colnames(df)) stop("Column 'p.adjust' not found!")

# --- FILTER TOP 20 DEPENDENT ---
# 1. Select Dependent
dep <- df %>% filter(Cluster == "Dependent")

# 2. Sort by significance (p.adjust asc)
dep <- dep %>% arrange(p.adjust)

# 3. Take Top 20
top20_terms <- head(dep$Description, 20)
cat("Top 20 UFM1-Dependent Terms:\n")
print(top20_terms)

# --- PREPARE PLOTTING DATA ---
# We want these 20 terms. We want their values for BOTH Dependent and Independent.
plot_df <- df %>% filter(Description %in% top20_terms)

# Ensure both groups exist for every term
# Create a skeleton
skeleton <- expand.grid(Description = top20_terms, Cluster = c("Dependent", "Independent"))
plot_df <- merge(skeleton, plot_df, by=c("Description", "Cluster"), all.x=TRUE)

# Handle Missing Values (Terms not enriched in Independent)
# If missing, it means p.adjust is effectively 1 (not significant)
plot_df$p.adjust[is.na(plot_df$p.adjust)] <- 1
plot_df$logP <- -log10(plot_df$p.adjust)

# Truncate very long names for plot
plot_df$Description_Short <- str_trunc(plot_df$Description, 50)

# Factor ordering: Maintain sort order of Dependent
# We want the terms with highest -logP in Dependent to be at top of plot
# Get order from top20 (which was sorted)
order_levels <- rev(top20_terms) # rev for ggplot coord_flip
plot_df$Description <- factor(plot_df$Description, levels = order_levels)

# Colors
cols <- c("Dependent" = "#E41A1C", "Independent" = "#377EB8")

# --- PLOT ---
p <- ggplot(plot_df, aes(x=Description, y=logP, fill=Cluster)) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() +
  scale_fill_manual(values=cols) +
  labs(title="UFM1-Specific GO Enrichment (Human)",
       subtitle="Top 20 UFM1-Dependent Terms vs Independent Background",
       y="-log10(Adjusted P-value)",
       x="") +
  theme_minimal() +
  theme(axis.text.y = element_text(size=10),
        legend.position="bottom") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="gray50") +
  annotate("text", x=1, y=-log10(0.05), label="FDR 0.05", vjust=-1, color="gray50", size=3)

# Save
ggsave(paste0(OUTPUT_PLOT, ".pdf"), p, width=8, height=10)
ggsave(paste0(OUTPUT_PLOT, ".png"), p, width=8, height=10)

# Save Data
write.table(plot_df, file.path(OUTPUT_DIR, "Differential_GO_Barplot_Data.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cat("Saved plots and data to", OUTPUT_PLOT, "\n")
