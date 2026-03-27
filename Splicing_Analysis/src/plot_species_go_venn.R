
library(VennDiagram)
library(dplyr)

# --- CONFIG ---
HUMAN_RESULTS <- "GO_Final_test/human/nucleus/step14_go_enrichment/Comparative_GO_Results.tsv"
MOUSE_RESULTS <- "GO_Final_test/mouse/total/step14_go_enrichment/Comparative_GO_Results.tsv"
OUTPUT_DIR <- "GO_Final_test/cross_species_visualization"
OUT_PDF <- file.path(OUTPUT_DIR, "GO_Overlap_Venn.pdf")
OUT_TSV <- file.path(OUTPUT_DIR, "GO_Overlap_Intersection_Terms.tsv")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- LOAD DATA ---
cat("Loading GO results...\n")
human_df <- read.table(HUMAN_RESULTS, header=TRUE, sep="\t", quote="", comment.char="")
mouse_df <- read.table(MOUSE_RESULTS, header=TRUE, sep="\t", quote="", comment.char="")

# --- FILTER SIGNIFICANT DEPENDENT TERMS ---
# threshold: pvalue < 0.01
human_sig <- human_df %>% 
  filter(Cluster == "Dependent" & pvalue < 0.01) %>% 
  pull(ID)

mouse_sig <- mouse_df %>% 
  filter(Cluster == "Dependent" & pvalue < 0.01) %>% 
  pull(ID)

cat("Human Significant Dependent GO Terms:", length(human_sig), "\n")
cat("Mouse Significant Dependent GO Terms:", length(mouse_sig), "\n")

# --- CALCULATE INTERSECTION ---
intersection_ids <- intersect(human_sig, mouse_sig)
cat("Intersection Count:", length(intersection_ids), "\n")

# Save the intersecting terms for reference
intersection_terms <- human_df %>%
  filter(ID %in% intersection_ids) %>%
  select(ID, Description) %>%
  distinct()

write.table(intersection_terms, OUT_TSV, sep="\t", quote=FALSE, row.names=FALSE)

# --- PLOT VENN ---
cat("Generating Venn Diagram...\n")

# VennDiagram writes to a log file by default, we disable it
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn.plot <- venn.diagram(
  x = list(Human = human_sig, Mouse = mouse_sig),
  filename = NULL,
  fill = c("#E41A1C", "#377EB8"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = c(-20, 20),
  cat.dist = 0.05,
  main = "Overlap of Significant UFM1-Dependent GO Terms",
  sub = "FDR < 0.05",
  main.cex = 1.5,
  sub.cex = 1
)

pdf(OUT_PDF, width=8, height=8)
grid.draw(venn.plot)
dev.off()

cat("Venn diagram saved to:", OUT_PDF, "\n")
cat("Intersection terms saved to:", OUT_TSV, "\n")
