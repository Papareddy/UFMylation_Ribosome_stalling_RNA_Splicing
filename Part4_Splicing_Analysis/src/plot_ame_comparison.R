
# ==============================================================================
# Comparison of Motif Enrichment: UFM1 Dependent vs Independent
# ==============================================================================

# Libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggrepel)
})

# Set Working Directory (as requested)
# Only run this if the directory exists (avoid errors on other systems)
wd_path <- "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis/"
if (dir.exists(wd_path)) {
    setwd(wd_path)
}

# ==========================================
# 1. SETUP & INPUTS
# ==========================================

# Default Paths (User Requested)
file_dep_default   <- "results/human/nucleus/step10_motif_analysis/RI_DeepDive_combined/ame_UFM1_dependent_vs_constitutive/ame.tsv"
file_indep_default <- "results/human/nucleus/step10_motif_analysis/RI_DeepDive_combined/ame_UFM1_independent_vs_constitutive/ame.tsv"
out_pdf_default    <- "results/human/nucleus/step10_motif_analysis/RI_DeepDive_combined/R_motif_comparison_scatter.pdf"

# Hybrid Argument Parsing: Check if CLI args are passed; if not, use defaults.
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(pattern) {
  match <- grep(pattern, args, value = TRUE)
  if (length(match) > 0) return(sub(pattern, "", match[1]))
  return(NULL)
}

if (length(args) >= 3) {
    message("Using CLI arguments...")
    file_dep   <- get_arg("--dep=")
    file_indep <- get_arg("--indep=")
    out_pdf    <- get_arg("--out=")
} else {
    message("No CLI arguments found. Using default hardcoded paths...")
    file_dep   <- file_dep_default
    file_indep <- file_indep_default
    out_pdf    <- out_pdf_default
}

message("Dependent TSV:   ", file_dep)
message("Independent TSV: ", file_indep)
message("Output PDF:      ", out_pdf)


# ==========================================
# 2. FUNCTIONS
# ==========================================

load_ame <- function(filepath) {
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(data.frame())
  }
  
  # Read TSV, skip comments handled by read_tsv automatically if comment="#"
  df <- read_tsv(filepath, comment = "#", show_col_types = FALSE)
  
  # Extract relevant columns
  if (!"E-value" %in% colnames(df)) return(data.frame())
  
  df %>%
    mutate(
      logE = -log10(`E-value`),
      # Create Display Name: Use Alt ID if available, else standard ID
      # AND append consensus sequence
      Display = paste0(ifelse(!is.na(motif_alt_ID) & motif_alt_ID != ".", motif_alt_ID, motif_ID), "_", consensus)
    ) %>%
    select(motif_ID, Display, `E-value`, logE)
}

# ==========================================
# 3. PROCESSING
# ==========================================

df_dep <- load_ame(file_dep) %>% 
  rename(logE_Dep = logE, E_Dep = `E-value`, Display_Dep = Display)

df_indep <- load_ame(file_indep) %>% 
  rename(logE_Indep = logE, E_Indep = `E-value`, Display_Indep = Display)

# Full join
merged <- full_join(df_dep, df_indep, by = "motif_ID")

# Fill NAs and Calculate Diff
merged <- merged %>%
  mutate(
    logE_Dep = replace_na(logE_Dep, 0),
    logE_Indep = replace_na(logE_Indep, 0),
    Display = coalesce(Display_Dep, Display_Indep),
    diff = logE_Dep - logE_Indep
  )

# ==========================================
# LABELING LOGIC (USER REQUESTED)
# ==========================================
# Criteria:
# 1. E-value < 0.01 (logE > 2)
# 2. Difference > 10x (abs(logE_Dep - logE_Indep) >= 1)

SIG_THRESHOLD_LOGE <- 2       # E-value < 0.01
DIFF_THRESHOLD     <- 1       # 10x difference in E-value

labels_to_show <- merged %>%
  filter(
    (logE_Dep > SIG_THRESHOLD_LOGE | logE_Indep > SIG_THRESHOLD_LOGE)
  )

# 1. Top 10 Dep-Specific (Highest Positive Diff)
top_dep <- labels_to_show %>% 
  filter(diff > 0) %>% 
  slice_max(order_by = diff, n = 10)

# 2. Top 10 Indep-Specific (Lowest Negative Diff)
top_indep <- labels_to_show %>% 
  filter(diff < 0) %>% 
  slice_min(order_by = diff, n = 10)

# 3. SRSF Logic: Top 5 per side (must meet >10x diff)
srsf_candidates <- labels_to_show %>%
  filter(
    grepl("SRSF", Display, ignore.case = TRUE),
    abs(diff) >= DIFF_THRESHOLD
  )

srsf_dep <- srsf_candidates %>%
  filter(diff > 0) %>%
  slice_max(order_by = diff, n = 5)

srsf_indep <- srsf_candidates %>%
  filter(diff < 0) %>%
  slice_min(order_by = diff, n = 5)

# Combine unique labels
labels_to_show <- bind_rows(top_dep, top_indep, srsf_dep, srsf_indep) %>% distinct(motif_ID, .keep_all = TRUE)

message(paste0("Number of labeled (Top 10 Gen + Top 5 SRSF per side): ", nrow(labels_to_show)))

# Start Plotting
max_val <- max(c(merged$logE_Dep, merged$logE_Indep)) * 1.1

# ==========================================
# 4. PLOTTING
# ==========================================

p <- ggplot(merged, aes(x = logE_Indep, y = logE_Dep, fill = diff)) +
  # Diagonal line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  
  # Scatter points
  geom_point(shape = 21, size = 3, color = "black", alpha = 0.8) +
  
  # Color Gradient (Blue -> White -> Red)
  scale_fill_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027", 
    midpoint = 0, name = "Enrichment Bias\n(Dep - Indep)"
  ) +
  
  # Labels
  geom_text_repel(
    data = labels_to_show,
    aes(label = Display),
    size = 3,
    box.padding = 0.5,
    max.overlaps = Inf,   
    segment.color = "grey50",
    min.segment.length = 0
  ) +
  
  # Scales & Themes
  scale_x_continuous(limits = c(0, max_val), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max_val), expand = c(0, 0)) +
  coord_fixed() + 
  theme_classic() +
  labs(
    title = "Motif Enrichment Comparison",
    subtitle = "Labels: Top 10 Bias + Top 5 SRSF (>10x diff)",
    x = "UFM1 Independent Enrichment (-log10 E-value)",
    y = "UFM1 Dependent Enrichment (-log10 E-value)"
  )

# print(p)

# Save
ggsave(out_pdf, plot = p, width = 10, height = 10)
message("Plot saved to: ", out_pdf)
