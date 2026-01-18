
# ==============================================================================
# Comparison of Motif Enrichment: UFM1 Dependent vs Independent
# ==============================================================================

# Libraries
library(tidyverse)
library(ggrepel)

# ==========================================
# 1. SETUP & INPUTS
# ==========================================

# Argument Parsing
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_ame_comparison.R --dep=<path> --indep=<path> --out=<path>")
}

get_arg <- function(pattern) {
  match <- grep(pattern, args, value = TRUE)
  if (length(match) > 0) return(sub(pattern, "", match[1]))
  return(NULL)
}

file_dep   <- get_arg("--dep=")
file_indep <- get_arg("--indep=")
out_pdf    <- get_arg("--out=")

if (is.null(file_dep) || is.null(file_indep) || is.null(out_pdf)) {
  stop("Missing required arguments. Need --dep=, --indep=, --out=")
}

message("Dependent TSV: ", file_dep)
message("Independent TSV: ", file_indep)
message("Output PDF: ", out_pdf)


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
  # Standardize column selection
  if (!"E-value" %in% colnames(df)) return(data.frame())
  
  df %>%
    mutate(
      logE = -log10(`E-value`),
      # Create Display Name: Use Alt ID if available, else standard ID
      Display = ifelse(!is.na(motif_alt_ID) & motif_alt_ID != ".", motif_alt_ID, motif_ID)
    ) %>%
    select(motif_ID, Display, `E-value`, logE)
}

# ==========================================
# 3. PROCESSING
# ==========================================

message("Loading Dependent: ", file_dep)
df_dep <- load_ame(file_dep) %>% 
  rename(logE_Dep = logE, E_Dep = `E-value`, Display_Dep = Display)

message("Loading Independent: ", file_indep)
df_indep <- load_ame(file_indep) %>% 
  rename(logE_Indep = logE, E_Indep = `E-value`, Display_Indep = Display)

# Merge
# Full join to keep motifs finding in either set
merged <- full_join(df_dep, df_indep, by = "motif_ID")

# Fill NAs
# If missing in one, it means E-value > threshold. Set logE to 0.
merged <- merged %>%
  mutate(
    logE_Dep = replace_na(logE_Dep, 0),
    logE_Indep = replace_na(logE_Indep, 0),
    # Combine Display Names
    Display = coalesce(Display_Dep, Display_Indep),
    # Calculate Difference (Positive = Dep Specific, Negative = Indep Specific)
    diff = logE_Dep - logE_Indep
  )

# Define Significance Threshold for Labeling
SIG_THRESHOLD_LOGE <- 2  # E-value < 0.01

# Identify top divergent motifs for labeling
top_dep_specific <- merged %>%
  filter(diff > 0, logE_Dep > SIG_THRESHOLD_LOGE) %>%
  slice_max(order_by = diff, n = 10)

top_indep_specific <- merged %>%
  filter(diff < 0, logE_Indep > SIG_THRESHOLD_LOGE) %>%
  slice_min(order_by = diff, n = 10)

labels_to_show <- bind_rows(top_dep_specific, top_indep_specific)

# Max limit for axes
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
    max.overlaps = 20,
    segment.color = "grey50"
  ) +
  
  # Scales & Themes
  scale_x_continuous(limits = c(0, max_val), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max_val), expand = c(0, 0)) +
  coord_fixed() + # Square plotting area
  theme_classic() +
  labs(
    title = "Motif Enrichment Comparison",
    subtitle = "UFM1 Dependent (Lost) vs Independent (Preserved) Introns",
    x = "UFM1 Independent Enrichment (-log10 E-value)",
    y = "UFM1 Dependent Enrichment (-log10 E-value)"
  )

print(p)

# Save
ggsave(out_pdf, plot = p, width = 8, height = 8)
message("Plot saved to: ", out_pdf)
