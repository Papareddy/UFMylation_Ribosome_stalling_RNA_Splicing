# Libraries
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(ggrepel)
    library(optparse)
})

# ==========================================
# 1. ARGUMENT PARSING
# ==========================================

option_list <- list(
  make_option(c("--dep"), type="character", default=NULL, 
              help="Path to UFM1-dependent AME TSV file", metavar="IS_DEP"),
  make_option(c("--indep"), type="character", default=NULL, 
              help="Path to UFM1-independent AME TSV file", metavar="IS_INDEP"),
  make_option(c("--out"), type="character", default="motif_enrichment_plot.pdf", 
              help="Output PDF filename [default= %default]", metavar="OUTPUT"),
  make_option(c("--max_val"), type="double", default=25, 
              help="Maximum axis value (logE) [default= %default]", metavar="MAX_VAL"),
  make_option(c("--min_val"), type="double", default=-1, 
              help="Minimum axis value (logE) [default= %default]", metavar="MIN_VAL"),
  make_option(c("--top_n"), type="integer", default=10, 
              help="Number of top motifs to label per side [default= %default]", metavar="TOP_N")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$dep) || is.null(opt$indep)) {
  print_help(opt_parser)
  stop("Both --dep and --indep input files must be provided.", call.=FALSE)
}

# ==========================================
# 2. FUNCTIONS
# ==========================================

load_ame <- function(filepath) {
  if (!file.exists(filepath)) {
    warning(paste("File not found:", filepath))
    return(data.frame())
  }
  
  df <- read_tsv(filepath, comment = "#", show_col_types = FALSE)
  
  if (!"E-value" %in% colnames(df)) return(data.frame())
  
  df %>%
    mutate(
      logE = -log10(`E-value`),
      Display = paste0(ifelse(!is.na(motif_alt_ID) & motif_alt_ID != ".", motif_alt_ID, motif_ID), "_", consensus)
    ) %>%
    select(motif_ID, Display, `E-value`, logE)
}

# ==========================================
# 3. PROCESSING
# ==========================================

message("Loading Dependent Data: ", opt$dep)
df_dep <- load_ame(opt$dep) %>% 
  rename(logE_Dep = logE, E_Dep = `E-value`, Display_Dep = Display)

message("Loading Independent Data: ", opt$indep)
df_indep <- load_ame(opt$indep) %>% 
  rename(logE_Indep = logE, E_Indep = `E-value`, Display_Indep = Display)

if (nrow(df_dep) == 0 && nrow(df_indep) == 0) {
  stop("No data loaded from input files.")
}

merged <- full_join(df_dep, df_indep, by = "motif_ID") %>%
  mutate(
    logE_Dep = ifelse(is.na(logE_Dep), 0, logE_Dep),
    logE_Indep = ifelse(is.na(logE_Indep), 0, logE_Indep),
    Display = coalesce(Display_Dep, Display_Indep),
    diff = logE_Dep - logE_Indep
  )

# ==========================================
# 4. LABELING LOGIC
# ==========================================

SIG_THRESHOLD_LOGE <- 2       # E-value < 0.01

# Filter candidates (must be significant in at least one)
candidates <- merged %>%
  filter(logE_Dep > SIG_THRESHOLD_LOGE | logE_Indep > SIG_THRESHOLD_LOGE)

# 1. Top N Dep-Specific (Highest Positive Diff)
top_dep <- candidates %>% 
  filter(diff > 0) %>% 
  slice_max(order_by = diff, n = opt$top_n)

# 2. Top N Indep-Specific (Lowest Negative Diff)
top_indep <- candidates %>% 
  filter(diff < 0) %>% 
  slice_min(order_by = diff, n = opt$top_n)

# 3. SRSF5 Specific Logic (User Requested: Always Add SRSF5)
srsf5_cand <- merged %>% 
  filter(grepl("SRSF5", Display, ignore.case = TRUE))

# Combine Labels
labels_to_show <- bind_rows(top_dep, top_indep, srsf5_cand) %>% 
  distinct(motif_ID, .keep_all = TRUE)

message(paste0("Labels generated: ", nrow(labels_to_show)))

# ==========================================
# 5. PLOTTING (GGPLOT2)
# ==========================================

# Create a Label Column in the main dataframe
merged <- merged %>%
  mutate(Label = ifelse(motif_ID %in% labels_to_show$motif_ID, Display, NA))

p <- ggplot(merged, aes(x = logE_Indep, y = logE_Dep)) +
  # Diagonal Line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  
  # Points with gradient color based on difference
  geom_point(aes(fill = diff), pch = 21, color = "black", size = 3, alpha = 0.8) +
  
  # Color Scale (Blue = Indep, Red = Dep)
  scale_fill_gradient2(
    low = "#1f77b4", 
    mid = "white", 
    high = "#d73027", 
    midpoint = 0, 
    name = "Enrichment\nBias"
  ) +
  
  # Labels with Repel
  geom_text_repel(
    aes(label = Label),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.5,
    min.segment.length = 0
  ) +
  
  # Highlight SRSF5 specifically if present
  geom_point(data = srsf5_cand, aes(x = logE_Indep, y = logE_Dep), 
             color = "purple", size = 4, stroke = 1.5, shape = 1) +
  
  # Scales
  scale_x_continuous(limits = c(opt$min_val, opt$max_val)) +
  scale_y_continuous(limits = c(opt$min_val, opt$max_val)) +
  coord_fixed(ratio = 1) +
  
  # Theme
  theme_classic() +
  labs(
    title = "Motif Enrichment Comparison",
    subtitle = "UFM1 Dependent vs Independent",
    x = "-log10(E-value) [Independent]",
    y = "-log10(E-value) [Dependent]"
  )

ggsave(opt$out, p, width = 8, height = 8)
message("Plot saved to: ", opt$out)
