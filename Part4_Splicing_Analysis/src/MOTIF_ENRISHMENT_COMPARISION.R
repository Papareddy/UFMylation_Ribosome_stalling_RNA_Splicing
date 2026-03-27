# ============================================================================
# STEP 1: SETUP
# ============================================================================
suppressPackageStartupMessages(library(tidyverse))

# Set working directory
setwd("/Users/ranjithpapareddy/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing/Part4_Splicing_Analysis")

# Define file paths
dep_file    <- "mammalian_RI_dpsi01_fdr05/human/nucleus/step08_motif_analysis/RI_DeepDive_combined/ame_intron_UFM1_dependent_vs_constitutive/ame.tsv"
indep_file  <- "mammalian_RI_dpsi01_fdr05/human/nucleus/step08_motif_analysis/RI_DeepDive_combined/ame_intron_UFM1_independent_vs_constitutive/ame.tsv"
output_pdf  <- "motif_enrichment_curves_final.pdf"

if (!file.exists(dep_file)) stop("Dependent file not found!")
if (!file.exists(indep_file)) stop("Independent file not found!")

# ============================================================================
# STEP 2: LOAD & PROCESS DATA
# ============================================================================
df_dep   <- read_tsv(dep_file, comment = "#", show_col_types = FALSE)
df_indep <- read_tsv(indep_file, comment = "#", show_col_types = FALSE)

# Helper to extract consensus
get_consensus <- function(df) {
  possible_names <- c("consensus", "cons", "motif_sequence", "sequence")
  match <- intersect(names(df), possible_names)
  if (length(match) > 0) return(df[[match[1]]])
  return(rep("", nrow(df)))
}

prep_df <- function(df, suffix) {
  df$temp_cons <- get_consensus(df)
  df %>%
    mutate(
      Display = ifelse(!is.na(motif_alt_ID) & motif_alt_ID != ".", motif_alt_ID, motif_ID),
      consensus = temp_cons,
      logE = -log10(pmax(`E-value`, 1e-300))
    ) %>%
    select(motif_ID, Display, consensus, `E-value`, logE) %>%
    rename_with(~paste0(., "_", suffix), -motif_ID)
}

d_dep   <- prep_df(df_dep, "Dep")
d_indep <- prep_df(df_indep, "Indep")

merged <- d_dep %>%
  full_join(d_indep, by = "motif_ID") %>%
  mutate(
    logE_Dep   = replace_na(logE_Dep, 0),
    logE_Indep = replace_na(logE_Indep, 0),
    Display    = coalesce(Display_Dep, Display_Indep),
    consensus  = coalesce(consensus_Dep, consensus_Indep),
    diff       = logE_Dep - logE_Indep,
    abs_diff   = abs(diff)
  )

# Create Labels
merged$label_text <- apply(merged, 1, function(row) {
  disp <- row["Display"]
  seq  <- row["consensus"]
  if(is.na(seq) || seq == "" || seq == ".") seq <- "" 
  if(seq != "") paste0(disp, "\n(", seq, ")") else disp
})

# ============================================================================
# STEP 3: PLOT SETTINGS (VARIABLES YOU REQUESTED)
# ============================================================================
# --- Styling ---
alpha_val   <- 0.9       # Transparency (0 to 1)
outline_col <- "white"   # Outline color (e.g., "white" or "black")

# --- Thresholds ---
min_enrichment <- 2    # Cutoff for Grey vs Color

# --- Curve Parameters ---
curve_c <- 1.0           # Minimum distance from diagonal (c)
curve_k <- 1.5           # Curvature (k)

# --- Size Scaling ---
min_pt_size <- 1       # Smallest point size
max_pt_size <- 3       # Largest point size

# ============================================================================
# STEP 4: ASSIGN LOGIC
# ============================================================================
# 1. Colors (Red/Blue/Grey) with Alpha
# Helper to add alpha to any color string
make_transparent <- function(col, alpha) {
  adjustcolor(col, alpha.f = alpha)
}

merged$pt_bg <- make_transparent("#929393", alpha_val) # Default
merged$pt_bg[merged$diff >= min_enrichment]  <- make_transparent("#d96d83", alpha_val)
merged$pt_bg[merged$diff <= -min_enrichment] <- make_transparent("#7991c9", alpha_val)

# 2. Sizes (Based on Abs Diff)
max_diff <- max(merged$abs_diff, na.rm=TRUE)
merged$pt_cex <- min_pt_size + (merged$abs_diff / max_diff) * (max_pt_size - min_pt_size)

# 3. Label Selection (Top 10 colored only)
sig_data <- merged[merged$abs_diff >= min_enrichment, ]
top_pos  <- sig_data %>% filter(diff > 0) %>% slice_max(diff, n = 25)
top_neg  <- sig_data %>% filter(diff < 0) %>% slice_min(diff, n = 25)
labels_to_plot <- bind_rows(top_pos, top_neg)

# ============================================================================
# STEP 5: PLOTTING
# ============================================================================


#par(mfrow=c(2,2), mar=c(5,5,4,2)+0.1, las=1, bty="n", pty="s")

# 1. Canvas (-0.5 to 30)
plot(1, type="n", 
     xlim=c(-0.5, 30), ylim=c(-0.5, 30),
     xlab="UFM1 Independent (-log10 E-value)", 
     ylab="UFM1 Dependent (-log10 E-value)",
     main="Motif Enrichment Comparison",
     axes=FALSE)

# 2. Grid & Diagonal
grid_locs <- seq(0, 30, 5)
abline(h=grid_locs, v=grid_locs, col="grey92", lwd=1)
abline(0, 1, col="grey50", lty=2, lwd=1.5)

# 3. DRAW CURVES (Hyperbolic)
# Formula: y = x + c + k/x  (and inverse)
s <- seq(0.1, 30, length.out=500)
# Upper Curve (Dep > Indep)
lines(s, s + curve_c + (curve_k/s), lty=2, col="black", lwd=1)
# Lower Curve (Indep > Dep)
lines(s + curve_c + (curve_k/s), s, lty=2, col="black", lwd=1)

# 4. Draw Points
# Sort so larger/colored points are plotted last (on top)
plot_order <- order(merged$abs_diff, decreasing = FALSE)
merged_ord <- merged[plot_order, ]

points(merged_ord$logE_Indep, merged_ord$logE_Dep, 
       pch = 21, 
       bg  = merged_ord$pt_bg, 
       col = outline_col,      # Variable outline color
       lwd = 0.5, 
       cex = merged_ord$pt_cex)

# 5. Labels
text(x = labels_to_plot$logE_Indep, 
     y = labels_to_plot$logE_Dep, 
     labels = labels_to_plot$label_text, 
     pos = 4, cex = 0.7, col = "black", offset = 0.5)

# 6. Axes
axis(1, at=seq(0, 30, 5), lwd=0, lwd.ticks=1)
axis(2, at=seq(0, 30, 5), lwd=0, lwd.ticks=1)
box(bty="o", col="black")
# 7. Legends
# A. Color Legend (Bottom Right)
legend("bottomright", 
       legend=c("Dependent", "Independent", "Non-Sig"), 
       pt.bg=c(make_transparent("red", alpha_val), 
               make_transparent("blue", alpha_val), 
               make_transparent("grey50", alpha_val)), 
       pch=21, col=outline_col, pt.cex=1.5, bty="n", 
       inset=c(0.02, 0.02))

# B. Curve Legend (Top Left)
legend("topleft", 
       legend = paste0("Significance Curve\n(c=", curve_c, ", k=", curve_k, ")"), 
       lty = 2, col = "black", bty = "n", cex = 0.8)

# C. Size Legend (Top Left, below Curve Legend)
size_vals <- c(5, 10, 20)
size_cex  <- min_pt_size + (size_vals / max_diff) * (max_pt_size - min_pt_size)

legend("topleft", 
       legend = paste0("Diff = ", size_vals), 
       pt.cex = size_cex, 
       pch = 21, pt.bg = "grey50", col = outline_col, 
       bty = "n", title = "Enrichment Size",
       inset = c(0, 0.12)) # 0.12 shifts it down below the curve legend