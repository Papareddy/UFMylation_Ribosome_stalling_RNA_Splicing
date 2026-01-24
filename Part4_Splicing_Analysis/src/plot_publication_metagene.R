#!/usr/bin/env Rscript

# ==============================================================================
# SCIENCE PUBLICATION QUALITY PROTEIN METAGENE PLOTTER (BASE R)
# ==============================================================================

# Output configuration
OUTPUT_PDF <- "protein_metagene_science_figure.pdf"
A4_WIDTH <- 8.27
A4_HEIGHT <- 11.69

# Species and Paths
SPECIES_ORDER <- c("Arabidopsis", "Human", "Mouse")
BASE_DIR   <- "mammalian_RI_dpsi01_fdr05"

PATHS <- list(
  Arabidopsis = file.path(BASE_DIR, "arabidopsis/nucleus/step12_protein_metagene"),
  Human       = file.path(BASE_DIR, "human/nucleus/step12_protein_metagene"),
  Mouse       = file.path(BASE_DIR, "mouse/total/step12_protein_metagene")
)

GROUPS <- c("Control", "Dependent", "Independent")
COLORS <- c(Control = "#9E9E9E", Dependent = "#E57373", Independent = "#64B5F6")

PROPERTIES <- c(
  "gravy", "kd_hydropathy", "aromaticity", "instability", 
  "isoelectric_point", "charge_at_pH7", "aliphatic_index",
  "flexibility", "helix_fraction", "sheet_fraction"
)

# Friendly names for plots
PROP_LABELS <- c(
  gravy = "GRAVY",
  kd_hydropathy = "Kyte-Doolittle",
  aromaticity = "Aromaticity",
  instability = "Instability",
  isoelectric_point = "pI",
  charge_at_pH7 = "Charge (pH 7)",
  aliphatic_index = "Aliphatic Idx",
  flexibility = "Flexibility",
  helix_fraction = "Alpha-Helix",
  sheet_fraction = "Beta-Sheet"
)

# Protein counts (N) derived from previous execution logs
PROTEIN_COUNTS <- list(
  Arabidopsis = c(Control = 360, Dependent = 128, Independent = 37),
  Human       = c(Control = 1399, Dependent = 489, Independent = 286),
  Mouse       = c(Control = 1394, Dependent = 159, Independent = 116)
)

# Start PDF
pdf(OUTPUT_PDF, width = A4_WIDTH, height = A4_HEIGHT)

# Layout: 10 rows (properties) x 3 columns (species)
par(mfcol = c(10, 3), 
    mar = c(1.5, 2.5, 1.2, 0.5), # bottom, left, top, right
    oma = c(3, 3, 5, 1), # outer margins
    mgp = c(1.2, 0.4, 0),
    tcl = -0.35, 
    bty = "n",
    cex.main = 0.8,
    cex.lab = 0.7,
    cex.axis = 0.6,
    family = "sans")

# Function to plot a single panel
plot_metagene <- function(property, data_list, species_name, row_idx, col_idx) {
  # Calculate global Y limits for this property across groups for the species
  y_min <- min(sapply(data_list, function(d) min(d[[property]] - d[[paste0(property, "_SEM")]], na.rm=TRUE)), na.rm=TRUE)
  y_max <- max(sapply(data_list, function(d) max(d[[property]] + d[[paste0(property, "_SEM")]], na.rm=TRUE)), na.rm=TRUE)
  
  # Ensure some padding
  y_range <- y_max - y_min
  if(y_range == 0) y_range <- 0.1
  ylim <- c(y_min - 0.05 * y_range, y_max + 0.15 * y_range)
  
  # Empty plot
  plot(NULL, xlim = c(0, 100), ylim = ylim, 
       main = if(row_idx == 1) "" else PROP_LABELS[property],
       xlab = "", ylab = "", axes = FALSE)
  
  # Custom header for top row
  if(row_idx == 1) {
    counts <- PROTEIN_COUNTS[[species_name]]
    header_text <- sprintf("%s\n(n: C=%d, D=%d, I=%d)", toupper(species_name), counts["Control"], counts["Dependent"], counts["Independent"])
    mtext(header_text, side = 3, line = 0.5, cex = 0.7, font = 2)
    # Put actual property label inside or above the first plot too
    mtext(PROP_LABELS[property], side = 3, line = -0.8, cex = 0.6, font = 1)
  }

  # Add axes (only left column for Y, only bottom row for X)
  axis(1, lwd = 0.5, labels = (row_idx == 10))
  axis(2, lwd = 0.5, labels = TRUE, las = 2)
  
  # Property label on leftmost plots as y-axis descriptor
  if(col_idx == 1) {
    mtext(PROP_LABELS[property], side = 2, line = 1.5, cex = 0.5, font = 2)
  }

  # Draw polygons for SEM
  for (grp in GROUPS[length(GROUPS):1]) { # Reversing order to put dependent/independent on top Visually
    df <- data_list[[grp]]
    if (is.null(df) || !(property %in% names(df))) next
    
    x <- df$Position_pct
    y <- df[[property]]
    sem <- df[[paste0(property, "_SEM")]]
    
    valid <- !is.na(y) & !is.na(sem)
    if(sum(valid) < 2) next
    
    polygon(c(x[valid], rev(x[valid])), 
            c(y[valid] + sem[valid], rev(y[valid] - sem[valid])),
            col = adjustcolor(COLORS[grp], alpha.f = 0.15), border = NA)
  }
  
  # Draw means
  for (grp in GROUPS) {
    df <- data_list[[grp]]
    if (is.null(df) || !(property %in% names(df))) next
    
    x <- df$Position_pct
    y <- df[[property]]
    valid <- !is.na(y)
    
    lines(x[valid], y[valid], col = COLORS[grp], lwd = 1)
  }
}

# Main Loop through Species (Columns)
for (col_idx in 1:length(SPECIES_ORDER)) {
  species <- SPECIES_ORDER[col_idx]
  
  # Load data for this species
  data_list <- list()
  for (grp in GROUPS) {
    file <- file.path(PATHS[[species]], paste0("metagene_", grp, ".tsv"))
    if (file.exists(file)) {
      data_list[[grp]] <- read.table(file, sep = "\t", header = TRUE)
    }
  }
  
  # Plot 10 properties (Rows)
  for (row_idx in 1:length(PROPERTIES)) {
    if (length(data_list) == 0) {
      plot.new()
    } else {
      plot_metagene(PROPERTIES[row_idx], data_list, species, row_idx, col_idx)
    }
  }
}

# Global Legend
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend = GROUPS, col = COLORS, lwd = 3, horiz = TRUE, bty = "n", cex = 0.8)

# Axis Labels
mtext("Protein Position (%)", side = 1, outer = TRUE, line = 1, cex = 1)
mtext("Biochemical Property Scores", side = 2, outer = TRUE, line = 1, cex = 1)
mtext("Protein Biochemical Properties along Metagene Sequence", side = 3, outer = TRUE, line = 2, cex = 1.2, font = 2)

dev.off()
cat("Figure generated: ", OUTPUT_PDF, "\n")

# Global Legend (outside margins)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend = GROUPS, col = COLORS, lwd = 3, horiz = TRUE, bty = "n", cex = 1)

# Super Title
mtext("Protein Biochemical Properties along Metagene Sequence", side = 3, outer = TRUE, line = -2, cex = 1.5, font = 2)

dev.off()
cat("Figure generated: ", OUTPUT_PDF, "\n")
