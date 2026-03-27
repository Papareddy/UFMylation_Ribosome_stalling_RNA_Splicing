
# ==============================================================================
# MANUAL PLOTTING SCRIPT: Eukaryotic Conserved Splicing Landscape
# ==============================================================================
# Use this script to recreate or customize the bubble plots from the 
# cross-species analysis (Human, Mouse, Arabidopsis).

library(ggplot2)
library(dplyr)
library(readr)

# --- 1. SET DATA PATHS ---
# Relative paths from project root. Change these if your files are elsewhere.
DEPENDENT_DATA   <- "GO_Final_test/cross_species_visualization/Eukaryotic_Conserved_Landscape_Data.tsv"
INDEPENDENT_DATA <- "GO_Final_test/cross_species_visualization/Eukaryotic_Conserved_Independent_Landscape_Data.tsv"

# --- 2. DEFINE PLOTTING FUNCTION ---
plot_eukaryotic_landscape <- function(file_path, title_text, min_count = 1, min_species = 1) {
  if(!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  # Load Data
  df <- read_tsv(file_path, show_col_types = FALSE)
  
  # Exclude rRNA terms as requested
  df <- df %>% filter(!grepl("rRNA", Description, ignore.case=TRUE))

  # Apply Count Filter
  df <- df %>% filter(Count >= min_count)
  
  # Apply Min Species Filter
  # We identify Labels that appear in at least X species
  valid_labels <- df %>%
    group_by(Label) %>%
    summarise(n_spec = n_distinct(Species)) %>%
    filter(n_spec >= min_species) %>%
    pull(Label)
  
  df <- df %>% filter(Label %in% valid_labels)
  
  if(nrow(df) == 0) {
    stop("No data left after filtering! Try lowering min_count or min_species.")
  }

  # Ensure Species order is consistent
  df$Species <- factor(df$Species, levels = c("Human", "Mouse", "Arabidopsis"))
  df$Label <- factor(df$Label, levels = rev(unique(df$Label)))

  # Generate Plot
  p <- ggplot(df, aes(x = Species, y = Label)) +
    geom_point(aes(size = Count, color = logP)) +
    scale_color_viridis_c(option = "viridis", name = "-log10(p-value)", na.value = "transparent", end = 0.9) +
    scale_size_continuous(name = "Gene Count") +
    theme_bw() +
    labs(title = title_text,
         subtitle = paste0("Filters: Count >= ", min_count, ", Min Species = ", min_species),
         x = "", y = "") +
    theme(
      axis.text.y = element_text(size = 9, family = "sans"),
      axis.text.x = element_text(size = 12, face = "bold"),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray95", linetype = "dashed")
    )
  
  return(p)
}

# --- 3. EXECUTE & SHOW ---

# EXAMPLE: Plot with filters (e.g., minimum 5 genes and present in at least 2 species)
p_dep <- plot_eukaryotic_landscape(DEPENDENT_DATA, 
                                   "Filtered Eukaryotic Conserved Dependent Landscape",
                                   min_count = 5, 
                                   min_species = 2)
print(p_dep)

# Plot UFM1-Independent (Core Splicing Maintenance)
# p_indep <- plot_eukaryotic_landscape(INDEPENDENT_DATA, "Eukaryotic Conserved Independent Landscape")
# print(p_indep)

# --- 4. OPTIONAL: SAVE CUSTOM VERSION ---
# ggsave("My_Custom_Conserved_Landscape.pdf", p_dep, width = 10, height = 8)
