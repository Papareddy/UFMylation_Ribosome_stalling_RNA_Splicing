library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(purrr)
library(biomaRt)

library(stringr)
mart <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "https://plants.ensembl.org")
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), mart = mart)

res <- res %>%
  dplyr::left_join(gene_info, by = c("gene_id" = "ensembl_gene_id")) %>%
  dplyr::mutate(label_to_show = ifelse(is.na(external_gene_name) | external_gene_name == "", gene_id, external_gene_name)) %>%
  cbind(expression_values[rownames(res), , drop = FALSE])




df <- list.files(path = "~/Desktop/Desktop-Ranjith-iMac/Export_MS/Draft1_Oct_2025/Supplimentary_Tables/", pattern = "\\.txt$", full.names = TRUE) %>%
  map_dfr(read_tsv)

df=read_tsv("~/Desktop/Desktop-Ranjith-iMac/Export_MS/Draft1_Oct_2025/Supplimentary_Tables/hsUFM_Coevolved_proteins_AF2_MM.tsv")
df <- df %>%
  tidyr::separate(NAME, into = c("number", "bait", "vs", "prey"), sep = "_") %>%
  dplyr::select(-vs) 

library(dplyr)
library(stringr)

k   <- 3      # number of top models
thr <- 0.7   # threshold (still usable elsewhere)

mean_top_k_from_colon_str <- function(x, k) {
  vals <- as.numeric(str_split(x, ":", simplify = TRUE))
  mean(sort(vals, decreasing = TRUE)[seq_len(k)], na.rm = TRUE)
}

df_max <- df %>%
  rowwise() %>%
  mutate(
    max_scaled_PEAK        = max(as.numeric(str_split(scaled_PEAK, ":", simplify = TRUE)), na.rm = TRUE),
    mean_top3_scaled_PEAK  = mean_top_k_from_colon_str(scaled_PEAK, k)
  ) %>%
  ungroup()

write_tsv(df_max, "~/Desktop/Desktop-Ranjith-iMac/Export_MS/Draft1_Oct_2025/Supplimentary_Tables/hsUFM_Coevolved_proteins_AF2_MM.tsv")

write_tsv(df_max, "~/Desktop/Desktop-Ranjith-iMac/forGIT/UFMylation_Ribosome_stalling_RNA_Splicing//SuplimentalTables//hsUFM_Coevolved_proteins_AF2_MM.tsv")

bait_col <- c(
  "C53"   = "#d69547",
  "DDRGK1"= "#d96d83",
  "UBA5"  = "#29558b",
  "UFC1"  = "#a65c4c",
  "UFL1"  = "#810000",
  "UFM1full"  = "#149a91",
  "UFSP2" = "#556052"
)


bait_col <- c(
  "hsCDK5RAP3"   = "#d69547",
  "hsDDRGK1"= "#d96d83",
  "hsUBA5"  = "#29558b",
  "hsUFC1"  = "#a65c4c",
  "hsUFL1"  = "#810000",
  "HsUFM1"  = "#149a91"
)

library(ggplot2)
library(ggrepel)

# Custom theme function
my_plot_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      axis.ticks.length = unit(0.2, "cm"),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      panel.grid = element_blank()
    )
}

df_max$prey=gsub("\\..*", "", df_max$prey)
df_max

df_max <- df_max %>%
  dplyr::inner_join(gene_info, by = c("prey" = "ensembl_gene_id")) 


# Plot
#df_max <- df_max %>%filter(bait %in% c("UFM1full", "UFL1", "C53", "DDRGK1")) %>% mutate(bait = factor(bait, levels = c("UFM1full", "UFL1", "C53", "DDRGK1")))
top=df_max%>%filter(max_scaled_PEAK > 0.75)

ggplot(df_max, aes(x = IPTMavg, y = scaled_PEAKavg)) +
  geom_point(
    size = 2.5,color="white",fill="#149a91",
    shape = 21,
    alpha = 0.75
  ) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey75", alpha = 0.25) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "grey75", alpha = 0.25) +
  scale_size_continuous(range = c(1, 5)) +
 
 
  ylim(0.2, 1.2) +xlim(0,1.2)+
  coord_cartesian(clip = "off") +
  my_plot_theme(base_size = 12)+
    geom_text_repel(
     data = top,
     aes(label = prey),
     size = 2,
     segment.curvature = -0.1,
     segment.ncp = 3,
     max.overlaps = 50,
     segment.angle = 20,
     box.padding = 0.5
   )
 

# Add labels for selected points
top3 <- df_max %>%
  group_by(bait) %>%
  slice_max(order_by = max_scaled_PEAK, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  # always include specific prey genes
  bind_rows(
    df_max %>% filter(prey %in% c("AT4G24690", "AT4G27120"))
  ) %>%
  distinct()
 
top=df_max%>%filter(max_scaled_PEAK > 0.75)%>%arrange(prey)


df=read_tsv("~/Desktop/Desktop-Ranjith-iMac/Export_MS/Draft1_Oct_2025/Supplimentary_Tables/hsUFM_Coevolved_proteins_AF2_MM.tsv")
 
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

bm <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol", "description"),
  filters = "uniprotswissprot",
  values = top$prey,
  mart = mart
)

bm

 
 
 library(ggplot2)
 library(ggbeeswarm)
 library(ggrepel)
 
 ggplot(df_max, aes(x = "All proteins", y = scaled_PEAKavg)) +
   geom_beeswarm(aes(size = 3, fill = bait),shape = 21,color = "white", alpha = 1) +

   
   geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey75", alpha = 0.25) +
   scale_y_continuous(limits = c(0, 1.0)) +
   scale_x_discrete(expand = c(0.1, 0.1)) + 
   coord_cartesian(clip = "off") +
   my_plot_theme(base_size = 12) +
   labs(x = NULL, y = "IPTMavg") +facet_wrap(~bait)+
   geom_text_repel(
     data = top,
     aes(x = "All proteins", y = max_scaled_PEAK, label = prey),
     size = 2.5,
     max.overlaps = 50,
     box.padding = 0.5,
     force = 2,             # stronger push away from points
     force_pull = 1,        # pulls labels radially outward
     segment.curvature = -0.1,
     segment.ncp = 3,
     segment.angle = 20
   )
 
 
 ggplot(df_max, aes(x = 1, y = max_scaled_PEAK)) +
   geom_beeswarm(
     size = 4, shape = 21,
     color = "white", fill = "#149a91", alpha = 0.75
   ) +
   geom_text_repel(
     data = top3,
     aes(x = 1, y = max_scaled_PEAK, label = external_gene_name),
     size = 2.5,
     max.overlaps = Inf,
     box.padding = 0.5,
     force_pull = 1
   ) +
   coord_polar(theta = "y") +
   theme_minimal() +
   labs(y = "IPTMavg")
 
 
 
 
 
 library(biomaRt)
 k   <- 3      # number of top models
 thr <- 0.7   # threshold (still usable elsewhere)
 
 mean_top_k_from_colon_str <- function(x, k) {
   vals <- as.numeric(str_split(x, ":", simplify = TRUE))
   mean(sort(vals, decreasing = TRUE)[seq_len(k)], na.rm = TRUE)
 }
 
 df_max <- df %>%
   rowwise() %>%
   mutate(
     max_scaled_PEAK        = max(as.numeric(str_split(scaled_PEAK, ":", simplify = TRUE)), na.rm = TRUE),
     mean_top3_scaled_PEAK  = mean_top_k_from_colon_str(scaled_PEAK, k)
   ) %>%
   ungroup()
 
 top <- df_max %>%
   filter(
     mean_top3_scaled_PEAK > 0.8
   )
 ids <- unique(top$prey)
 
 mart_hs <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
 
 
 bm <- getBM(
   attributes = c("uniprotswissprot", "external_gene_name", "description"),
   filters    = "uniprotswissprot",
   values     = ids,
   mart       = mart_hs
 )
 
 bm
 
 top_annot <- top %>%
   left_join(bm, by = c("prey" = "uniprotswissprot"))
 
 
 library(dplyr)
 library(stringr)
 
 # --- 1) Build gene sets (explicit lists + prefix rules) ---
 genes_rna_degradation <- c("SMG9", "ZCCHC7","MEX3A")
 
 genes_dna_damage <- c(
   "FANCI","FANCL","FANCD2","INTS7","INTS3","NABP2","NABP1",
   "NHEJ1","TREX1","TREX2","CCAR2"
 )
 
 genes_chromatin <- c("SETD7")
 
 genes_rna_processing_explicit <- c(
   "RNPC3","RBM41","NUP153","THOC6","PNN","ZNF830","TTF2",
   "GPATCH8","CSTF2","CSTF2T","ZCCHC13"
 )
 
 genes_ufmylation <- c(
   "UBA5","UFC1","UFL1","UFM1","UFSP1","UFSP2","ODR4","DDRGK1","CDK5RAP3"
 )
 
 # Prefix-based membership check
 is_pref <- function(x, pref) str_detect(toupper(x), paste0("^", toupper(pref)))
 
 # --- 2) Categorize and filter top_annot by external_gene_name ---
 top_annot_cat <- top_annot %>%
   mutate(
     gene_name = external_gene_name,
     category = case_when(
       is_pref(gene_name, "ERI") | gene_name %in% genes_rna_degradation ~ "RNA degradation",
       is_pref(gene_name, "FANC") | gene_name %in% genes_dna_damage     ~ "DNA damage",
       is_pref(gene_name, "PRDM") | is_pref(gene_name, "PCGF") | gene_name %in% genes_chromatin ~ "Chromatin modifiers",
       is_pref(gene_name, "RBM")  | is_pref(gene_name, "ERI") |
         gene_name %in% c(genes_rna_processing_explicit, genes_rna_degradation) ~ "RNA processing",
       gene_name %in% genes_ufmylation ~ "UFMylation",
       TRUE ~ NA_character_
     )
   ) %>%
   filter(!is.na(category))
 

 
 
 
 ggplot(df_max, aes(x = "All proteins", y = mean_top3_scaled_PEAK)) +
   geom_beeswarm(aes(size = 3, fill = bait),shape = 21,color = "white", alpha = 1) +
   
   
   geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey75", alpha = 0.5) +
   scale_y_continuous(limits = c(0, 1.0)) +
   scale_x_discrete(expand = c(0.1, 0.1)) + 
   coord_cartesian(clip = "off") +
   my_plot_theme(base_size = 12) +
   labs(x = NULL, y = "IPTMavg") +facet_wrap(~bait)+
   geom_text_repel(
     data = top_annot_cat,
     aes(x = "All proteins", y = max_scaled_PEAK, label = external_gene_name),
     size = 2.5,
     max.overlaps = 50,
     box.padding = 0.5,
     force = 2,             # stronger push away from points
     force_pull = 1,        # pulls labels radially outward
     segment.curvature = -0.1,
     segment.ncp = 3,
     segment.angle = 20
   )
 