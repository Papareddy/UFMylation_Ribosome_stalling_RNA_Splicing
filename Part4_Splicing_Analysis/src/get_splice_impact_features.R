suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(biomaRt)
  library(biomaRt)
  library(stats)
  library(digest)
})

# --- Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (grepl("^--", arg)) {
      parts <- strsplit(sub("^--", "", arg), "=")[[1]]
      if (length(parts) == 2) out[[parts[1]]] <- parts[2]
    }
  }
  return(out)
}
opt <- parse_args(args)

if (is.null(opt$gtf)) opt$gtf <- "data/human/pcg_gencode.v45.annotation.gtf.gz"

# --- Main Logic ---

# 1. Load Data
lost_df <- read_tsv(opt$lost, show_col_types = FALSE)
preserved_df <- read_tsv(opt$preserved, show_col_types = FALSE)

# 2. Prepare TxDb
# 2. Prepare TxDb (with Caching)
message("Preparing TxDb from GTF...")
cache_dir <- if (!is.null(opt$cache_dir)) opt$cache_dir else "data/cache"
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

txdb_cache_file <- file.path(cache_dir, paste0(basename(opt$gtf), ".sqlite"))

if (file.exists(txdb_cache_file)) {
    message("Loading cached TxDb: ", txdb_cache_file)
    txdb <- suppressWarnings(AnnotationDbi::loadDb(txdb_cache_file))
} else {
    message("Creating TxDb from GTF (this may take a while)...")
    txdb <- suppressWarnings(makeTxDbFromGFF(opt$gtf, format = "gtf"))
    AnnotationDbi::saveDb(txdb, txdb_cache_file)
    message("Saved cached TxDb to ", txdb_cache_file)
}

# 3. Setup biomaRt
# 3. Setup biomaRt
message("Connecting to biomaRt...")
dataset <- "hsapiens_gene_ensembl"
mart_name <- "genes"
host_url <- NULL # Default for main Ensembl

if (!is.null(opt$species)) {
  if (opt$species == "mouse") {
    dataset <- "mmusculus_gene_ensembl"
  } else if (opt$species == "arabidopsis") {
    dataset <- "athaliana_eg_gene"
    mart_name <- "plants_mart"
    host_url <- "https://plants.ensembl.org"
  }
}
message("Using Ensembl dataset: ", dataset, " (mart: ", mart_name, ")")

connect_biomart <- function(dataset, mart, host) {
    if (!is.null(host)) {
        # Ensembl Genomes (Plants, etc.)
        tryCatch({
            message("Connecting to ", host, "...")
            return(useMart(biomart = mart, dataset = dataset, host = host))
        }, error = function(e) { 
             message("Connection failed: ", e$message)
             stop("Could not connect to Ensembl Plants.")
        })
    } else {
        # Main Ensembl (Human, Mouse)
        # Try 1: Default
        tryCatch({
            message("Attempt 1: Connecting to main Ensembl site...")
            return(useEnsembl(biomart = mart, dataset = dataset, mirror = "www"))
        }, error = function(e) { message("Attempt 1 failed: ", e$message) })

        # Try 2: US East
        tryCatch({
            message("Attempt 2: Connecting to US East mirror...")
            return(useEnsembl(biomart = mart, dataset = dataset, mirror = "useast"))
        }, error = function(e) { message("Attempt 2 failed: ", e$message) })
        
        # Try 3: Asia
        tryCatch({
            message("Attempt 3: Connecting to Asia mirror...")
            return(useEnsembl(biomart = mart, dataset = dataset, mirror = "asia"))
        }, error = function(e) { message("Attempt 3 failed: ", e$message) })

        # Try 4: Ensembl Version specific (occasionally helps)
        tryCatch({
            message("Attempt 4: Connecting to Ensembl (no mirror spec)...")
            return(useEnsembl(biomart = mart, dataset = dataset))
        }, error = function(e) { message("Attempt 4 failed: ", e$message) })

        stop("Could not connect to any BioMart mirror after 4 attempts.")
    }
}

ensembl <- connect_biomart(dataset, mart_name, host_url)

# --- Helper: Map Events to Genes ---
get_gene_ids_from_events <- function(df, txdb) {
  if (nrow(df) == 0) return(character(0))
  
  # Construct GRanges
  get_col_name <- function(base, df_names) {
    if (base %in% df_names) return(base)
    if (paste0(base, ".WT") %in% df_names) return(paste0(base, ".WT"))
    return(NA)
  }
  
  starts <- numeric(nrow(df))
  ends   <- numeric(nrow(df))
  chrs   <- character(nrow(df))
  strands <- character(nrow(df))
  
  if (!"EventType" %in% names(df)) return(character(0))
  
  unique_types <- unique(df$EventType)
  for (ev in unique_types) {
     idx <- which(df$EventType == ev)
     c_s <- NA; c_e <- NA
     if (ev == "SE") { c_s <- get_col_name("exonStart_0base", names(df)); c_e <- get_col_name("exonEnd", names(df)) }
     else if (ev == "RI") { c_s <- get_col_name("riExonStart_0base", names(df)); c_e <- get_col_name("riExonEnd", names(df)) }
     else if (ev %in% c("A3SS", "A5SS")) { c_s <- get_col_name("longExonStart_0base", names(df)); c_e <- get_col_name("longExonEnd", names(df)) }
     else if (ev == "MXE") { c_s <- get_col_name("1stExonStart_0base", names(df)); c_e <- get_col_name("1stExonEnd", names(df)) }
     c_chr <- get_col_name("chr", names(df))
     c_str <- get_col_name("strand", names(df))
     
     if (!is.na(c_s) && !is.na(c_e)) {
         starts[idx] <- df[[c_s]][idx]
         ends[idx]   <- df[[c_e]][idx]
         chrs[idx] <- df[[c_chr]][idx]
         strands[idx] <- df[[c_str]][idx]
     }
  }
  
  valid_mask <- !is.na(starts)
  if (sum(valid_mask) == 0) return(character(0))
  
  # Suppress deprecation warnings from S4Vectors during GRanges construction/seqlevels manipulation
  suppressWarnings({
      gr <- GRanges(seqnames = chrs[valid_mask], ranges = IRanges(starts[valid_mask] + 1, ends[valid_mask]), strand = strands[valid_mask])
      
      # Fix seqlevels
      target_style <- seqlevelsStyle(txdb)
      data_style <- seqlevelsStyle(gr)
      if (!any(data_style %in% target_style)) {
         tryCatch({ seqlevelsStyle(gr) <- target_style[1] }, error = function(e){})
      }
    
      genes_gr <- genes(txdb)
      ov <- findOverlaps(gr, genes_gr)
  })
  return(unique(names(genes_gr)[subjectHits(ov)]))
}

message("Extracting gene lists...")
lost_genes <- get_gene_ids_from_events(lost_df, txdb)
preserved_genes <- get_gene_ids_from_events(preserved_df, txdb)
background_genes <- names(genes(txdb)) # All genes in annotation

# Clean IDs (remove .version)
clean_id <- function(x) sub("\\..*", "", x)
lost_genes_clean <- unique(clean_id(lost_genes))
preserved_genes_clean <- unique(clean_id(preserved_genes))
bg_genes_clean <- unique(clean_id(background_genes))

message("Found Genes: Lost=", length(lost_genes_clean), " Preserved=", length(preserved_genes_clean), " Background=", length(bg_genes_clean))

# --- Fetch Domains for Background ---
# --- Load Custom Background if provided (Move to top for consistent hashing/fetching) ---
if ("background_genes" %in% names(opt)) {
  message("Using custom background genes from: ", opt$background_genes)
  bg_genes_custom <- readLines(opt$background_genes)
  bg_genes_clean <- clean_id(bg_genes_custom)
} else {
  message("Using default genomic background (all genes in TxDb).")
}

# --- Fetch Domains for Background ---
message("Fetching background domains...")

# Create content-aware cache key based on the gene list
bg_hash <- digest::digest(sort(bg_genes_clean), algo="md5")
domain_cache_file <- file.path(cache_dir, paste0("domains_bg_", dataset, "_", bg_hash, ".rds"))
biophys_cache_file <- file.path(cache_dir, paste0("biophys_bg_", dataset, "_", bg_hash, ".rds"))

if (file.exists(domain_cache_file)) {
    message("Loading cached domains from: ", domain_cache_file)
    bg_domains_df <- readRDS(domain_cache_file)
} else {
    message("Fetching background domains from Biomart (this may take a while)...")
    get_domains_bulk <- function(gene_ids, mart) {
        res <- data.frame(ensembl_gene_id=character(), interpro_short_description=character())
        chunk_size <- 500
        ids <- unique(gene_ids)
        for(i in seq(1, length(ids), by=chunk_size)) {
            chunk <- ids[i : min(i+chunk_size-1, length(ids))]
            tryCatch({
                bm <- getBM(attributes = c("ensembl_gene_id", "interpro_short_description"), 
                            filters = "ensembl_gene_id", 
                            values = chunk, 
                            mart = mart)
                res <- rbind(res, bm)
            }, error = function(e) { message("Chunk error: ", e$message)})
            if (i %% 10000 == 1) message("Fetched ", i, "/", length(ids))
        }
        return(res[res$interpro_short_description != "", ])
    }

    bg_domains_df <- get_domains_bulk(bg_genes_clean, ensembl)
    saveRDS(bg_domains_df, domain_cache_file)
    message("Saved cached domains to ", domain_cache_file)
}

# --- Enrichment Analysis Function ---
calc_enrichment <- function(target_genes, bg_genes, domain_df, set_name) {
    # Universe: Background Genes
    # Target: Set Genes
    
    # Map Check
    # Total Universe Count
    N_bg <- length(bg_genes)
    N_target <- length(target_genes)
    N_complement <- N_bg - N_target
    
    # Domains present in target
    target_doms_df <- domain_df[domain_df$ensembl_gene_id %in% target_genes, ]
    dom_counts_target <- table(target_doms_df$interpro_short_description)
    
    # Domains present in background
    # Note: domain_df encompasses background (as we queried bg_genes)
    dom_counts_bg <- table(domain_df$interpro_short_description)
    
    results <- data.frame()
    
    # Test only domains present in the target set (for enrichment)
    # Can also test all domains for strict depletion, but simpler loop over observed
    test_domains <- names(dom_counts_target)
    
    for (d in test_domains) {
        # a = Genes in Target with Domain
        # b = Genes in Target without Domain
        # c = Genes in Complement with Domain
        # d = Genes in Complement without Domain
        
        # Count unique genes with domain d
        genes_with_d <- unique(domain_df$ensembl_gene_id[domain_df$interpro_short_description == d])
        
        a <- length(intersect(genes_with_d, target_genes))
        b <- N_target - a
        
        # c = (Total genes with D) - a
        c_val <- length(intersect(genes_with_d, setdiff(bg_genes, target_genes)))
        d_val <- N_complement - c_val
        
        # Fisher Test
        mat <- matrix(c(a, c_val, b, d_val), nrow=2)
        ft <- fisher.test(mat, alternative = "two.sided")
        
        results <- rbind(results, data.frame(
            domain = d,
            pval = ft$p.value,
            odds_ratio = ft$estimate,
            log2_odds = log2(ft$estimate),
            count_in_set = a,
            count_in_bg = length(genes_with_d),
            set = set_name
        ))
    }
    
    if (nrow(results) > 0) {
        results$fdr <- p.adjust(results$pval, method = "BH")
    }
    
    return(results)
}




# --- Enrichment Analysis Wrapper ---
run_enrichment_sets <- function(sets_list, bg_genes, domain_df) {
    res_list <- list()
    for (set_name in names(sets_list)) {
        message("Calculating enrichment for ", set_name, "...")
        genes <- sets_list[[set_name]]
        if (length(genes) > 0) {
            res_list[[set_name]] <- calc_enrichment(genes, bg_genes, domain_df, set_name)
        }
    }
    return(do.call(rbind, res_list))
}

# --- Define Sets (Directional or Not) ---
sets_to_analyze <- list()
direction_enabled <- !is.null(opt$direction) && opt$direction == "TRUE" # Passed as string?

# Check if flag is present in opt (it's a flag, so if parse_args handles flags without value, 
# but run_pipeline probably passes --direction without value?
# Wait, parse_args logic split only by =. If run_pipeline passes just --direction, 
# my parse_args might fail or put it as key without value.
# Let's check parse_args implementation above.
# It does: parts <- strsplit(sub("^--", "", arg), "=")[[1]]; if (length(parts) == 2) ...
# So boolean flags like --direction won't be in `opt` list if passed as just --direction !
# I need to fix arg parsing or pass --direction=TRUE.
# Safest: check raw 'args' for "--direction".

# Safest check: Use the parsed opt value
use_direction <- FALSE
if (!is.null(opt$direction)) {
  if (opt$direction == "TRUE") use_direction <- TRUE
}

if (use_direction) {
    message("Directionality enabled: Splitting by dPSI...")
    
    # Helper to get dPSI column. Usually dPSI_num.WT or similar.
    # We need to look at column names.
    get_dpsi <- function(df) {
        if ("dPSI_num.WT" %in% names(df)) return(df$dPSI_num.WT)
        if ("dPSI_num" %in% names(df)) return(df$dPSI_num)
        return(rep(0, nrow(df))) # Fallback
    }
    
    # Function to clean genes for a subset
    get_genes_subset <- function(df, condition_idx) {
        sub_df <- df[condition_idx, ]
        g <- get_gene_ids_from_events(sub_df, txdb)
        return(unique(clean_id(g)))
    }
    
    # Lost
    dpsi_l <- get_dpsi(lost_df)
    sets_to_analyze[["Lost_psI"]] <- get_genes_subset(lost_df, dpsi_l > 0)
    sets_to_analyze[["Lost_psE"]] <- get_genes_subset(lost_df, dpsi_l < 0)
    
    # Preserved
    dpsi_p <- get_dpsi(preserved_df)
    sets_to_analyze[["Preserved_psI"]] <- get_genes_subset(preserved_df, dpsi_p > 0)
    sets_to_analyze[["Preserved_psE"]] <- get_genes_subset(preserved_df, dpsi_p < 0)
    
} else {
    sets_to_analyze[["Lost"]] <- lost_genes_clean
    sets_to_analyze[["Preserved"]] <- preserved_genes_clean
}

final_res <- run_enrichment_sets(sets_to_analyze, bg_genes_clean, bg_domains_df)

out_file <- file.path(opt$outdir, "domain_enrichment.tsv")
write_tsv(final_res, out_file)
message("Saved enrichment results to ", out_file)


# Also create the annotated files as requested previously (combining logic)
# (Re-using parts of previous logic for file output if needed, but primary goal here is enrichment)
# To preserve the 'impact_domains' column in the TSVs, we can quickly map back.
# But for now, focusing on the Enrichment TSV output.


# --- Biophysical Feature Enrichment (SignalP, TMHMM, NCOILS) ---
# Function defined globally for clarity/reuse
get_biophys_bulk <- function(gene_ids, mart) {
    # chunking
    res <- data.frame(ensembl_gene_id=character(), tmhmm_start=integer(), signalp_start=integer(), ncoils_start=integer())
    chunk_size <- 500
    ids <- unique(gene_ids)
    
    # Attributes: using start positions as proxy for existence
    attrs <- c("ensembl_gene_id", "tmhmm_start", "signalp_start", "ncoils_start")
    
    for(i in seq(1, length(ids), by=chunk_size)) {
        chunk <- ids[i : min(i+chunk_size-1, length(ids))]
        tryCatch({
            bm <- getBM(attributes = attrs, 
                        filters = "ensembl_gene_id", 
                        values = chunk, 
                        mart = mart)
            res <- rbind(res, bm)
        }, error = function(e) { message("Biophys chunk error: ", e$message) })
        if (i %% 10000 == 1) message("Fetched biophys ", i, "/", length(ids))
    }
    return(res)
}

if (file.exists(biophys_cache_file)) {
    message("Loading cached biophysical attributes from: ", biophys_cache_file)
    biophys_raw <- readRDS(biophys_cache_file)
} else {
    message("Fetching biophysical features (SignalP, TMHMM, NCOILS)...")
    biophys_raw <- get_biophys_bulk(bg_genes_clean, ensembl)
    saveRDS(biophys_raw, biophys_cache_file)
    message("Saved cached biophysical attributes to ", biophys_cache_file)
}

# (Removed old custom background check block from here as it is moved up)

# Process into binary flags per gene
# If a gene appears in the result with a non-NA start value for a feature, it has that feature.
# Note: getBM might return NA or just omit rows.
# Aggregating by gene
# Handle empty biophys_raw
if (nrow(biophys_raw) == 0) {
    message("[WARN] No biophysical features fetched. Creating dummy flags.")
    biophys_flags <- data.frame(ensembl_gene_id = character(), has_tmhmm = logical(), has_signalp = logical(), has_ncoils = logical())
} else {
    biophys_flags <- biophys_raw %>%
        group_by(ensembl_gene_id) %>%
        summarise(
            has_tmhmm = any(!is.na(tmhmm_start)),
            has_signalp = any(!is.na(signalp_start)),
            has_ncoils = any(!is.na(ncoils_start))
        )
}

# Ensure all background genes are present (defaults to FALSE if not found in BM result but in bg list?)
# Actually, if not in BM result, we might lack info or it lacks features. Assume lacks features.
biophys_bg <- data.frame(ensembl_gene_id = bg_genes_clean) %>%
    left_join(biophys_flags, by = "ensembl_gene_id") %>%
    mutate(
        has_tmhmm = ifelse(is.na(has_tmhmm), FALSE, has_tmhmm),
        has_signalp = ifelse(is.na(has_signalp), FALSE, has_signalp),
        has_ncoils = ifelse(is.na(has_ncoils), FALSE, has_ncoils)
    )

calc_feature_enrichment <- function(target_ids, bg_df, feature_col, label) {
    # Target in Background
    # We restrict background to genes where we have info (which is all bg_df)
    
    # Counts
    # q = Genes in Target with Feature
    # m = Total Genes in Background with Feature
    # n = Total Genes in Background WITHOUT Feature
    # k = Size of Target
    
    # Fisher matrix:
    #      HasFeat  NoFeat
    # InSet    a       b
    # NotInSet c       d
    
    in_set <- bg_df$ensembl_gene_id %in% target_ids
    has_feat <- bg_df[[feature_col]]
    
    tbl <- table(factor(in_set, levels=c(TRUE, FALSE)), factor(has_feat, levels=c(TRUE, FALSE)))
    # tbl:
    #       TRUE (Has)   FALSE (No)
    # TRUE     a            b
    # FALSE    c            d
    
    ft <- fisher.test(tbl, alternative = "two.sided")
    
    return(data.frame(
        feature = label,
        odds_ratio = ft$estimate,
        pval = ft$p.value,
        log2_odds = log2(ft$estimate),
        count_in_set = tbl[1,1],
        set_size = sum(in_set),
        bg_count = sum(has_feat),
        bg_size = nrow(bg_df)
    ))
}

features <- c("has_signalp", "has_tmhmm", "has_ncoils")
labels <- c("SignalP", "TMHMM", "NCOILS")

res_biophys <- data.frame()


# Start Biophysical Analysis
# User Request: Directionality should NOT apply to Protein Attributes.
# We always analyze the lumped "Lost" and "Preserved" sets here, regardless of opt$direction.

biophys_sets <- list()
biophys_sets[["Lost"]] <- lost_genes_clean
biophys_sets[["Preserved"]] <- preserved_genes_clean

# Loop over defined sets (Always Lost/Preserved)
for(set_name in names(biophys_sets)) {
    target_ids <- biophys_sets[[set_name]]
    if(length(target_ids) == 0) next
    
    # Compare vs Genome
    for(i in 1:3) {
        r <- calc_feature_enrichment(target_ids, biophys_bg, features[i], labels[i])
        r$comparison <- paste0(set_name, "_vs_Genome")
        res_biophys <- rbind(res_biophys, r)
    }
}

# Optional: Direct comparisons (e.g. Lost_Inc vs Preserved_Inc) can be complex to automate generically.
# For now, we stick to vs Genome for all subgroups.
# If direction is NOT enabled, we can keeping the old "Lost vs Preserved" comparison for backward compatibility?
# Optional: Direct comparison Lost vs Preserved (Always done now since sets are fixed)
pop_lp <- biophys_bg[biophys_bg$ensembl_gene_id %in% c(lost_genes_clean, preserved_genes_clean), ]
for(i in 1:3) {
    r <- calc_feature_enrichment(lost_genes_clean, pop_lp, features[i], labels[i])
    r$comparison <- "Lost_vs_Preserved"
    res_biophys <- rbind(res_biophys, r)
}

res_biophys$fdr <- p.adjust(res_biophys$pval, method = "BH")

write_tsv(res_biophys, file.path(opt$outdir, "biophysical_enrichment.tsv"))
message("Saved biophysical enrichment to ", file.path(opt$outdir, "biophysical_enrichment.tsv"))

message("Done.")
