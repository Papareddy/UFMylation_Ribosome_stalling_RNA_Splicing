
library(MotifDb)

# 1. Fetch Motifs
print("Querying MotifDb for Hsapiens...")
all_human <- query(MotifDb, "Hsapiens")

# 2. Filter for RBP sources
# Common RBP sources: "CisBP-RNA", "Ray2013", "ATtRACT", "oRNAment", "RBPDB"
# Check what's available
sources <- table(values(all_human)$dataSource)
print("Available sources:")
print(sources)

rbp_sources <- c("Ray2013", "cisbp_1.02", "CisBP-RNA", "ATtRACT", "oRNAment", "RBPDB")
# Note: "cisbp_1.02" is likely DNA. "CisBP-RNA" is distinct.
# Let's inspect sources table to get exact names if possible.
# But for now, let's filter by regex "RNA" or specific names.

# Grep for RNA in dataSource or organism?
# Usually RBP motifs explicitly state RNA or come from these DBs.
# Let's try to get Ray2013 and CisBP-RNA specifically.
target_sources <- c("Ray2013", "CisBP-RNA")
motifs <- split(all_human, values(all_human)$dataSource)

final_list <- list()

if ("Ray2013" %in% names(motifs)) {
    print(paste("Adding", length(motifs[["Ray2013"]]), "Ray2013 motifs"))
    final_list <- c(final_list, as.list(motifs[["Ray2013"]]))
}

if ("CisBP-RNA" %in% names(motifs)) {
    # Verify exact name in next run, but assuming standard
     print(paste("Adding", length(motifs[["CisBP-RNA"]]), "CisBP-RNA motifs"))
     final_list <- c(final_list, as.list(motifs[["CisBP-RNA"]]))
}

# If no explict RBP sources found, we might need to look closer.
# But assuming they exist (Ray2013 is standard in MotifDb).

if (length(final_list) == 0) {
    # Regex fallback
    indicators <- grep("RNA", names(motifs), value=TRUE)
    for (src in indicators) {
         print(paste("Adding", length(motifs[[src]]), "motifs from", src))
         final_list <- c(final_list, as.list(motifs[[src]]))
    }
}

print(paste("Total RBP motifs collected:", length(final_list)))

# 3. Write to MEME format manually
write_meme_manual <- function(motif_list, filename) {
    con <- file(filename, "w")
    writeLines("MEME version 4", con)
    writeLines("", con)
    writeLines("ALPHABET= ACGU", con) # RBP motifs
    writeLines("", con)
    writeLines("strands: + -", con)
    writeLines("", con)
    writeLines("Background letter frequencies", con)
    writeLines("A 0.25 C 0.25 G 0.25 U 0.25", con)
    writeLines("", con)
    
    for (i in seq_along(motif_list)) {
        m <- motif_list[[i]]
        # m is a matrix (4 x width) usually A C G T
        # Need to check rownames. If DNA (T), convert to U for RNA context?
        # AME allows T in DNA alphabet mode, but if we say ACGU...
        # Let's check headers.
        
        # MotifDb names are like "Hsapiens-Ray2013-RNCMPT00123-SRSF1"
        full_name <- names(motif_list)[i]
        # Clean ID
        id <- gsub("[:; ]", "_", full_name) 
        
        # Transpose to rows=positions, cols=bases (MEME standard)
        # S4 matrix from MotifDb is typically 4 rows (ACGT) and N cols.
        pwm <- t(m) 
        width <- nrow(pwm)
        
        writeLines(paste("MOTIF", id, full_name), con)
        writeLines(paste("letter-probability matrix: alength= 4 w=", width, "nsites= 20 E= 0"), con)
        
        # Write rows
        for (r in 1:width) {
            # Ensure order A C G U (or T)
            # MotifDb rows are sorted? usually yes.
            probs <- pwm[r, ]
            # Handle T -> U output if needed, but numbers are same.
            # Output format: probA probC probG probT/U
            # Check column names of m (or now colnames(pwm) which are rownames(m))
            # If they are not named, assume ACGT.
            
            line <- paste(formatC(probs, format="f", digits=6), collapse=" ")
            writeLines(line, con)
        }
        writeLines("", con)
    }
    
    close(con)
    print(paste("Written to", filename))
}

if (length(final_list) > 0) {
    outfile <- "data/motifs/Human_RBPs_MotifDb.meme"
    write_meme_manual(final_list, outfile)
} else {
    print("No RBP motifs found to export.")
}
