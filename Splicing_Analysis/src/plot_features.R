suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

# Get args
args <- commandArgs(trailingOnly = TRUE)
# Expected: --input=master_file.tsv --out=out.pdf
input_file <- NA
out_file <- "features.pdf"

for(arg in args){
    if(startsWith(arg, "--input=")) input_file <- sub("^--input=", "", arg)
    if(startsWith(arg, "--out=")) out_file <- sub("^--out=", "", arg)
}

if(is.na(input_file)) stop("Usage: Rscript plot_features.R --input=master.tsv --out=out.pdf")

df <- read.delim(input_file, stringsAsFactors=FALSE)

# Ensure groups are factors
df$Group <- factor(df$Group, levels=c("Constitutive", "UFM1_independent", "UFM1_dependent"))

# Features to plot
features <- c("Length", "GC", "MaxEnt5", "MaxEnt3")

pdf(out_file, width=12, height=12)
par(mfrow=c(4,4), las=2, tcl=-0.3, bty="n", mar=c(4,4,2,1))

# For each feature (row), plot all groups?
# User requested: "For each metric, plot the three groups together ... panel PDF ... par(mfrow=c(4,4))"
# 4x4 grid suggests 16 plots? But we only have 4 metrics?
# Maybe the user implies 4 rows, 4 columns?
# Or maybe just use the parameter settings for style?
# "Visualize: Create a single multi-panel PDF. For each metric, plot the three groups together... Use the R parameters: par(mfrow=c(4,4)..."
# If I only have 4 metrics, I will fill 4 slots.

cols <- c("Constitutive"="grey", "UFM1_independent"="#66C2A5", "UFM1_dependent"="#FC8D62")

for (feat in features) {
    if (!feat %in% names(df)) next
    
    vals <- df[[feat]]
    grps <- df$Group
    
    # 1. Boxplot
    boxplot(vals ~ grps, col=cols[levels(grps)], main=paste0(feat, " Distribution"), ylab=feat)
    
    # 2. Density
    # Plot first density, add lines
    d_list <- split(vals, grps)
    
    # Get range
    x_lims <- range(vals, na.rm=TRUE)
    y_max <- 0
    densities <- list()
    for(g in names(d_list)){
        if(length(d_list[[g]]) > 1) {
            d <- density(d_list[[g]], na.rm=TRUE)
            densities[[g]] <- d
            y_max <- max(y_max, d$y)
        }
    }
    
    plot(NA, xlim=x_lims, ylim=c(0, y_max), xlab=feat, ylab="Density", main=paste0(feat, " Density"))
    for(g in names(densities)){
        lines(densities[[g]], col=cols[g], lwd=2)
    }
    legend("topright", legend=names(densities), col=cols[names(densities)], lwd=2, bty="n", cex=0.8)
    
    # 3. Barplot of Means? Or just skip to fill grid?
    # User asked for "each metric... overlapping density plots or grouped boxplots"
    # I have done both.
    
    # Let's add Statistical Comparison P-values as text plot?
    plot(0,0, type="n", axes=FALSE, xlab="", ylab="")
    
    # Calc stats
    # Safe format function
    safe_fmt <- function(val) {
        if(is.null(val) || is.na(val) || !is.numeric(val)) return("NA")
        formatC(val, format="e", digits=2)
    }

    res_l_c <- tryCatch(wilcox.test(d_list[["UFM1_dependent"]], d_list[["Constitutive"]])$p.value, error=function(e) NA)
    res_p_c <- tryCatch(wilcox.test(d_list[["UFM1_independent"]], d_list[["Constitutive"]])$p.value, error=function(e) NA)
    res_l_p <- tryCatch(wilcox.test(d_list[["UFM1_dependent"]], d_list[["UFM1_independent"]])$p.value, error=function(e) NA)
    
    text(0,0, paste0(
        "Wilcoxon Stats:\n",
        "UFM1_Dep vs Const: p=", safe_fmt(res_l_c), "\n",
        "UFM1_Indep vs Const: p=", safe_fmt(res_p_c), "\n",
        "UFM1_Dep vs UFM1_Indep:  p=", safe_fmt(res_l_p)
    ), cex=1.2)
    
    # 4. Empty or Cumulative?
    # CDF
    plot(NA, xlim=x_lims, ylim=c(0,1), xlab=feat, ylab="CDF", main="Cumulative Distribution")
    for(g in names(d_list)){
       if(length(d_list[[g]]) > 0) {
           ec <- ecdf(d_list[[g]])
           plot(ec, add=TRUE, col=cols[g], verticals=TRUE, pch=NA)
       }
    }
}

dev.off()
