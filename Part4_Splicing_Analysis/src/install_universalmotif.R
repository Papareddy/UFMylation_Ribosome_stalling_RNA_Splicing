
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="http://cran.us.r-project.org")

# Try installing universalmotif again
BiocManager::install("universalmotif", update=FALSE, ask=FALSE)
