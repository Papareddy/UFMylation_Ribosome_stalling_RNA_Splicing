
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="http://cran.us.r-project.org")

BiocManager::install("MotifDb", ask=FALSE)
BiocManager::install("universalmotif", ask=FALSE)
BiocManager::install("AnnotationHub", ask=FALSE)
