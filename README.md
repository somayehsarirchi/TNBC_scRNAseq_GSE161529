# TNBC scRNA-seq (GSE161529)

End-to-end pipeline: QC → integration (RPCA) → clustering/annotation → DE (pseudobulk edgeR) → enrichment (enrichR) → CNV (CopyKAT) → (CellChat/trajectory optional).

## Quick start
```r
# Reproducible env
install.packages("renv")
renv::init()          # creates renv.lock
# run scripts in R/ in order
