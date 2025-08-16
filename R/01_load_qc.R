# ============================================================
# 01_qc_load.R  |  Load raw 10x, basic QC, save QC’d objects
# Project: TNBC scRNA-seq (GSE161529)
# Author: Somayeh Sarirchi  |  License: MIT
# ------------------------------------------------------------
# What this does:
# - Reads 10x matrices for TNBC and Normal-Epi
# - Basic QC filters (nFeature_RNA, MT%)
# - Normalization + HVG
# - Saves TNBC_QC_HVG.rds and NormalEpi_QC_HVG.rds
# Reproducibility: set.seed(), sessionInfo() at end.
# ============================================================


source("R/00_utils.R")
suppressPackageStartupMessages({
  library(Matrix)
})

# ---- Edit these to your local paths to 10x dirs ----
TNBC_SAMPLES <- c(
  "TNBC/GSM4909281/10X/","TNBC/GSM4909282/10X/","TNBC/GSM4909283/10X/",
  "TNBC/GSM4909284/10X/","TNBC/GSM4909285/10X/","TNBC/GSM4909286/10X/",
  "TNBC/GSM4909287/10X/","TNBC/GSM4909288/10X/"
)
NEPI_SAMPLES <- c(
  "Normal-Epi/GSM4909255/10X/","Normal-Epi/GSM4909256/10X/","Normal-Epi/GSM4909258/10X/",
  "Normal-Epi/GSM4909259/10X/","Normal-Epi/GSM4909260/10X/","Normal-Epi/GSM4909264/10X/",
  "Normal-Epi/GSM4909267/10X/","Normal-Epi/GSM4909273/10X/","Normal-Epi/GSM4909275/10X/"
)

# Helper to read one 10x folder
read_one_10x <- function(path, project){
  m <- Read10X(file.path(BASE_DIR, path))
  CreateSeuratObject(m, min.cells = 3, min.features = 200, project = project)
}

# ---- TNBC ----
t_list <- mapply(function(path, prj){
  read_one_10x(path, prj)
}, TNBC_SAMPLES, gsub(".*/(GSM\\d+).*", "\\1", TNBC_SAMPLES), SIMPLIFY = FALSE)

TNBC <- Reduce(function(a,b) merge(a,b, add.cell.ids = c(a@project.name, b@project.name)), t_list)
TNBC[["MTpercent"]] <- PercentageFeatureSet(TNBC, pattern = "^MT-")
TNBC <- subset(TNBC, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & MTpercent < 20)
TNBC <- NormalizeData(TNBC) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

# ---- Normal Epithelium ----
n_list <- mapply(function(path, prj){
  read_one_10x(path, prj)
}, NEPI_SAMPLES, gsub(".*/(GSM\\d+).*", "\\1", NEPI_SAMPLES), SIMPLIFY = FALSE)

NormalEpi <- Reduce(function(a,b) merge(a,b, add.cell.ids = c(a@project.name, b@project.name)), n_list)
NormalEpi[["MTpercent"]] <- PercentageFeatureSet(NormalEpi, pattern = "^MT-")
NormalEpi <- subset(NormalEpi, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & MTpercent < 20)
NormalEpi <- NormalizeData(NormalEpi) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

# Save QC’d objects
saveRDS(TNBC,      file.path(RES_DIR, "TNBC_QC_HVG.rds"),      compress = FALSE)
saveRDS(NormalEpi, file.path(RES_DIR, "NormalEpi_QC_HVG.rds"), compress = FALSE)

