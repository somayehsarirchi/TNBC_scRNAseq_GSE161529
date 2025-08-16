# 00_utils.R — shared helpers & global config

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(Seurat)
  library(SeuratObject)
})

# ---- Paths ----
BASE_DIR <- normalizePath(".", winslash = "/")  # project root
RES_DIR  <- file.path(BASE_DIR, "results")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RES_DIR, "dimplots"), recursive = TRUE, showWarnings = FALSE)

# Optional: where your breast SingleR reference lives (edit if you have it)
BREAST_REF_RDS <- file.path(BASE_DIR, "data", "breast_SingleR_ref.rds")

# ---- Safe layer accessor (Seurat v5+ compatible) ----
safe_get_layer <- function(obj, assay = "RNA", layer = c("counts","data","scale.data")){
  layer <- match.arg(layer)
  if (!assay %in% names(obj@assays)) stop("Assay not found: ", assay)
  a <- obj[[assay]]
  if ("LayerData" %in% ls(getNamespace("SeuratObject"))) {
    SeuratObject::LayerData(a, layer = layer)
  } else {
    slot <- switch(layer, counts = "counts", data = "data", `scale.data` = "scale.data")
    GetAssayData(obj, assay = assay, slot = slot)
  }
}

# ---- Small ggplot theme ----
theme_clean <- function(base_size = 12){
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = element_blank())
}
