# ============================================================
# 03_annotation_SingleR.R  |  Hierarchical SingleR annotation
# - Major labels via HPCA + Blueprint
# - Epithelial subtypes via breast reference
# - Builds FinalLabel and Family
# Outputs: tables + updated object (optional save)
# ============================================================

source("R/00_utils.R")
suppressPackageStartupMessages({
  if (!requireNamespace("celldex", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("celldex")
  }
  library(celldex)
  library(SingleR)
  library(SummarizedExperiment)
})

srobj <- readRDS(file.path(RES_DIR, "Integrated_TNBC_NormalEpi.rds"))

# Ensure RNA data layer available
srobj <- JoinLayers(srobj, assay = "RNA")
DefaultAssay(srobj) <- "RNA"
test_log <- safe_get_layer(srobj, assay="RNA", layer="data")

# Major-type using celldex references
hpca <- celldex::HumanPrimaryCellAtlasData()
blue <- celldex::BlueprintEncodeData()

inter_hpca <- intersect(rownames(test_log), rownames(hpca))
pred_hpca  <- SingleR(test = test_log[inter_hpca,], ref = hpca[inter_hpca,],
                      labels = hpca$label.main, assay.type.test = "logcounts", assay.type.ref = "logcounts")
inter_blue <- intersect(rownames(test_log), rownames(blue))
pred_blue  <- SingleR(test = test_log[inter_blue,], ref = blue[inter_blue,],
                      labels = blue$label.main, assay.type.test = "logcounts", assay.type.ref = "logcounts")

cells <- colnames(srobj)
hpca_lab <- setNames(pred_hpca$labels, rownames(pred_hpca))
blue_lab <- setNames(pred_blue$labels,  rownames(pred_blue))
lab_major <- data.frame(hpca = hpca_lab[cells], blue = blue_lab[cells], row.names = cells, check.names = FALSE)
lab_major$consensus <- apply(lab_major[, c("hpca","blue")], 1, function(x){
  tab <- table(x); names(sort(tab, decreasing = TRUE))[1]
})
srobj <- AddMetaData(srobj, metadata = data.frame(MajorType = lab_major$consensus, row.names = cells))

# Epithelial subtyping using user breast reference (if provided)
if (file.exists(BREAST_REF_RDS)) {
  myref <- readRDS(BREAST_REF_RDS)
  labcol_my <- if ("label.fine" %in% colnames(colData(myref))) "label.fine"
  else if ("label.main" %in% colnames(colData(myref))) "label.main"
  else stop("No usable label column in breast reference.")

  is_epi    <- grepl("epithel|keratin", srobj$MajorType, ignore.case = TRUE)
  cells_epi <- colnames(srobj)[is_epi]

  srobj$EpithelialSubtype <- NA_character_
  if (length(cells_epi) > 0) {
    epi_obj <- subset(srobj, cells = cells_epi) %>%
      NormalizeData() %>% FindVariableFeatures(nfeatures = 3000, selection.method = "vst")
    hvgs_epi <- VariableFeatures(epi_obj)

    genes_use <- Reduce(intersect, list(hvgs_epi, rownames(test_log), rownames(myref)))

    Idents(srobj) <- srobj$seurat_clusters
    cl_vec <- as.character(Idents(srobj)[cells_epi]); names(cl_vec) <- cells_epi

    pred_epi_clu <- SingleR(
      test     = test_log[genes_use, cells_epi, drop = FALSE],
      ref      = myref[genes_use, ],
      labels   = colData(myref)[, labcol_my],
      clusters = cl_vec,
      assay.type.test = "logcounts",
      assay.type.ref  = "logcounts",
      de.n = 100
    )

    cl2lab_epi <- setNames(pred_epi_clu$labels, rownames(pred_epi_clu))
    srobj$EpithelialSubtype[cells_epi] <- cl2lab_epi[as.character(Idents(srobj)[cells_epi])]
    srobj$FinalLabel      <- ifelse(!is.na(srobj$EpithelialSubtype), srobj$EpithelialSubtype, srobj$MajorType)
    srobj$FinalLabel_main <- ifelse(!is.na(srobj$EpithelialSubtype),
                                    sub("-.*", "", srobj$EpithelialSubtype),
                                    srobj$MajorType)
  } else {
    srobj$FinalLabel      <- srobj$MajorType
    srobj$FinalLabel_main <- srobj$MajorType
  }
} else {
  warning("Breast reference not found; using MajorType only.")
  srobj$FinalLabel      <- srobj$MajorType
  srobj$FinalLabel_main <- srobj$MajorType
}

# Map FinalLabel to families for downstream
family_map <- function(x){
  x <- as.character(x)
  dplyr::case_when(
    grepl("Luminal", x, TRUE)                       ~ "Epithelial-Luminal",
    grepl("\\bBasal\\b", x, TRUE)                   ~ "Epithelial-Basal",
    grepl("\\bMyo\\b", x, TRUE)                     ~ "Epithelial-Myo",
    grepl("Endothel|Pericyte", x, TRUE)             ~ "Endothelial/Pericyte",
    grepl("Fibro|Stromal|MSC|Chondro|Myocyte|Smooth", x, TRUE) ~ "Fibro/MSC",
    grepl("Macroph|Mono|DC|Neutro|GMP", x, TRUE)    ~ "Myeloid",
    grepl("CD4|CD8|NK|T[ _-]?cell", x, TRUE)        ~ "T/NK",
    grepl("\\bB[ _-]?cell\\b|Plasma", x, TRUE)      ~ "B/Plasma",
    grepl("Erythro|Megakaryo", x, TRUE)             ~ "Ery/MK",
    TRUE                                            ~ "Other"
  )
}
srobj$Family <- factor(family_map(srobj$FinalLabel),
                       levels = c("Epithelial-Luminal","Epithelial-Basal","Epithelial-Myo",
                                  "Endothelial/Pericyte","Fibro/MSC","Myeloid","T/NK",
                                  "B/Plasma","Ery/MK","Other"))

# Save
saveRDS(srobj, file.path(RES_DIR, "Annotated_TNBC_NormalEpi.rds"), compress = TRUE)

# Plots
DefaultAssay(srobj) <- "integrated"
fam_cols <- c("pink","#ff7f0e","#2ca","#9467bd","#8c5","#e377c2","#17becf","#bcbd22","tomato","thistle3")
p_fam <- DimPlot(srobj, group.by = "Family", cols = fam_cols,
                 label = TRUE, repel = TRUE, label.size = 3) + NoLegend() + theme_clean()
ggsave(file.path(RES_DIR,"dimplots/umap_by_Family.jpg"), p_fam, width=9, height=7, dpi=300)

p_split <- DimPlot(srobj, group.by="Family", split.by="Class", cols=fam_cols, ncol=2, label = TRUE) +
  NoLegend() + theme_clean()
ggsave(file.path(RES_DIR,"dimplots/umap_Family_split_Class.jpg"), p_split, width=12, height=6, dpi=300)
