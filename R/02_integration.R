# 02_integration.R — Integrate (RPCA) + UMAP + clusters

source("R/00_utils.R")

TNBC      <- readRDS(file.path(RES_DIR, "TNBC_QC_HVG.rds"))
NormalEpi <- readRDS(file.path(RES_DIR, "NormalEpi_QC_HVG.rds"))

# make layers explicit (v5 safety)
TNBC      <- JoinLayers(TNBC, assay = "RNA")
NormalEpi <- JoinLayers(NormalEpi, assay = "RNA")
DefaultAssay(TNBC) <- "RNA"; DefaultAssay(NormalEpi) <- "RNA"

TNBC      <- NormalizeData(TNBC) %>% FindVariableFeatures(nfeatures = 2000)
NormalEpi <- NormalizeData(NormalEpi) %>% FindVariableFeatures(nfeatures = 2000)

features <- SelectIntegrationFeatures(object.list = list(NormalEpi, TNBC), nfeatures = 3000)
TNBC      <- ScaleData(TNBC, features = features, verbose = FALSE) %>% RunPCA(features = features, verbose = FALSE)
NormalEpi <- ScaleData(NormalEpi, features = features, verbose = FALSE) %>% RunPCA(features = features, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = list(NormalEpi, TNBC),
                                  anchor.features = features, reduction = "rpca",
                                  dims = 1:30, k.anchor = 5, k.filter = 200)

saveRDS(anchors, file.path(RES_DIR, "Anchors_joined_rpca.rds"), compress = TRUE)

srobj <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(srobj) <- "integrated"

srobj <- ScaleData(srobj, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
srobj <- FindNeighbors(srobj, dims = 1:30) %>% FindClusters(resolution = 0.5)

# Class from cell barcode prefixes (add.cell.ids)
srobj$Class <- ifelse(grepl("^NormalEpi", colnames(srobj)), "NormalEpi", "TNBC")

saveRDS(srobj, file.path(RES_DIR, "Integrated_TNBC_NormalEpi.rds"), compress = TRUE)

# quick plots
p1 <- DimPlot(srobj, group.by = "Class") + theme_clean()
p2 <- DimPlot(srobj, label = TRUE) + NoLegend() + theme_clean()
ggsave(file.path(RES_DIR, "dimplots/umap_class.jpg"),    p1, width=8, height=6, dpi=300)
ggsave(file.path(RES_DIR, "dimplots/umap_clusters.jpg"), p2, width=8, height=6, dpi=300)
