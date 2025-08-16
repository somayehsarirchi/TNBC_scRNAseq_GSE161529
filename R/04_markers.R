# 04_markers.R — Top markers per cluster + DotPlots

source("R/00_utils.R")
suppressPackageStartupMessages({
  library(ggplot2)
})

srobj <- readRDS(file.path(RES_DIR, "Annotated_TNBC_NormalEpi.rds"))

DefaultAssay(srobj) <- "RNA"
markers <- FindAllMarkers(srobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fwrite(markers, file.path(RES_DIR, "AllMarkers_all.csv"))

top20_tbl <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  select(cluster, gene, avg_log2FC, p_val_adj)
fwrite(top20_tbl, file.path(RES_DIR, "Top20-Markers_byCluster.csv"))

# Canonical gene panels
genes <- c(
  "EPCAM","KRT8","KRT18","KRT19","MUC1",
  "KRT5","KRT14","KRT17","ACTA2","TAGLN","MYL9","CALD1",
  "VWF","PECAM1","KDR",
  "CD3D","CD3E","NKG7","KLRD1",
  "CD79A","CD79B","CD74","MZB1","XBP1",
  "LYZ","S100A8","S100A9","LST1"
)
genes <- intersect(genes, rownames(srobj))
stopifnot(length(genes) > 0)

p <- DotPlot(srobj, features = genes, group.by = "Family", assay = "RNA", dot.scale = 4) +
  RotatedAxis() + theme_clean()
ggsave(file.path(RES_DIR, "dimplots/dotplot_canonical_family.jpg"), p, width = 10, height = 5, dpi = 300)

dp <- DotPlot(srobj, features = genes, group.by = "Family", split.by = "Class", assay = "RNA", dot.scale = 5)
df <- dp$data
df$Class  <- sub(".*_", "", df$id)
df$Family <- sub("_.*", "", df$id)
group_levels <- unique(df$Family)
split_levels <- as.vector(sapply(group_levels, function(fam) c(paste0(fam,"_NormalEpi"), paste0(fam,"_TNBC"))))
df$id <- factor(df$id, levels = rev(split_levels))

p_split <- ggplot(df, aes(x = features.plot, y = id, size = pct.exp, colour = avg.exp.scaled)) +
  geom_point() +
  scale_colour_gradient(low = "grey90", high = "royalblue3", name = "Average Expression") +
  scale_size(name = "Percent Expressed", range = c(0, 6), breaks = c(0,25,50,75,100)) +
  theme_clean(12) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Features", y = "Split Identity")
ggsave(file.path(RES_DIR, "dimplots/dotplot_canonical_family_split_ordered.jpg"), p_split, width = 12, height = 6, dpi = 300)
