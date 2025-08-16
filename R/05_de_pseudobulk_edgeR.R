# ============================================================
# 05_pseudobulk_edgeR.R  |  edgeR DE (TNBC vs NormalEpi) by Family
# Outputs per family under results/tables/DE_by_family/*
# Figures: volcano plots under results/figures/volcano/*
# ============================================================

source("R/00_utils.R")
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("edgeR", ask = FALSE, update = FALSE)
  library(edgeR)
  library(ggplot2)
  library(Matrix)
})

srobj <- readRDS(file.path(RES_DIR, "Annotated_TNBC_NormalEpi.rds"))

# Count family sizes per class and keep those with enough cells
ct_counts <- as.data.frame.matrix(table(srobj$Family, srobj$Class))
ct_counts$Family <- rownames(ct_counts)
min_cells <- 200
eligible <- ct_counts %>% filter(NormalEpi >= min_cells, TNBC >= min_cells) %>% pull(Family)
dir.create(file.path(RES_DIR,"DE_by_family"), recursive = TRUE, showWarnings = FALSE)
fwrite(ct_counts[,c("Family","NormalEpi","TNBC")], file.path(RES_DIR,"celltype_counts_by_class_Family.csv"))

run_pseudobulk_edgeR <- function(obj, family_name, outdir){
  sub <- subset(obj, subset = Family == family_name & !is.na(Class) & Class %in% c("NormalEpi","TNBC"))
  if (ncol(sub) < 50) return(invisible(NULL))

  counts <- tryCatch(GetAssayData(sub, assay="RNA", slot="counts"), error = function(e) NULL)
  if (is.null(counts) || nrow(counts) == 0) counts <- SeuratObject::LayerData(sub[["RNA"]], layer = "counts")

  samp <- sub$orig.ident
  grp  <- sub$Class
  key  <- paste(samp, grp, sep="|")
  split_cols <- split(seq_len(ncol(sub)), key)
  pb_mat <- do.call(cbind, lapply(split_cols, function(idx) Matrix::rowSums(counts[, idx, drop = FALSE])))
  colnames(pb_mat) <- names(split_cols)

  md <- do.call(rbind, strsplit(colnames(pb_mat), "\\|"))
  md <- data.frame(sample = md[,1], group = md[,2], row.names = colnames(pb_mat), check.names = FALSE)
  if (length(unique(md$group)) < 2) return(invisible(NULL))

  y <- DGEList(pb_mat, group = md$group)
  keep <- filterByExpr(y, group = md$group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ md$group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)

  tt <- topTags(qlf, n = Inf)$table
  tt$gene <- rownames(tt)

  safe_family <- gsub("[^A-Za-z0-9_-]", "_", family_name)
  odir <- file.path(outdir, safe_family)
  dir.create(odir, recursive = TRUE, showWarnings = FALSE)

  fwrite(tt, file.path(odir, sprintf("DE_%s_TNBC_vs_NormalEpi.csv", safe_family)))

  tt$FDRsig <- tt$FDR < 0.05
  g <- ggplot(tt, aes(logFC, -log10(FDR))) +
    geom_point(aes(alpha = FDRsig)) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = paste("TNBC vs NormalEpi –", family_name),
         x = "log2FC (TNBC / NormalEpi)", y = "-log10(FDR)") +
    theme_clean(12)
  ggsave(file.path(odir, "volcano.jpg"), g, width = 6, height = 5, dpi = 300)

  fwrite(head(tt[order(tt$FDR), ], 50), file.path(odir, "top50_by_FDR.csv"))
  fwrite(head(tt[order(-tt$logFC), ], 50), file.path(odir, "top50_up_in_TNBC.csv"))
  fwrite(head(tt[order( tt$logFC), ], 50), file.path(odir, "top50_down_in_TNBC.csv"))

  invisible(tt)
}

outdir <- file.path(RES_DIR, "DE_by_family")
res_list <- lapply(eligible, function(f) run_pseudobulk_edgeR(srobj, f, outdir))
