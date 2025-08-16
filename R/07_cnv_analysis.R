# ------------------------------------------------------------
# Script: 07_cnv_analysis.R
# Goal:
#   1) Run CopyKAT on TNBC cells only and attach 'Malignancy' labels
#   2) Compute DEGs: aneuploid vs diploid (within TNBC cells)
#   3) Run enrichR on UP-regulated genes and save plots/tables
#
# Outputs (under results/ and results/plots/):
#   - results/copykat_TNBC_only.rds
#   - results/DEGs_TNBC_aneuploid_vs_diploid.csv
#   - results/DEGs_TNBC_aneuploid_vs_diploid-Sig.csv
#   - results/enrich_GO_TNBC_up.csv
#   - results/enrich_KEGG_TNBC_up.csv
#   - results/plots/enrich_GO_up_top10_padj.pdf/png
#   - results/plots/enrich_KEGG_up_top10_padj.pdf/png
#   - results/plots/enrich_GO_up_top10_combined.pdf/png
#   - results/plots/enrich_KEGG_up_top10_combined.pdf/png
#
# Notes:
#   - This script mirrors the logic you used during development, with English comments.
#   - Uses Seurat v5-friendly access (JoinLayers + layer='counts').
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(copykat)
  library(dplyr)
  library(ggplot2)
})

# ------------------ CONFIG ------------------
# Base path = project root (assumes script is in 'R/' folder)
base <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)

obj_path <- file.path(base, "Integrated_TNBC_NormalEpi.rds")
out_dir  <- file.path(base, "results")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)


# TNBC / Normal sample IDs
tnbc_ids <- c("GSM4909281","GSM4909282","GSM4909283","GSM4909284",
              "GSM4909285","GSM4909286","GSM4909287","GSM4909288")
normal_ids <- c("GSM4909255","GSM4909256","GSM4909258","GSM4909259",
                "GSM4909260","GSM4909264","GSM4909267","GSM4909273","GSM4909275")

# ------------------ LOAD OBJECT ------------------
srobj <- readRDS(obj_path)
DefaultAssay(srobj) <- "RNA"
# Ensure v5-safe access to counts
try(srobj <- JoinLayers(srobj, assay = "RNA"), silent = TRUE)
DefaultAssay(srobj) <- "RNA"

# Label classes (NormalEpi vs TNBC)
srobj$Class <- ifelse(srobj$orig.ident %in% normal_ids, "NormalEpi",
                      ifelse(srobj$orig.ident %in% tnbc_ids, "TNBC", NA))
srobj$Class <- factor(srobj$Class, levels = c("NormalEpi","TNBC"))
stopifnot(!any(is.na(srobj$Class)))

# ------------------ 1) COPYKAT (TNBC ONLY) ------------------
tnbc_cells_all <- WhichCells(srobj, expression = Class == "TNBC")
var.genes <- VariableFeatures(srobj)

# counts layer (Seurat v5)
counts_mat <- GetAssayData(srobj, assay = "RNA", layer = "counts")[var.genes, tnbc_cells_all]

ck <- copykat(
  rawmat   = as.matrix(counts_mat),
  id.type  = "S",   # human SYMBOL
  ngene.chr= 5,
  win.size = 25,
  KS.cut   = 0.1,
  n.cores  = 1
)

saveRDS(ck, file.path(out_dir, "copykat_TNBC_only.rds"))
pred <- ck$prediction

# Map CopyKAT predictions back to Seurat cells via barcode
# (assumes cell names like 'TNBC<idx>_<BARCODE>-1')
tnbc_cells <- colnames(srobj)[srobj$Class == "TNBC"]
seurat_map <- data.frame(
  full    = tnbc_cells,
  barcode = gsub("^TNBC\\d+_([ACGT]+)-1$", "\\1", tnbc_cells)
)

pred$barcode <- gsub("^.*?_([ACGT]+)-1$", "\\1", pred$cell.names)
matched <- merge(pred, seurat_map, by = "barcode", all.x = FALSE, all.y = FALSE)

srobj$Malignancy <- NA
srobj$Malignancy[matched$full] <- matched$copykat.pred
srobj$Malignancy <- factor(srobj$Malignancy, levels = c("diploid","aneuploid"))

# Basic check & quick UMAP
print(table(srobj$Malignancy, useNA = "always"))
if ("umap" %in% names(srobj@reductions)) {
  p_umap <- DimPlot(srobj, group.by = "Malignancy", pt.size = 0.3) +
    ggtitle("CopyKAT CNV Prediction — diploid vs aneuploid (TNBC only)")
  ggsave(file.path(plot_dir, "umap_malignancy.pdf"), p_umap, width = 7, height = 6)
  ggsave(file.path(plot_dir, "umap_malignancy.png"), p_umap, width = 7, height = 6, dpi = 300)
}

# ------------------ 2) DEGs: aneuploid vs diploid (TNBC) ------------------
tnbc_cells_valid <- WhichCells(srobj, expression = Class == "TNBC" & !is.na(Malignancy))
subobj <- subset(srobj, cells = tnbc_cells_valid)
Idents(subobj) <- "Malignancy"

deg_cnv <- FindMarkers(
  subobj,
  ident.1 = "aneuploid",
  ident.2 = "diploid",
  logfc.threshold = 0.5,
  min.pct = 0.25
)

deg_cnv$gene <- rownames(deg_cnv)
write.csv(deg_cnv, file.path(out_dir, "DEGs_TNBC_aneuploid_vs_diploid.csv"), row.names = FALSE)

deg_sig <- deg_cnv %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)
write.csv(deg_sig, file.path(out_dir, "DEGs_TNBC_aneuploid_vs_diploid-Sig.csv"), row.names = FALSE)

# ------------------ 3) EnrichR on UP genes ------------------
suppressPackageStartupMessages({
  if (!requireNamespace("enrichR", quietly = TRUE)) install.packages("enrichR")
  library(enrichR)
})

deg_sig <- read.csv(file.path(out_dir, "DEGs_TNBC_aneuploid_vs_diploid-Sig.csv"), check.names = FALSE)

up_genes <- deg_sig %>%
  filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
  pull(gene) %>% unique()

dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human")
enr_up <- if (length(up_genes)) enrichr(up_genes, dbs) else NULL

if (!is.null(enr_up)) {
  write.csv(enr_up[["GO_Biological_Process_2023"]],
            file.path(out_dir, "enrich_GO_TNBC_up.csv"), row.names = FALSE)
  write.csv(enr_up[["KEGG_2021_Human"]],
            file.path(out_dir, "enrich_KEGG_TNBC_up.csv"), row.names = FALSE)

  # --- Plot helpers (Top 10 by adj.P and by Combined Score) ---
  plot_enr <- function(df, title, outfile, by = c("padj","combined")) {
    if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
    by <- match.arg(by)
    if (by == "padj") {
      d <- df %>% arrange(Adjusted.P.value) %>% slice_head(n = 10)
      p <- ggplot(d, aes(x = reorder(Term, -log10(Adjusted.P.value)), y = -log10(Adjusted.P.value))) +
        geom_col() + coord_flip() + theme_minimal() +
        labs(title = title, x = NULL, y = "-log10(adj. P)")
    } else {
      d <- df %>% arrange(desc(Combined.Score)) %>% slice_head(n = 10)
      p <- ggplot(d, aes(x = reorder(Term, Combined.Score), y = Combined.Score)) +
        geom_col() + coord_flip() + theme_minimal() +
        labs(title = title, x = NULL, y = "Combined Score")
    }
    ggsave(file.path(plot_dir, paste0(outfile, ".pdf")), p, width = 7, height = 5)
    ggsave(file.path(plot_dir, paste0(outfile, ".png")), p, width = 7, height = 5, dpi = 300)
    invisible(p)
  }

  plot_enr(enr_up[["GO_Biological_Process_2023"]],
           "Top 10 GO BP (Upregulated in Aneuploid, by adj.P)",
           "enrich_GO_up_top10_padj", by = "padj")
  plot_enr(enr_up[["KEGG_2021_Human"]],
           "Top 10 KEGG (Upregulated in Aneuploid, by adj.P)",
           "enrich_KEGG_up_top10_padj", by = "padj")
  plot_enr(enr_up[["GO_Biological_Process_2023"]],
           "Top 10 GO BP (Upregulated in Aneuploid, by Combined Score)",
           "enrich_GO_up_top10_combined", by = "combined")
  plot_enr(enr_up[["KEGG_2021_Human"]],
           "Top 10 KEGG (Upregulated in Aneuploid, by Combined Score)",
           "enrich_KEGG_up_top10_combined", by = "combined")
} else {
  message("No UP genes passed the thresholds; enrichment skipped.")
}

message("DONE: 07_cnv_analysis.R")
