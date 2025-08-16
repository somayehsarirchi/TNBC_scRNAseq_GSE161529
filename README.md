# TNBC scRNA-seq (GSE161529) — End-to-end Reproducible Workflow

**Goal:** An end-to-end single-cell RNA-seq analysis of TNBC vs Normal epithelium:
QC → Integration (RPCA) → Annotation (SingleR) → Markers → Pseudobulk DE (edgeR) → Enrichment (enrichR) → Pathway Prioritization (A/B) → (optional) CNV (CopyKAT).

## Data
- GEO: GSE161529 (raw 10x matrices; not included).  
- Reference for epithelial subtyping: Breast reference RDS (path noted in scripts).

## Environment
- R version: `x.y.z`
- Main packages: Seurat (v5), SingleR, celldex, edgeR, enrichR, data.table, dplyr, ggplot2.
- Optional: `renv` lockfile for exact versions.

## How to run (reproducible order)
1. `Rscript scripts/01_qc_load.R`
2. `Rscript scripts/02_integration.R`
3. `Rscript scripts/03_annotation_SingleR.R`
4. `Rscript scripts/04_markers.R`
5. `Rscript scripts/05_pseudobulk_edgeR.R`
6. `Rscript scripts/06_enrichment_prioritization.R`
7. `Rscript scripts/07_cnv_copykat.R` *(optional; long)*

Outputs go to `results/figures/` and `results/tables/`.

## Results preview
- UMAPs, violin plots: `results/figures/`
- All markers and Top20 by cluster: `results/tables/`
- DE by family + enrichment summaries and prioritized pathways.

## Citation
Please cite:
- **Dataset:** GSE161529 (link)
- **Tools:** Seurat, SingleR, celldex, edgeR, enrichR, CopyKAT.

## License
MIT
