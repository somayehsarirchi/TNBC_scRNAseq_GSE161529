# Triple-Negative Breast Cancer (TNBC) scRNA-seq Analysis — GSE161529

This repository contains a fully reproducible **single-cell RNA-seq workflow** for Triple-Negative Breast Cancer (TNBC) versus normal epithelial tissue, using the publicly available **GSE161529** dataset.

The pipeline integrates:

- quality control and filtering  
- data integration and batch correction  
- automated cell type annotation  
- differential expression and functional enrichment  
- pseudobulk analysis  
- copy number variation (CNV) analysis  

into a clear, modular **R-based framework**.

---

## 1. Biological Context

Triple-Negative Breast Cancer (TNBC) is one of the most aggressive breast cancer subtypes, characterized by:

- lack of ER, PR, and HER2 expression  
- high heterogeneity and poor prognosis  
- frequent genomic instability and CNV events  

By comparing **TNBC vs. normal epithelial cells** at single-cell resolution, this project aims to:

- identify **cell type–specific transcriptional changes**  
- map **enriched pathways** linked to tumour aggressiveness  
- detect **CNV patterns** that may drive tumour progression  

---

## 2. Analysis Goals

The main goals of this workflow are to:

1. Build a **reproducible scRNA-seq pipeline** for TNBC datasets.  
2. Characterize **transcriptional programs and pathways** associated with malignant epithelial subpopulations.  
3. Explore **CNV-driven changes** using inferred copy number profiles at the single-cell level.  

This repository can be used as a **template** for other cancer single-cell projects (e.g., breast, lung, colorectal).

---

## 3. Repository Structure

```text
TNBC_scRNAseq_GSE161529/
├── TNBC_scRNAseq_report.Rmd     # Main R Markdown report (analysis + figures)
├── R/                           # R scripts for each analysis step
│   ├── 01_load_qc.R
│   ├── 02_integration.R
│   ├── 03_annotation_SingleR.R
│   ├── 04_markers.R
│   ├── 05_de_pseudobulk_edgeR.R
│   ├── 06_enrichment_analysis.R
│   └── 07_cnv_analysis.R
├── results/                     # Auto-generated outputs (figures, tables, etc.)
├── LICENSE
└── README.md
Note: The results/ directory is intended for auto-generated outputs and can be safely excluded from version control using .gitignore.
All figures and tables shown in the report are reproducible by running the scripts and R Markdown file.

4. Pipeline Summary

The workflow follows a modular structure to ensure clarity and reproducibility.

Order of execution:

01_load_qc.R

Load raw 10X data, perform quality control, and filter cells/genes.

02_integration.R

Integrate TNBC and normal epithelial datasets using Seurat (RPCA or similar methods).

Compute dimensionality reduction (PCA, UMAP) and clustering.

03_annotation_SingleR.R

Automated cell type annotation with SingleR using a breast epithelial reference.

04_markers.R

Identify cluster-specific markers for TNBC vs. normal epithelial cells.

05_de_pseudobulk_edgeR.R

Pseudobulk differential expression analysis using edgeR.

06_enrichment_analysis.R

Functional enrichment (GO / KEGG) using enrichR or equivalent tools.

07_cnv_analysis.R

CNV detection and visualization with CopyKAT, highlighting malignant subpopulations.

The R Markdown report TNBC_scRNAseq_report.Rmd sources these scripts and compiles the figures/tables into a single, publication-ready document.

5. Example Outputs

The report includes, among others:

UMAP visualizations of TNBC vs. normal epithelial cells

Cluster-level marker expression and heatmaps

Pathway enrichment barplots (e.g., immune response, antigen presentation)

CNV maps distinguishing malignant vs. non-malignant epithelial cells

These outputs are automatically regenerated when the pipeline is re-run.

6. Software Requirements

This workflow was developed and tested with:

R ≥ 4.3.x

Key R packages:

Seurat (v5)

SingleR

celldex

edgeR

data.table, dplyr, ggplot2

enrichR

CopyKAT

Please install the required packages before running the analysis.

7. How to Run

Clone the repository:

git clone https://github.com/somayehsarrichi/TNBC_scRNAseq_GSE161529.git
cd TNBC_scRNAseq_GSE161529

Option A – Run step-by-step scripts

From the project root:

Rscript R/01_load_qc.R
Rscript R/02_integration.R
Rscript R/03_annotation_SingleR.R
Rscript R/04_markers.R
Rscript R/05_de_pseudobulk_edgeR.R
Rscript R/06_enrichment_analysis.R
Rscript R/07_cnv_analysis.R


All outputs will be written to the results/ directory (figures, tables, CNV calls, etc.).

Option B – Knit the full report

Open TNBC_scRNAseq_report.Rmd in RStudio (or run via command line) and knit to HTML:

rmarkdown::render("TNBC_scRNAseq_report.Rmd")


This will:

source the scripts in R/

run the full analysis

generate a unified report with all figures and tables

8. Dataset

The analysis is based on the public dataset:

GSE161529 – Triple-Negative Breast Cancer (TNBC) scRNA-seq

Raw data can be obtained from GEO and pre-processed into 10X-like matrices as described in the dataset’s original documentation.

9. Citation

If you use this repository or scripts in your work, please cite:

Sarirchi, S. (2025). TNBC scRNA-seq Analysis – GSE161529 (v1.0.0). Zenodo.
DOI: 10.5281/zenodo.17127154

10. Contact

For questions, collaborations, or consulting inquiries:

Email: somayeh.sarirchi@gmail.com

LinkedIn: https://www.linkedin.com/in/somayeh-sarirchi-9b25b9171

GitHub: https://github.com/somayehsarrichi
