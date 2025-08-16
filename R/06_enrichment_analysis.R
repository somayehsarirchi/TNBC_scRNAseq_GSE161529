# ------------------------------------------------------------
# Script: 06_enrichment_analysis.R
# Goal:
#   1) Build gene sets from pseudobulk DE (step 05 outputs)
#   2) Run enrichR (GO/KEGG) and save full + filtered tables
#   3) Prioritize pathways with two profiles (A neutral, B TNBC-aware)
# Output:
#   results/summary_DE/
#     ├─ enrich_full/               (all enrichment tables)
#     ├─ enrich_sig/                (filtered by padj & combined score)
#     └─ enrich_sig/prioritized/    (A/B prioritized CSV + Top20 plot)
# Notes:
#   - Designed to be idempotent: set cfg$overwrite = TRUE to overwrite.
#   - Safe to re-run after adding new DE families.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(ggplot2)
  if (!requireNamespace("enrichR", quietly = TRUE)) install.packages("enrichR")
  library(enrichR)
})

## ---------------- CONFIG ----------------
cfg <- list(
  base         = "C:/Users/asus/Desktop/Job/Work/GSE161529",
  de_dir       = "results/DE_by_family",
  out_root     = "results/summary_DE",
  dbs          = c("GO_Biological_Process_2023","KEGG_2021_Human"),
  de_logfc     = 1,        # DEG threshold (|logFC| >= 1)
  de_fdr       = 0.05,     # DEG threshold (FDR < 0.05)
  enr_padj     = 0.05,     # adjusted p-value for enrichment filtering
  enr_comb_min = 20,       # min Combined Score; set 0 to disable
  min_genes    = 10,       # skip enrichment when smaller than this
  overwrite    = FALSE     # set TRUE to overwrite existing outputs
)

## --------- paths ---------
cfg$de_dir   <- file.path(cfg$base, cfg$de_dir)
cfg$out_root <- file.path(cfg$base, cfg$out_root)
dir.create(cfg$out_root, recursive = TRUE, showWarnings = FALSE)
dir_full <- file.path(cfg$out_root, "enrich_full")
dir_sig  <- file.path(cfg$out_root, "enrich_sig")
dir_prio <- file.path(dir_sig, "prioritized")
dir.create(dir_full, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_sig,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_prio, recursive = TRUE, showWarnings = FALSE)

## --------------- helpers ---------------
to_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9eE+\\-\\.]", "", gsub(",", ".", trimws(x)))))
scale01 <- function(x){ if (all(!is.finite(x))) return(rep(0, length(x)))
  r <- range(x[is.finite(x)], na.rm=TRUE); if (diff(r)==0) return(rep(0.5,length(x))); (x-r[1])/(r[2]-r[1]) }
parse_genes <- function(s){ s <- gsub("\\s+","",as.character(s)); toupper(unlist(strsplit(s,"[,;]"))) }
write_if_needed <- function(dt, path, overwrite=FALSE){ if (!overwrite && file.exists(path)) return(invisible(FALSE)); fwrite(dt, path); TRUE }
find_col <- function(nms, pats, fallback = NA_character_){
  for (p in pats){ hit <- grep(p, nms, ignore.case=TRUE, value=TRUE); if (length(hit)) return(hit[1]) }
  if (!is.na(fallback) && fallback %in% nms) return(fallback); NA_character_
}
fix_overlap <- function(dt){
  dt <- as.data.table(dt)
  if ("Overlap" %in% names(dt)) {
    ov <- as.character(dt[["Overlap"]])
    m  <- tstrsplit(ov, "/", fixed = TRUE)
    if (length(m) == 2) {
      dt[, Overlap_k := suppressWarnings(as.integer(m[[1]]))]
      dt[, Overlap_N := suppressWarnings(as.integer(m[[2]]))]
    }
    dt[, Overlap := NULL]
  } else {
    if (!("Overlap_k" %in% names(dt)) && "Genes" %in% names(dt)) {
      dt[, Overlap_k := vapply(Genes, function(x) length(parse_genes(x)), integer(1))]
    }
    if (!("Overlap_N" %in% names(dt))) dt[, Overlap_N := NA_integer_]
  }
  dt
}

## ---------------- 1) Load DE and build gene sets ----------------
read_de_folder <- function(fdir){
  f <- list.files(fdir, pattern="^DE_.*_TNBC_vs_NormalEpi\\.csv$", full.names=TRUE)
  if (!length(f)) return(NULL)
  dt <- fread(f[1]); req <- c("gene","logFC","FDR")
  if (!all(req %in% names(dt))) return(NULL)
  dt[, Family := basename(fdir)]
  dt[]
}
fam_dirs <- list.dirs(cfg$de_dir, full.names = TRUE, recursive = FALSE)
all_DE <- rbindlist(lapply(fam_dirs, read_de_folder), use.names = TRUE, fill = TRUE)
stopifnot(nrow(all_DE) > 0)

sig <- all_DE %>% filter(FDR < cfg$de_fdr, abs(logFC) >= cfg$de_logfc)
write_if_needed(sig, file.path(cfg$out_root, "DE_sig_all_families.csv"), cfg$overwrite)

up_list   <- split(sig %>% filter(logFC > 0) %>% pull(gene),   sig %>% filter(logFC > 0) %>% pull(Family))
down_list <- split(sig %>% filter(logFC < 0) %>% pull(gene),   sig %>% filter(logFC < 0) %>% pull(Family))

epi_fams <- sort(unique(grep("^Epithelial", names(up_list), value = TRUE)))
common_up_epi   <- if (length(epi_fams) >= 2) Reduce(intersect, up_list[epi_fams]) else character(0)
common_down_epi <- if (length(epi_fams) >= 2) Reduce(intersect, down_list[epi_fams]) else character(0)

## ---------------- 2) Enrichment (full + filtered) ----------------
filter_and_write_enrich <- function(dt, prefix, cfg){
  dt <- fix_overlap(dt)
  fwrite(dt, file.path(dir_full, paste0(prefix, ".csv")))  # full table
  
  padj_col <- if ("Adjusted.P.value" %in% names(dt)) "Adjusted.P.value" else find_col(names(dt), c("adjust.*p.?val"))
  if (is.na(padj_col)) { warning("No Adjusted p-value in ", prefix); return(invisible(NULL)) }
  
  dt[, (padj_col) := to_num(get(padj_col))]
  if ("Combined.Score" %in% names(dt)) dt[, Combined.Score := to_num(Combined.Score)]
  
  if ("Combined.Score" %in% names(dt) && cfg$enr_comb_min > 0) {
    dts <- dt[get(padj_col) < cfg$enr_padj & Combined.Score >= cfg$enr_comb_min]
    if (nrow(dts)) data.table::setorderv(dts, c(padj_col, "Combined.Score"), c(1,-1))
  } else {
    dts <- dt[get(padj_col) < cfg$enr_padj]
    if (nrow(dts)) data.table::setorderv(dts, padj_col, 1)
  }
  if (nrow(dts)) fwrite(dts, file.path(dir_sig, paste0(prefix, "_sig.csv")))
  invisible(NULL)
}

run_enrich_set <- function(genes, set_name, cfg){
  genes <- unique(na.omit(genes))
  if (length(genes) < cfg$min_genes) return(invisible(NULL))
  enr <- enrichr(genes, cfg$dbs)
  for (nm in names(enr)) {
    prefix <- paste0("enrich_", set_name, "_", nm)
    sig_path  <- file.path(dir_sig,  paste0(prefix, "_sig.csv"))
    full_path <- file.path(dir_full, paste0(prefix, ".csv"))
    if (!cfg$overwrite && file.exists(sig_path) && file.exists(full_path)) next
    filter_and_write_enrich(enr[[nm]], prefix, cfg)
  }
}

# epithelial common sets
run_enrich_set(common_up_epi,   "common_UP_epithelial",   cfg)
run_enrich_set(common_down_epi, "common_DOWN_epithelial", cfg)

# per-family (UP only is typical; add DOWN if you like)
for (fam in names(up_list)) {
  set_name <- paste0("UP_", gsub("[^A-Za-z0-9_-]", "_", fam))
  run_enrich_set(up_list[[fam]], set_name, cfg)
}

## ---------------- 3) Prioritization (A vs B) ----------------
toupper_all <- function(x) unique(toupper(x))
TNBC_EPITHELIAL  <- toupper_all(c("EPCAM","KRT8","KRT18","KRT19","MUC1"))
TNBC_BASAL_MYO   <- toupper_all(c("KRT5","KRT14","KRT17","ACTA2","TAGLN","MYL9","EGFR"))
EMT_METASTASIS   <- toupper_all(c("VIM","SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","TWIST2","FN1","ITGA5","ITGB1","CDH1","CDH2","MMP2","MMP9"))
DDR_REPAIR       <- toupper_all(c("BRCA1","BRCA2","RAD51","CHEK1","CHEK2","ATR","ATM","PARP1"))
IMMUNE_CHECKPT   <- toupper_all(c("CD274","PDCD1","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2","IDO1","HLA-DRA","HLA-DRB1","HLA-DPB1"))
PROLIFERATION    <- toupper_all(c("MKI67","TOP2A","PCNA","AURKA","AURKB","PLK1","CCNB1","CCNB2","CDK1"))
ANGIO_VASC       <- toupper_all(c("VEGFA","KDR","FLT1","ANGPT2","TEK","HIF1A"))
EGFR_PI3K_MAPK   <- toupper_all(c("EGFR","ERBB2","PIK3CA","PIK3CB","AKT1","AKT2","PTEN","KRAS","NRAS","BRAF","MAPK1","MAPK3"))

panel_list <- list(
  Epithelial   = TNBC_EPITHELIAL,  Basal_Myo    = TNBC_BASAL_MYO,
  EMT_Meta     = EMT_METASTASIS,   DDR          = DDR_REPAIR,
  Immune_CP    = IMMUNE_CHECKPT,   Proliferation= PROLIFERATION,
  Angiogenesis = ANGIO_VASC,       EGFR_PI3K    = EGFR_PI3K_MAPK
)

kw <- c("cancer","carcinoma","breast","pi3k","akt","mapk","erk","p53",
        "apoptosis","cell cycle","dna repair","ecm","integrin","tgf-?b","jak-?stat")

find_col2 <- function(nms, patterns, fallback = NA_character_){
  for (p in patterns) {
    hit <- grep(p, nms, ignore.case=TRUE, value=TRUE)
    if (length(hit)) return(hit[1])
  }
  if (!is.na(fallback) && fallback %in% nms) return(fallback)
  NA_character_
}
score_core <- function(dt){
  nms <- names(dt)
  padj_col <- find_col2(nms, c("^Adjusted.P.value$","adjust.*p.?val"))
  p_col    <- find_col2(nms, c("^P\\.value$","\\bp\\.?value\\b"))
  comb_col <- find_col2(nms, c("^Combined.Score$","combined\\s*score"))
  genes_col<- find_col2(nms, c("^Genes$","gene.*list"), "Genes")
  k_col    <- find_col2(nms, c("^Overlap_k$","overlap_k"))
  N_col    <- find_col2(nms, c("^Overlap_N$","overlap_n"))
  if (is.na(padj_col)) stop("Adjusted p-value column not found.")
  
  dt[, (padj_col) := to_num(get(padj_col))]
  if (!is.na(p_col))    dt[, (p_col)    := to_num(get(p_col))]
  if (!is.na(comb_col)) dt[, (comb_col) := to_num(get(comb_col))]
  
  if (is.na(k_col)) { dt[, Overlap_k := vapply(get(genes_col), function(x) length(parse_genes(x)), integer(1))]; k_col <- "Overlap_k" }
  if (is.na(N_col)) { dt[, Overlap_N := NA_integer_]; N_col <- "Overlap_N" }
  
  neglog_padj <- -log10(dt[[padj_col]])
  frac_ovlp   <- if (all(is.na(dt[[N_col]]))) rep(NA_real_, nrow(dt)) else dt[[k_col]] / pmax(dt[[N_col]], 1)
  spec_comp   <- if (all(is.na(dt[[N_col]]))) rep(0.5, nrow(dt)) else 1 - scale01(log10(pmax(dt[[N_col]], 2)))
  
  # panel hits
  panel_hits <- matrix(0L, nrow=nrow(dt), ncol=length(panel_list)); colnames(panel_hits) <- names(panel_list)
  for (j in seq_along(panel_list)) {
    gs <- panel_list[[j]]
    panel_hits[, j] <- vapply(dt[[genes_col]], function(x) length(intersect(parse_genes(x), gs)), integer(1))
  }
  panel_hits_dt <- as.data.table(panel_hits)
  
  data.table(
    idx = seq_len(nrow(dt)),
    Term = if ("Term" %in% names(dt)) dt$Term else dt[[1]],
    pval = if (!is.na(p_col)) dt[[p_col]] else NA_real_,
    padj = dt[[padj_col]],
    combined = if (!is.na(comb_col)) dt[[comb_col]] else NA_real_,
    genes = dt[[genes_col]],
    Overlap_k = dt[[k_col]], Overlap_N = dt[[N_col]],
    neglog10_padj = neglog_padj, frac_overlap = frac_ovlp, spec_comp = spec_comp
  )[, (names(panel_hits_dt)) := panel_hits_dt]
}

apply_profile <- function(core_dt, profile=c("A","B")){
  profile <- match.arg(profile)
  w <- list(w_padj=0.60, w_comb=0.30, w_ovlp=0.10, w_tnbc=0.00, w_spec=0.00,
            panel_weight=setNames(rep(0,length(panel_list)), names(panel_list)), kw_boost=0)
  if (profile=="B"){
    w <- list(w_padj=0.55, w_comb=0.25, w_ovlp=0.10, w_tnbc=0.08, w_spec=0.02,
              panel_weight=setNames(c(1,1,1,0.9,0.8,0.7,0.7,0.9), names(panel_list))*0.5,
              kw_boost=0.05)
  }
  if (all(!is.finite(core_dt$combined))){
    s <- w$w_padj + w$w_ovlp + w$w_tnbc + w$w_spec
    w$w_padj <- w$w_padj/s; w$w_comb <- 0; w$w_ovlp <- w$w_ovlp/s; w$w_tnbc <- w$w_tnbc/s; w$w_spec <- w$w_spec/s
  }
  tnbc_hits <- rep(0, nrow(core_dt)); for (nm in names(panel_list))
    if (nm %in% names(core_dt)) tnbc_hits <- tnbc_hits + w$panel_weight[[nm]] * core_dt[[nm]]
  tnbc_norm <- scale01(tnbc_hits)
  kw_bonus  <- w$kw_boost * as.numeric(grepl(paste(kw, collapse="|"), core_dt$Term, ignore.case=TRUE))
  
  score <- w$w_padj*scale01(core_dt$neglog10_padj) +
    w$w_comb*scale01(core_dt$combined) +
    w$w_ovlp*scale01(core_dt$frac_overlap) +
    w$w_tnbc*tnbc_norm +
    w$w_spec*scale01(core_dt$spec_comp) + kw_bonus
  
  out <- copy(core_dt)[, c("priority","rank","tnbc_hits_weighted") := .(score, rank(-score,"min"), tnbc_hits)]
  out[]
}

sig_files <- list.files(dir_sig, pattern="^enrich_.*_sig\\.csv$", full.names=TRUE)
if (length(sig_files) == 0) {
  message("No filtered enrichment files found in: ", dir_sig, "\nRun enrichment first.")
} else {
  for (ff in sig_files){
    out_csv <- file.path(dir_prio, sub("_sig\\.csv$", "_prioritized_compare.csv", basename(ff), ignore.case=TRUE))
    if (!cfg$overwrite && file.exists(out_csv)) next
    
    raw_dt <- fread(ff, colClasses = list(character="Genes"))
    core   <- score_core(raw_dt)
    
    resA <- apply_profile(core, "A"); setnames(resA, c("priority","rank","tnbc_hits_weighted"), c("priority_A","rank_A","tnbc_hits_A"))
    resB <- apply_profile(core, "B"); setnames(resB, c("priority","rank","tnbc_hits_weighted"), c("priority_B","rank_B","tnbc_hits_B"))
    
    merged <- merge(resA, resB[, .(idx, priority_B, rank_B, tnbc_hits_B)], by="idx", all.x=TRUE)
    keep_cols <- c("Term","pval","padj","combined","Overlap_k","Overlap_N",
                   "neglog10_padj","frac_overlap","spec_comp",
                   names(panel_list),
                   "tnbc_hits_A","priority_A","rank_A",
                   "tnbc_hits_B","priority_B","rank_B")
    keep_cols <- intersect(keep_cols, names(merged))
    out <- merged[, ..keep_cols]; data.table::setorder(out, -priority_B, rank_A)
    fwrite(out, out_csv)
    
    p <- ggplot(head(out, 20L), aes(x = reorder(Term, priority_B), y = priority_B)) +
      geom_col() + coord_flip() +
      labs(x="Pathway", y="Priority (B)", title=paste0("Top20 – ", basename(out_csv))) +
      theme_minimal(base_size = 11)
    ggsave(file.path(dir_prio, sub("\\.csv$", "_top20_B.png", basename(out_csv))), p, width=9, height=7, dpi=150)
  }
  message("\nEnrichment finished.\n- Full:      ", dir_full,
          "\n- Filtered:  ", dir_sig,
          "\n- Prioritized:", dir_prio)
