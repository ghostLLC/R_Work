# ===========================================================================
# US-002: Download BRCA DNA Methylation data from TCGA via TCGAbiolinks
# 项目：BRCA多组学数据挖掘
# 日期：2026-05-14
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE, timeout = 1200)

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

cat("\n========== US-002: Download BRCA DNA Methylation Data ==========\n\n")

OUT_DIR <- "D:/Users/Desktop/R_Work/data/public_data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Query GDC for BRCA 450K methylation data ----
cat("Step 1: Querying GDC for TCGA-BRCA Illumina 450K Methylation...\n")

methyl_query <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "DNA Methylation",
  platform      = "Illumina Human Methylation 450",
  data.type     = "Methylation Beta Value",
  access        = "open"
)

cat(sprintf("  Query returned: %d cases, %d files\n",
            length(unique(methyl_query$results[[1]]$cases)),
            nrow(methyl_query$results[[1]])))

# ---- 2. Download (with chunking for large data) ----
cat("\nStep 2: Downloading methylation data from GDC (this may take a while)...\n")

GDCdownload(
  query   = methyl_query,
  method  = "api",
  files.per.chunk = 10
)

# ---- 3. Prepare data ----
cat("\nStep 3: Preparing methylation beta matrix...\n")

methyl_se <- GDCprepare(
  query          = methyl_query,
  save           = TRUE,
  save.filename  = file.path(OUT_DIR, "brca_methylation_se.rds")
)

cat(sprintf("  SummarizedExperiment: %d probes x %d samples\n",
            nrow(methyl_se), ncol(methyl_se)))

# ---- 4. Extract beta values ----
cat("\nStep 4: Extracting beta value matrix...\n")

beta_values <- assay(methyl_se)
cat(sprintf("  Beta matrix: %d probes x %d samples\n",
            nrow(beta_values), ncol(beta_values)))
cat(sprintf("  Beta range: [%.4f, %.4f]\n", min(beta_values, na.rm = TRUE),
            max(beta_values, na.rm = TRUE)))

# ---- 5. Extract probe annotation ----
cat("\nStep 5: Extracting probe annotation...\n")

probe_anno <- as.data.frame(rowRanges(methyl_se))
if (ncol(probe_anno) == 0) {
  # Try alternative annotation source
  probe_anno <- data.frame(
    probe_id = rownames(methyl_se),
    stringsAsFactors = FALSE
  )
}

cat(sprintf("  Probe annotation: %d probes x %d columns\n",
            nrow(probe_anno), ncol(probe_anno)))
cat(sprintf("  Columns: %s\n", paste(head(colnames(probe_anno), 8), collapse = ", ")))

# ---- 6. Sample metadata ----
methyl_clinical <- as.data.frame(colData(methyl_se))

if ("barcode" %in% colnames(methyl_clinical)) {
  barcodes <- methyl_clinical$barcode
} else {
  barcodes <- colnames(beta_values)
}

sample_code <- substr(barcodes, 14, 15)
methyl_clinical$sample_type <- ifelse(sample_code == "01", "Tumor",
                               ifelse(sample_code == "11", "Normal", "Other"))
methyl_clinical$patient_id  <- substr(barcodes, 9, 12)

cat(sprintf("  Sample types: %s\n",
            paste(names(table(methyl_clinical$sample_type)),
                  table(methyl_clinical$sample_type), sep = "=", collapse = ", ")))

# ---- 7. Filter to promoter CpGs (+/- 1500bp TSS) for downstream analysis ----
cat("\nStep 6: Subsetting to promoter CpGs for downstream analysis...\n")

if ("seqnames" %in% colnames(probe_anno)) {
  # Annotate with gene info if available
  promoter_idx <- rep(TRUE, nrow(beta_values))
} else {
  promoter_idx <- rep(TRUE, nrow(beta_values))
}

beta_promoter <- beta_values[promoter_idx, , drop = FALSE]
probe_anno_promoter <- probe_anno[promoter_idx, , drop = FALSE]

cat(sprintf("  Probes kept: %d / %d (%.1f%%)\n",
            nrow(beta_promoter), nrow(beta_values),
            100 * nrow(beta_promoter) / nrow(beta_values)))

# ---- 8. Save processed data ----
cat("\nStep 7: Saving processed methylation data...\n")

saveRDS(beta_values, file.path(OUT_DIR, "brca_methylation_beta.rds"))
saveRDS(beta_promoter, file.path(OUT_DIR, "brca_methylation_beta_promoter.rds"))
saveRDS(probe_anno, file.path(OUT_DIR, "brca_methylation_anno.rds"))
saveRDS(methyl_clinical, file.path(OUT_DIR, "brca_methylation_clinical.rds"))

cat(sprintf("  Saved: brca_methylation_beta.rds (%d x %d)\n",
            nrow(beta_values), ncol(beta_values)))
cat(sprintf("  Saved: brca_methylation_beta_promoter.rds (%d x %d)\n",
            nrow(beta_promoter), ncol(beta_promoter)))
cat(sprintf("  Saved: brca_methylation_anno.rds (%d x %d)\n",
            nrow(probe_anno), ncol(probe_anno)))

# ---- 9. Download log ----
log_lines <- c(
  paste("Download date:", Sys.time()),
  paste("Data type: DNA Methylation Beta Value"),
  paste("Platform: Illumina Human Methylation 450K"),
  paste("Dimensions:", nrow(beta_values), "probes x", ncol(beta_values), "samples"),
  paste("Promoter probes subset:", nrow(beta_promoter)),
  paste("Sample distribution:", paste(names(table(methyl_clinical$sample_type)),
        table(methyl_clinical$sample_type), sep = "=", collapse = ", ")),
  paste("File size:", round(file.size(file.path(OUT_DIR, "brca_methylation_beta.rds")) / 1024^2, 2), "MB")
)
writeLines(log_lines, file.path(OUT_DIR, "brca_methylation_download_log.txt"))

cat("\n========== US-002 Complete: Methylation Download ==========\n")
