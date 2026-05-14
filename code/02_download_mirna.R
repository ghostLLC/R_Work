# ===========================================================================
# US-001: Download BRCA miRNA-seq data from TCGA via TCGAbiolinks
# 项目：BRCA多组学数据挖掘
# 日期：2026-05-14
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE, timeout = 1200)

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

cat("\n========== US-001: Download BRCA miRNA-seq Data ==========\n\n")

OUT_DIR <- "D:/Users/Desktop/R_Work/data/public_data"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Query GDC for BRCA miRNA data ----
cat("Step 1: Querying GDC for TCGA-BRCA miRNA-seq data...\n")
cat("  Data Category: Transcriptome Profiling\n")
cat("  Data Type: miRNA Expression Quantification\n")

mirna_query <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  access        = "open"
)

cat(sprintf("  Query returned: %d cases, %d files\n",
            length(unique(mirna_query$results[[1]]$cases)),
            nrow(mirna_query$results[[1]])))

# ---- 2. Download data ----
cat("\nStep 2: Downloading miRNA data from GDC...\n")

GDCdownload(
  query   = mirna_query,
  method  = "api",
  files.per.chunk = 20
)

# ---- 3. Prepare data ----
cat("\nStep 3: Preparing miRNA expression matrix...\n")

mirna_data <- GDCprepare(
  query   = mirna_query,
  save           = TRUE,
  save.filename  = file.path(OUT_DIR, "brca_mirna_se.rds")
)

cat(sprintf("  Data dimensions: %d rows x %d columns\n",
            nrow(mirna_data), ncol(mirna_data)))
cat(sprintf("  Class: %s\n", paste(class(mirna_data), collapse = ", ")))

# ---- 4. Extract count matrix ----
# GDC miRNA data may return data.frame or SummarizedExperiment
cat("\nStep 4: Extracting read count matrix...\n")

if (inherits(mirna_data, "SummarizedExperiment")) {
  mirna_counts <- assay(mirna_data, "read_count")
  if (is.null(mirna_counts)) {
    avail_assays <- assayNames(mirna_data)
    cat(sprintf("  Available assays: %s\n", paste(avail_assays, collapse = ", ")))
    mirna_counts <- assay(mirna_data, 1)
  }
  mirna_clinical <- as.data.frame(colData(mirna_data))
  barcodes <- colnames(mirna_counts)

} else if (is.data.frame(mirna_data)) {
  # Data.frame format: col1=miRNA_ID, then 3 cols per sample (read_count, RPM, cross-mapped)
  cat("  Data is data.frame — extracting read_count columns\n")

  mirna_id_col <- mirna_data[, 1, drop = TRUE]
  all_cols <- colnames(mirna_data)[-1]

  # Each sample has 3 quantification types; extract read_count columns
  read_count_cols <- grep("read_count", all_cols, ignore.case = TRUE, value = TRUE)

  if (length(read_count_cols) == 0) {
    # Alternative: every 3rd column starting from col 2
    sample_cols <- all_cols[seq(1, length(all_cols), by = 3)]
    mirna_counts <- as.matrix(mirna_data[, sample_cols, drop = FALSE])
    colnames(mirna_counts) <- gsub("_read_count|_reads_per_million|_cross.mapped", "", sample_cols)
  } else {
    mirna_counts <- as.matrix(mirna_data[, read_count_cols, drop = FALSE])
    colnames(mirna_counts) <- gsub("_read_count", "", read_count_cols)
  }

  rownames(mirna_counts) <- mirna_id_col

  # Build clinical from barcodes
  barcodes <- colnames(mirna_counts)
  mirna_clinical <- data.frame(
    barcode = barcodes,
    stringsAsFactors = FALSE
  )

} else {
  stop("Unexpected data format from GDCprepare")
}

cat(sprintf("  Count matrix: %d miRNAs x %d samples\n",
            nrow(mirna_counts), ncol(mirna_counts)))

# ---- 5. Parse TCGA barcode ----
cat("\nStep 5: Parsing sample metadata...\n")

sample_code  <- substr(barcodes, 14, 15)
mirna_clinical$sample_type <- ifelse(sample_code == "01", "Tumor",
                              ifelse(sample_code == "11", "Normal", "Other"))
mirna_clinical$patient_id  <- substr(barcodes, 9, 12)

cat(sprintf("  Sample types: %s\n",
            paste(names(table(mirna_clinical$sample_type)),
                  table(mirna_clinical$sample_type), sep = "=", collapse = ", ")))

# ---- 6. Save processed data ----
cat("\nStep 5: Saving processed miRNA data...\n")

saveRDS(mirna_counts, file.path(OUT_DIR, "brca_mirna_counts.rds"))
saveRDS(mirna_clinical, file.path(OUT_DIR, "brca_mirna_clinical.rds"))

cat(sprintf("  Saved: brca_mirna_counts.rds (%d x %d)\n",
            nrow(mirna_counts), ncol(mirna_counts)))
cat(sprintf("  Saved: brca_mirna_clinical.rds (%d x %d)\n",
            nrow(mirna_clinical), ncol(mirna_clinical)))

# ---- 7. Download log ----
log_lines <- c(
  paste("Download date:", Sys.time()),
  paste("Data type: miRNA Expression Quantification"),
  paste("Platform: Illumina HiSeq (BCGSC miRNA Profiling)"),
  paste("Dimensions:", nrow(mirna_counts), "miRNAs x", ncol(mirna_counts), "samples"),
  paste("Sample distribution:", paste(names(table(mirna_clinical$sample_type)),
        table(mirna_clinical$sample_type), sep = "=", collapse = ", ")),
  paste("File size:", round(file.size(file.path(OUT_DIR, "brca_mirna_counts.rds")) / 1024^2, 2), "MB")
)
writeLines(log_lines, file.path(OUT_DIR, "brca_mirna_download_log.txt"))

cat("\n========== US-001 Complete: miRNA Download ==========\n")
