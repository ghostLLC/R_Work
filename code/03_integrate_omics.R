# ===========================================================================
# US-004: BRCA Multi-Omics Data Integration by Patient ID
# 项目：BRCA多组学数据挖掘
# 日期：2026-05-14
# 说明：对齐mRNA、miRNA数据按TCGA患者ID，建立统一样本映射表
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("\n========== US-004: Multi-Omics Integration ==========\n\n")

INPUT_DIR  <- "D:/Users/Desktop/R_Work/data/processed"
MIRNA_DIR  <- "D:/Users/Desktop/R_Work/data/public_data"
OUTPUT_DIR <- "D:/Users/Desktop/R_Work/data/processed"
FIG_DIR    <- "D:/Users/Desktop/R_Work/results/figures"
TBL_DIR    <- "D:/Users/Desktop/R_Work/results/tables"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load mRNA data ----
cat("Step 1: Loading mRNA expression data...\n")
brca <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
mRNA_tpm    <- brca$tpm
mRNA_clinical <- brca$clinical
mRNA_anno   <- brca$gene_annotation

cat(sprintf("  mRNA TPM: %d genes x %d samples\n", nrow(mRNA_tpm), ncol(mRNA_tpm)))

# Extract patient IDs from mRNA column names (TCGA barcode)
mRNA_barcodes <- colnames(mRNA_tpm)
mRNA_patient_ids <- substr(mRNA_barcodes, 9, 12)
mRNA_patient_full <- paste0("TCGA-", substr(mRNA_barcodes, 6, 7), "-", mRNA_patient_ids)

mRNA_patients <- data.frame(
  barcode    = mRNA_barcodes,
  patient_id = mRNA_patient_ids,
  patient_full = mRNA_patient_full,
  sample_type = substr(mRNA_barcodes, 14, 15),
  omics_mRNA = TRUE,
  stringsAsFactors = FALSE
)
cat(sprintf("  mRNA unique patients: %d\n", length(unique(mRNA_patient_ids))))

# ---- 2. Load miRNA data ----
cat("\nStep 2: Loading miRNA expression data...\n")
mirna_counts <- readRDS(file.path(MIRNA_DIR, "brca_mirna_counts.rds"))
mirna_clin   <- readRDS(file.path(MIRNA_DIR, "brca_mirna_clinical.rds"))

cat(sprintf("  miRNA counts: %d miRNAs x %d samples\n", nrow(mirna_counts), ncol(mirna_counts)))

# Extract patient IDs from miRNA column names (strip "read_count_" prefix)
mirna_barcodes <- colnames(mirna_counts)
# Remove prefix like "read_count_" if present
mirna_barcodes <- gsub("^read_count_", "", mirna_barcodes)
mirna_patient_ids <- substr(mirna_barcodes, 9, 12)

mirna_patients <- data.frame(
  barcode    = mirna_barcodes,
  patient_id = mirna_patient_ids,
  sample_type = substr(mirna_barcodes, 14, 15),
  omics_miRNA = TRUE,
  stringsAsFactors = FALSE
)
cat(sprintf("  miRNA unique patients: %d\n", length(unique(mirna_patient_ids))))

# ---- 3. Build unified sample mapping ----
cat("\nStep 3: Building unified multi-omics sample map...\n")

# Get all unique patient IDs
all_patients <- sort(unique(c(mRNA_patient_ids, mirna_patient_ids)))

sample_map <- data.frame(
  patient_id = all_patients,
  mRNA  = all_patients %in% mRNA_patient_ids,
  miRNA = all_patients %in% mirna_patient_ids,
  stringsAsFactors = FALSE
)

# Count sample distribution
mRNA_tumor_ids  <- unique(mRNA_patients$patient_id[mRNA_patients$sample_type == "01"])
miRNA_tumor_ids <- unique(mirna_patients$patient_id[mirna_patients$sample_type == "01"])

cat(sprintf("  Total unique patients across all omics: %d\n", nrow(sample_map)))
cat(sprintf("  Patients with mRNA:  %d (tumor: %d)\n",
            sum(sample_map$mRNA), length(mRNA_tumor_ids)))
cat(sprintf("  Patients with miRNA: %d (tumor: %d)\n",
            sum(sample_map$miRNA), length(miRNA_tumor_ids)))
cat(sprintf("  Patients with BOTH:  %d\n", sum(sample_map$mRNA & sample_map$miRNA)))

# ---- 4. Create aligned multi-omics matrices (patient-level) ----
cat("\nStep 4: Creating patient-aligned multi-omics matrices...\n")

# For patients with both mRNA and miRNA data
common_patients <- sample_map$patient_id[sample_map$mRNA & sample_map$miRNA]
cat(sprintf("  Common patients (both mRNA + miRNA): %d\n", length(common_patients)))

# Build aligned matrices using first tumor sample per patient
get_patient_matrix <- function(counts, patient_ids, target_patients) {
  # For each target patient, take first matching column
  result <- matrix(NA, nrow = nrow(counts), ncol = length(target_patients))
  rownames(result) <- rownames(counts)
  colnames(result) <- target_patients

  for (i in seq_along(target_patients)) {
    pid <- target_patients[i]
    cols <- which(patient_ids == pid)
    if (length(cols) > 0) {
      result[, i] <- counts[, cols[1]]
    }
  }
  return(result)
}

mRNA_aligned <- get_patient_matrix(mRNA_tpm, mRNA_patient_ids, common_patients)
miRNA_aligned <- get_patient_matrix(mirna_counts, mirna_patient_ids, common_patients)

cat(sprintf("  Aligned mRNA matrix:  %d genes x %d patients\n",
            nrow(mRNA_aligned), ncol(mRNA_aligned)))
cat(sprintf("  Aligned miRNA matrix: %d miRNAs x %d patients\n",
            nrow(miRNA_aligned), ncol(miRNA_aligned)))

# ---- 5. Generate integration statistics ----
cat("\nStep 5: Computing integration summary...\n")

integration_summary <- data.frame(
  Omics_Layer = c("mRNA", "miRNA", "Both"),
  N_Samples   = c(sum(sample_map$mRNA), sum(sample_map$miRNA),
                  sum(sample_map$mRNA & sample_map$miRNA)),
  N_Patients_Unique = c(length(unique(mRNA_patient_ids)),
                        length(unique(mirna_patient_ids)),
                        length(common_patients)),
  Data_Type = c("RNA-seq TPM", "miRNA-seq read_count", "Aligned multi-omics"),
  Dimensions = c(sprintf("%d genes x %d samples", nrow(mRNA_tpm), ncol(mRNA_tpm)),
                 sprintf("%d miRNAs x %d samples", nrow(mirna_counts), ncol(mirna_counts)),
                 sprintf("mRNA:%dx%d + miRNA:%dx%d", nrow(mRNA_aligned),
                         ncol(mRNA_aligned), nrow(miRNA_aligned), ncol(miRNA_aligned))),
  stringsAsFactors = FALSE
)

write.csv(integration_summary, file.path(TBL_DIR, "brca_integration_summary.csv"), row.names = FALSE)
cat("  Integration summary saved.\n")

# ---- 6. Save aligned data ----
cat("\nStep 6: Saving aligned multi-omics data...\n")

saveRDS(sample_map, file.path(OUTPUT_DIR, "brca_multimics_sample_map.rds"))
saveRDS(mRNA_aligned, file.path(OUTPUT_DIR, "brca_mRNA_aligned.rds"))
saveRDS(miRNA_aligned, file.path(OUTPUT_DIR, "brca_miRNA_aligned.rds"))

write.csv(sample_map, file.path(OUTPUT_DIR, "brca_multimics_sample_map.csv"), row.names = FALSE)

cat("  Saved: brca_multimics_sample_map.rds/csv\n")
cat("  Saved: brca_mRNA_aligned.rds\n")
cat("  Saved: brca_miRNA_aligned.rds\n")

# ---- 7. UpSet plot of omics overlap ----
cat("\nStep 7: Generating omics overlap visualization...\n")

# Use UpSetR if available, otherwise base R
if (requireNamespace("UpSetR", quietly = TRUE)) {
  library(UpSetR)

  # Build binary matrix for UpSet
  upset_data <- as.data.frame(sample_map[, c("mRNA", "miRNA")])
  upset_data[upset_data == TRUE] <- 1
  upset_data[upset_data == FALSE] <- 0
  colnames(upset_data) <- c("mRNA (Transcriptome)", "miRNA (miRNA-seq)")

  pdf(file.path(FIG_DIR, "omics_overlap_upset.pdf"), width = 8, height = 5)
  upset(upset_data,
        sets = colnames(upset_data),
        keep.order = TRUE,
        sets.bar.color = c("#E41A1C", "#377EB8"),
        mainbar.y.label = "Patient Intersection Size",
        sets.x.label = "Patients per Omics Layer")
  dev.off()
  cat("  UpSet plot saved.\n")

} else {
  # Fallback: simple venn-like barplot
  cat("  UpSetR not installed, using base R barplot...\n")

  pdf(file.path(FIG_DIR, "omics_overlap_bar.pdf"), width = 7, height = 5)
  counts <- c(
    mRNA_Only  = sum(sample_map$mRNA & !sample_map$miRNA),
    miRNA_Only = sum(!sample_map$mRNA & sample_map$miRNA),
    Both       = sum(sample_map$mRNA & sample_map$miRNA)
  )
  bp <- barplot(counts, col = c("#E41A1C", "#377EB8", "#4DAF4A"),
                main = "BRCA Multi-Omics Sample Overlap",
                ylab = "Number of Patients", ylim = c(0, max(counts) * 1.2))
  text(bp, counts + max(counts) * 0.03, labels = counts, cex = 1.2)
  dev.off()
  cat("  Bar plot saved.\n")
}

# ---- 8. Generate patient-level multi-omics overview ----
cat("\nStep 8: Multi-omics data overview...\n")

omics_overview <- data.frame(
  Metric = c(
    "Total unique patients",
    "Patients with mRNA data",
    "Patients with miRNA data",
    "Patients with both mRNA + miRNA",
    "Patients with mRNA only",
    "Patients with miRNA only",
    "mRNA features (genes)",
    "miRNA features",
    "Data type: mRNA",
    "Data type: miRNA",
    "Integration strategy"
  ),
  Value = c(
    nrow(sample_map),
    sum(sample_map$mRNA),
    sum(sample_map$miRNA),
    sum(sample_map$mRNA & sample_map$miRNA),
    sum(sample_map$mRNA & !sample_map$miRNA),
    sum(!sample_map$mRNA & sample_map$miRNA),
    nrow(mRNA_tpm),
    nrow(mirna_counts),
    "RNA-seq TPM (tumor samples)",
    "miRNA-seq read counts",
    "Patient-level alignment by TCGA patient ID (chars 9-12 of barcode)"
  ),
  stringsAsFactors = FALSE
)

write.csv(omics_overview, file.path(TBL_DIR, "brca_omics_overview.csv"), row.names = FALSE)
write.csv(sample_map, file.path(OUTPUT_DIR, "brca_multimics_sample_map.csv"), row.names = FALSE)

cat("\n========== US-004 Complete: Multi-Omics Integration ==========\n")
cat(sprintf("  Common patients: %d | mRNA only: %d | miRNA only: %d\n",
            sum(sample_map$mRNA & sample_map$miRNA),
            sum(sample_map$mRNA & !sample_map$miRNA),
            sum(!sample_map$mRNA & sample_map$miRNA)))
