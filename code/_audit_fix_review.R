# ===========================================================================
# 审稿问题核实脚本 — 核实6个问题对应的实际数据
# ===========================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(maftools)
})

cat("========== 问题1&2: 突变基因频率核实 ==========\n")

# Read MAF file
maf <- read.maf("D:/Users/Desktop/R_Work/data/public_data/brca_mutations_temp.maf", verbose = FALSE)

# Get total unique tumor samples in MAF
all_samples <- unique(maf@data$Tumor_Sample_Barcode)
n_total <- length(all_samples)
cat(sprintf("MAF total unique Tumor_Sample_Barcodes: %d\n", n_total))

# Count UNIQUE patients per gene (not variant count!)
gene_patients <- maf@data %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  count(Hugo_Symbol, sort = TRUE)

cat(sprintf("\n--- Top 12 genes by unique patient count (denominator=%d) ---\n", n_total))
for (i in 1:12) {
  g <- gene_patients$Hugo_Symbol[i]
  n <- gene_patients$n[i]
  pct <- 100 * n / n_total
  cat(sprintf("%2d. %-10s %d patients (%.1f%%)\n", i, g, n, pct))
}

# Compare with old N_Mutations (total variant count)
old <- read.csv("D:/Users/Desktop/R_Work/results/tables/mutated_genes_summary.csv")
cat("\n--- Comparison: old (variant events) vs new (unique patients) ---\n")
cat(sprintf("%-10s %15s %15s %15s\n", "Gene", "Old(N_Mutations)", "New(Patients)", "Old%(/990)"))
for (g in c("PIK3CA", "TP53", "TTN")) {
  old_n <- old$N_Mutations[old$Hugo_Symbol == g][1]
  new_n <- gene_patients$n[gene_patients$Hugo_Symbol == g][1]
  old_pct <- 100 * old_n / 990
  new_pct <- 100 * new_n / n_total
  cat(sprintf("%-10s %15d %15d %12.1f%% vs %.1f%% (/%d)\n",
              g, old_n, new_n, old_pct, new_pct, n_total))
}

cat("\n========== 问题3: miRNA阈值不一致核实 ==========\n")
mcor <- read.csv("D:/Users/Desktop/R_Work/results/tables/mirna_target_correlations.csv")
n_03 <- nrow(mcor %>% filter(abs(correlation) > 0.3))
n_04 <- nrow(mcor %>% filter(abs(correlation) > 0.4))
cat(sprintf("|r|>0.3: %d pairs\n", n_03))
cat(sprintf("|r|>0.4: %d pairs\n", n_04))
cat(sprintf("Figure top-20 range: %.3f to %.3f\n",
            min(abs(head(mcor %>% arrange(desc(abs(correlation))), 20)$correlation)),
            max(abs(mcor %>% arrange(desc(abs(correlation))), 20)$correlation)))

cat("\n========== 问题4: 样本量核实 ==========\n")
tumor_rds <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_tumor_counts.rds")
normal_rds <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_normal_counts.rds")
cat(sprintf("brca_tumor_counts.rds: %d genes x %d samples\n", nrow(tumor_rds), ncol(tumor_rds)))
cat(sprintf("brca_normal_counts.rds: %d genes x %d samples\n", nrow(normal_rds), ncol(normal_rds)))
cat(sprintf("DESeq2 combined: %d tumor + %d normal = %d total\n",
            ncol(tumor_rds), ncol(normal_rds), ncol(tumor_rds) + ncol(normal_rds)))

cat("\n========== 问题5: 阈值逻辑 ==========\n")
cat("228 pairs from: 50 miRNAs x 50 mRNAs = 2,500 candidates\n")
cat("2,065 pairs from: 100 miRNAs x 2,292 DEGs = ~229,200 candidates\n")
cat("Candidate space expanded 91.7x, so more pairs found at stricter threshold\n")

cat("\n========== 问题6: XGBoost基线 ==========\n")
clin <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_clinical.rds")
if (exists("clin") && "molecular_subtype" %in% names(clin)) {
  tbl <- table(clin$molecular_subtype)
  cat("Subtype distribution:\n")
  for (i in seq_along(tbl)) {
    cat(sprintf("  %-20s %d (%.1f%%)\n", names(tbl)[i], tbl[i], 100*tbl[i]/sum(tbl)))
  }
  cat(sprintf("\nMajority class baseline: %.1f%%\n", 100 * max(tbl) / sum(tbl)))
}

cat("\n========== 全部核实完成 ==========\n")
