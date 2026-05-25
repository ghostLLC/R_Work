# ===========================================================================
# 审稿v2: 核实7个新问题
# ===========================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
})

cat("========== 问题1: Cox Age变量核实 ==========\n")
cox_res <- read.csv("D:/Users/Desktop/R_Work/results/tables/cox_regression.csv")
print(cox_res)

cat("\n========== 问题4: 突变样本量核实 ==========\n")
cat("MAF unique Tumor_Sample_Barcodes: 990\n")
# Check integration summary
int <- read.csv("D:/Users/Desktop/R_Work/results/tables/brca_integration_summary.csv")
print(int)
cat("983 = patients with mRNA+miRNA+clinical+mutation (四层交集)\n")
cat("990 = all patients with mutation data (可能无miRNA或临床)\n")

cat("\n========== 问题5: HER2样本量核实 ==========\n")
clin <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_clinical.rds")
cat(sprintf("Columns: %s\n", paste(names(clin), collapse=", ")))

# Find subtype column
subtype_col <- NULL
for (cname in c("molecular_subtype", "subtype", "PAM50", "Subtype")) {
  if (cname %in% names(clin)) {
    subtype_col <- cname
    break
  }
}
if (!is.null(subtype_col)) {
  tbl <- table(clin[[subtype_col]])
  cat(sprintf("\nSubtype column: %s\n", subtype_col))
  for (i in seq_along(tbl)) {
    cat(sprintf("  %-20s %d (%.1f%%)\n", names(tbl)[i], tbl[i], 100*tbl[i]/sum(tbl)))
  }
  cat(sprintf("Total: %d\n", sum(tbl)))

  # HER2-enriched specifically
  her2_count <- tbl[grepl("HER2", names(tbl), ignore.case=TRUE)]
  cat(sprintf("\nHER2-enriched: %d\n", sum(her2_count)))
}

cat("\n========== 问题6: 重叠率核实 ==========\n")
mRNA_tumor <- 1105
miRNA_total <- 1079
both <- 1076
cat(sprintf("%d/%d = %.1f%%\n", both, mRNA_tumor, 100*both/mRNA_tumor))
cat(sprintf("%d/%d = %.1f%%\n", both, miRNA_total, 100*both/miRNA_total))
cat(sprintf("Neither gives 98.4%%; need to recalculate\n"))

cat("\n========== Done ==========\n")
