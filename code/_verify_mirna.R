# Verify miRNA correlation numbers are consistent
mirna_cor <- read.csv("D:/Users/Desktop/R_Work/results/tables/mirna_target_correlations.csv")
cat("File:", nrow(mirna_cor), "rows\n")
cat("|r| range:", range(abs(mirna_cor$correlation)), "\n")

# Count at different thresholds
for (th in c(0.3, 0.4, 0.5, 0.6)) {
  n <- sum(abs(mirna_cor$correlation) > th)
  neg <- sum(mirna_cor$correlation < -th)
  pos <- sum(mirna_cor$correlation > th)
  cat(sprintf("|r|>%.1f: %d pairs (%d neg, %d pos)\n", th, n, neg, pos))
}

# Check if this is the 228-pair file or the 2,065-pair file
cat("\nFirst few rows:\n")
print(head(mirna_cor[order(-abs(mirna_cor$correlation)), 1:3], 5))
