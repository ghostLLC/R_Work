# Fix miRNA figure: change subtitle from |r|>0.4 to |r|>0.3, 228 pairs
suppressPackageStartupMessages({library(tidyverse); library(ggsci)})
npg10 <- pal_npg("nrc")(10)
mcor <- read.csv("D:/Users/Desktop/R_Work/results/tables/mirna_target_correlations.csv")
cat(sprintf("Total pairs: %d\n", nrow(mcor)))
cat(sprintf("|r|>0.3: %d pairs\n", nrow(mcor %>% filter(abs(correlation) > 0.3))))
cat(sprintf("|r|>0.4: %d pairs\n", nrow(mcor %>% filter(abs(correlation) > 0.4))))

tp <- mcor %>%
  arrange(desc(abs(correlation))) %>%
  head(20) %>%
  mutate(miRNA_s = gsub("^read_count_|^hsa-", "", miRNA),
         mRNA_s  = gsub("_ENSG\\d+$", "", mRNA))

OUT <- "D:/Users/Desktop/R_Work/results/figures_pub"

p10 <- ggplot(tp, aes(mRNA_s, miRNA_s)) +
  geom_tile(aes(fill = correlation), color = "white", linewidth = 1) +
  scale_fill_gradient2(low = npg10[5], mid = "white", high = npg10[3],
                       midpoint = 0, name = "r") +
  geom_text(aes(label = sprintf("%.2f", correlation)),
            size = 5, fontface = "bold") +
  labs(title = "miRNA-mRNA Regulatory Pairs (Top 20)",
       subtitle = paste0("|r| > 0.3  |  Total: ",
                         nrow(mcor %>% filter(abs(correlation) > 0.3)), " pairs"),
       x = "mRNA", y = "miRNA") +
  theme_minimal(base_size = 17) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title    = element_text(size = 21, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30"),
    legend.title  = element_text(size = 18),
    legend.text   = element_text(size = 16),
    panel.grid    = element_blank()
  )

ggsave(file.path(OUT, "fig10_mirna_mrna_corr.png"), p10,
       width = 13, height = 10, dpi = 300, bg = "white")
cat("  Saved: fig10_mirna_mrna_corr.png\n")

# Also verify sample sizes
cat("\n=== Sample size verification ===\n")
tumor_rds <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_tumor_counts.rds")
normal_rds <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_normal_counts.rds")
cat(sprintf("brca_tumor_counts.rds: %d genes x %d samples\n", nrow(tumor_rds), ncol(tumor_rds)))
cat(sprintf("brca_normal_counts.rds: %d genes x %d samples\n", nrow(normal_rds), ncol(normal_rds)))
cat(sprintf("DESeq2 combined: %d tumor + %d normal = %d\n", ncol(tumor_rds), ncol(normal_rds), ncol(tumor_rds)+ncol(normal_rds)))

# Check clinical subtypes for XGBoost baseline
clin <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_clinical.rds")
if ("molecular_subtype" %in% names(clin)) {
  tbl <- table(clin$molecular_subtype)
  cat("\nSubtype distribution:\n")
  for (i in seq_along(tbl)) {
    cat(sprintf("  %-20s %d (%.1f%%)\n", names(tbl)[i], tbl[i], 100*tbl[i]/sum(tbl)))
  }
  cat(sprintf("Majority class baseline: %.1f%%\n", 100 * max(tbl) / sum(tbl)))
}

cat("\nDone.\n")
