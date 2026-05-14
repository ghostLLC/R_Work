# ===========================================================================
# US-006: BRCA Differential Expression Analysis (Tumor vs Normal)
# 项目：BRCA多组学数据挖掘
# 日期：2026-05-14
# 方法：DESeq2识别肿瘤vs正常差异表达基因
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
})

cat("\n========== US-006: Differential Expression Analysis ==========\n\n")

INPUT_DIR  <- "D:/Users/Desktop/R_Work/data/processed"
FIG_DIR    <- "D:/Users/Desktop/R_Work/results/figures"
TBL_DIR    <- "D:/Users/Desktop/R_Work/results/tables"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load data ----
cat("Step 1: Loading mRNA count data...\n")

brca <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
counts_all <- brca$counts
cat(sprintf("  Full counts matrix: %d genes x %d samples\n", nrow(counts_all), ncol(counts_all)))

# Get all sample info including normal
sample_info <- read.csv(file.path(INPUT_DIR, "sample_info.csv"))
cat(sprintf("  Sample info: %d samples\n", nrow(sample_info)))
cat(sprintf("  Sample types: %s\n", paste(names(table(sample_info$sample_type)),
      table(sample_info$sample_type), sep="=", collapse=", ")))

# ---- 2. Build tumor vs normal comparison ----
cat("\nStep 2: Building Tumor vs Normal DESeq2 dataset...\n")

# Load full counts (including normal samples)
counts_normal <- readRDS(file.path(INPUT_DIR, "brca_normal_counts.rds"))
counts_tumor  <- readRDS(file.path(INPUT_DIR, "brca_tumor_counts.rds"))

# Take common genes
common_genes <- intersect(rownames(counts_tumor), rownames(counts_normal))
cat(sprintf("  Common genes: %d\n", length(common_genes)))

# Combine tumor + normal
counts_combined <- cbind(
  counts_tumor[common_genes, , drop = FALSE],
  counts_normal[common_genes, , drop = FALSE]
)

# Build colData
col_data <- data.frame(
  sample_id   = colnames(counts_combined),
  condition   = c(rep("Tumor", ncol(counts_tumor)),
                  rep("Normal", ncol(counts_normal))),
  stringsAsFactors = FALSE
)
rownames(col_data) <- col_data$sample_id
col_data$condition <- factor(col_data$condition, levels = c("Normal", "Tumor"))

cat(sprintf("  Combined: %d genes x %d samples (%d Tumor + %d Normal)\n",
            nrow(counts_combined), ncol(counts_combined),
            ncol(counts_tumor), ncol(counts_normal)))

# ---- 3. Run DESeq2 ----
cat("\nStep 3: Running DESeq2...\n")

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_combined),
  colData   = col_data,
  design    = ~ condition
)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= min(10, ncol(dds) / 10)
dds <- dds[keep, ]
cat(sprintf("  After filtering: %d genes\n", sum(keep)))

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)
res <- res[order(res$padj), ]

summary(res)
cat(sprintf("  DEGs (|log2FC|>1 & padj<0.05): %d up, %d down\n",
            sum(res$log2FoldChange > 1 & res$padj < 0.05, na.rm = TRUE),
            sum(res$log2FoldChange < -1 & res$padj < 0.05, na.rm = TRUE)))

# ---- 4. Save results ----
cat("\nStep 4: Saving DEG results...\n")

deg_table <- data.frame(
  gene      = rownames(res),
  baseMean  = res$baseMean,
  log2FC    = res$log2FoldChange,
  lfcSE     = res$lfcSE,
  stat      = res$stat,
  pvalue    = res$pvalue,
  padj      = res$padj,
  stringsAsFactors = FALSE
)
deg_table$regulation <- ifelse(deg_table$padj < 0.05 & deg_table$log2FC > 1, "Up",
                        ifelse(deg_table$padj < 0.05 & deg_table$log2FC < -1, "Down", "NS"))

write.csv(deg_table, file.path(TBL_DIR, "brca_degs_deseq2.csv"), row.names = FALSE)
cat(sprintf("  Saved: brca_degs_deseq2.csv (%d genes)\n", nrow(deg_table)))

# Top DEGs
top_degs <- deg_table %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj) %>%
  head(100)
write.csv(top_degs, file.path(TBL_DIR, "brca_degs_top100.csv"), row.names = FALSE)

# ---- 5. Volcano plot ----
cat("\nStep 5: Generating volcano plot...\n")

pdf(file.path(FIG_DIR, "volcano_brca.pdf"), width = 8, height = 7)

deg_table_plot <- deg_table %>%
  mutate(log10p = -log10(pvalue),
         sig = case_when(
           padj < 0.05 & log2FC > 1  ~ "Up",
           padj < 0.05 & log2FC < -1 ~ "Down",
           TRUE ~ "NS"
         ))

# Cap extreme values for plotting
deg_table_plot$log2FC <- pmax(pmin(deg_table_plot$log2FC, 8), -8)
deg_table_plot$log10p <- pmin(deg_table_plot$log10p, 50)

cols <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")

ggplot(deg_table_plot, aes(x = log2FC, y = log10p, color = sig)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  labs(title = "BRCA: Tumor vs Normal Differential Expression",
       subtitle = sprintf("%d DEGs (|log2FC|>1, padj<0.05)", sum(deg_table_plot$sig != "NS")),
       x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top") +
  annotate("text", x = 5, y = 2, label = paste("Up:", sum(deg_table_plot$sig == "Up")),
           color = "#E41A1C", size = 5) +
  annotate("text", x = -5, y = 2, label = paste("Down:", sum(deg_table_plot$sig == "Down")),
           color = "#377EB8", size = 5)

dev.off()
cat("  Volcano plot saved.\n")

# ---- 6. MA plot ----
cat("\nStep 6: Generating MA plot...\n")

pdf(file.path(FIG_DIR, "ma_plot_brca.pdf"), width = 8, height = 7)
plotMA(res, ylim = c(-5, 5), alpha = 0.05,
       main = "BRCA: MA Plot (Tumor vs Normal)")
dev.off()
cat("  MA plot saved.\n")

# ---- 7. Top DEG heatmap ----
cat("\nStep 7: Generating top DEG heatmap...\n")

top50 <- rownames(res)[order(res$padj)][1:50]
top50 <- top50[!is.na(top50)]

if (length(top50) >= 20) {
  # Use variance-stabilized counts
  vsd <- vst(dds, blind = FALSE)
  vsd_mat <- assay(vsd)[top50, ]

  # Annotate by condition
  annotation_col <- data.frame(Condition = col_data$condition)
  rownames(annotation_col) <- colnames(vsd_mat)
  ann_colors <- list(Condition = c(Tumor = "#E41A1C", Normal = "#377EB8"))

  pdf(file.path(FIG_DIR, "deg_heatmap_top50.pdf"), width = 12, height = 10)
  pheatmap(vsd_mat,
           scale = "row",
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           show_rownames = TRUE,
           show_colnames = FALSE,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           fontsize_row = 6,
           main = "BRCA Top 50 DEGs: Tumor vs Normal")
  dev.off()
  cat("  Heatmap saved.\n")
}

# ---- 8. Summary ----
cat("\n========== US-006 Complete: Differential Expression ==========\n")
cat(sprintf("  Total DEGs (|log2FC|>1, padj<0.05): %d\n", sum(deg_table$regulation != "NS")))
cat(sprintf("  Up-regulated: %d | Down-regulated: %d\n",
            sum(deg_table$regulation == "Up"), sum(deg_table$regulation == "Down")))
