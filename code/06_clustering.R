# ===========================================================================
# US-008: BRCA Clustering Analysis (Unsupervised Sample Stratification)
# 项目：BRCA多组学数据挖掘
# 日期：2026-05-14
# 方法：PCA, t-SNE, UMAP + Consensus Clustering
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

cat("\n========== US-008: Clustering Analysis ==========\n\n")

INPUT_DIR  <- "D:/Users/Desktop/R_Work/data/processed"
FIG_DIR    <- "D:/Users/Desktop/R_Work/results/figures"
TBL_DIR    <- "D:/Users/Desktop/R_Work/results/tables"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load data ----
cat("Step 1: Loading mRNA expression data...\n")

brca <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
tpm_mat <- brca$tpm
clinical <- brca$clinical

cat(sprintf("  TPM matrix: %d genes x %d samples\n", nrow(tpm_mat), ncol(tpm_mat)))

# ---- 2. Select top variable genes for clustering ----
cat("\nStep 2: Selecting top variable genes...\n")

log_tpm <- log2(tpm_mat + 1)
gene_var <- apply(log_tpm, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:2000]

log_tpm_top <- log_tpm[top_genes, ]
cat(sprintf("  Using top %d variable genes\n", length(top_genes)))

# ---- 3. PCA Analysis ----
cat("\nStep 3: PCA analysis...\n")

pca <- prcomp(t(log_tpm_top), center = TRUE, scale. = TRUE)
pca_var <- summary(pca)$importance[2, ] * 100

# Prepare plot data
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  sample_id = rownames(pca$x),
  stringsAsFactors = FALSE
)

# Add clinical annotation
expr_patients <- substr(pca_df$sample_id, 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)

match_idx <- match(expr_patients, clinical$pid_short)
pca_df$stage <- clinical$stage_simple[match_idx]
pca_df$subtype <- clinical$molecular_subtype[match_idx]

cat(sprintf("  PC1: %.1f%% | PC2: %.1f%% | PC3: %.1f%% variance\n",
            pca_var[1], pca_var[2], pca_var[3]))

# PCA by subtype
pdf(file.path(FIG_DIR, "pca_subtype.pdf"), width = 9, height = 7)

subtype_colors <- c("Luminal A" = "#1B9E77", "Luminal B" = "#D95F02",
                    "HER2-enriched" = "#7570B3", "Triple Negative" = "#E7298A")

ggplot(pca_df %>% filter(!is.na(subtype)),
       aes(x = PC1, y = PC2, color = subtype)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = subtype_colors, name = "Molecular Subtype") +
  labs(title = "BRCA PCA: Colored by Molecular Subtype",
       x = sprintf("PC1 (%.1f%%)", pca_var[1]),
       y = sprintf("PC2 (%.1f%%)", pca_var[2])) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  stat_ellipse(level = 0.95, linewidth = 1)

dev.off()
cat("  PCA by subtype saved.\n")

# PCA by stage
pdf(file.path(FIG_DIR, "pca_stage.pdf"), width = 9, height = 7)

stage_colors <- c("Stage I" = "#1B9E77", "Stage II" = "#D95F02",
                  "Stage III" = "#7570B3", "Stage IV" = "#E7298A")

ggplot(pca_df %>% filter(!is.na(stage)),
       aes(x = PC1, y = PC2, color = stage)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = stage_colors, name = "Stage") +
  labs(title = "BRCA PCA: Colored by Pathologic Stage",
       x = sprintf("PC1 (%.1f%%)", pca_var[1]),
       y = sprintf("PC2 (%.1f%%)", pca_var[2])) +
  theme_bw(base_size = 14)

dev.off()
cat("  PCA by stage saved.\n")

# ---- 4. t-SNE Analysis ----
cat("\nStep 4: t-SNE dimensionality reduction...\n")

has_tsne <- requireNamespace("Rtsne", quietly = TRUE)

if (has_tsne) {
  library(Rtsne)
  set.seed(42)
  tsne_res <- Rtsne(t(log_tpm_top), perplexity = 30, max_iter = 1000, check_duplicates = FALSE)

  tsne_df <- data.frame(
    tSNE1 = tsne_res$Y[, 1],
    tSNE2 = tsne_res$Y[, 2],
    sample_id = colnames(log_tpm_top),
    subtype = pca_df$subtype,
    stringsAsFactors = FALSE
  )

  pdf(file.path(FIG_DIR, "tsne_brca.pdf"), width = 9, height = 7)
  p <- ggplot(tsne_df %>% filter(!is.na(subtype)),
              aes(x = tSNE1, y = tSNE2, color = subtype)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = subtype_colors, name = "Molecular Subtype") +
    labs(title = "BRCA t-SNE: Top 2000 Variable Genes",
         x = "t-SNE 1", y = "t-SNE 2") +
    theme_bw(base_size = 14)
  print(p)
  dev.off()
  cat("  t-SNE plot saved.\n")
} else {
  cat("  Rtsne not installed, skipping t-SNE.\n")
}

# ---- 5. UMAP Analysis ----
cat("\nStep 5: UMAP dimensionality reduction...\n")

has_umap <- requireNamespace("uwot", quietly = TRUE)

if (has_umap) {
  set.seed(42)
  umap_res <- uwot::umap(t(log_tpm_top), n_neighbors = 30, min_dist = 0.3, n_components = 2)

  umap_df <- data.frame(
    UMAP1 = umap_res[, 1],
    UMAP2 = umap_res[, 2],
    sample_id = colnames(log_tpm_top),
    subtype = pca_df$subtype,
    stringsAsFactors = FALSE
  )

  pdf(file.path(FIG_DIR, "umap_brca.pdf"), width = 9, height = 7)
  p <- ggplot(umap_df %>% filter(!is.na(subtype)),
              aes(x = UMAP1, y = UMAP2, color = subtype)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = subtype_colors, name = "Molecular Subtype") +
    labs(title = "BRCA UMAP: Top 2000 Variable Genes",
         x = "UMAP 1", y = "UMAP 2") +
    theme_bw(base_size = 14)
  print(p)
  dev.off()
  cat("  UMAP plot saved.\n")
} else {
  cat("  uwot not installed, skipping UMAP.\n")
}

# ---- 6. Hierarchical Clustering ----
cat("\nStep 6: Hierarchical clustering with heatmap...\n")

# Use top 500 variable genes for heatmap
top500 <- names(sort(gene_var, decreasing = TRUE))[1:500]
log_tpm_500 <- log_tpm[top500, ]

# Sample annotation
annot_col <- data.frame(
  Subtype = pca_df$subtype,
  Stage = pca_df$stage,
  stringsAsFactors = FALSE
)
rownames(annot_col) <- colnames(log_tpm_500)
annot_col[is.na(annot_col)] <- "Unknown"

pdf(file.path(FIG_DIR, "hclust_heatmap_top500.pdf"), width = 14, height = 10)
pheatmap(log_tpm_500,
         scale = "row",
         annotation_col = annot_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         main = "BRCA Hierarchical Clustering: Top 500 Variable Genes")
dev.off()
cat("  Hierarchical clustering heatmap saved.\n")

# ---- 7. K-means and Silhouette ----
cat("\nStep 7: K-means clustering with silhouette analysis...\n")

sil_scores <- numeric(9)
names(sil_scores) <- as.character(2:10)

set.seed(42)
for (k in 2:10) {
  km <- kmeans(t(log_tpm_top), centers = k, nstart = 25)
  ss <- cluster::silhouette(km$cluster, dist(t(log_tpm_top)))
  sil_scores[as.character(k)] <- mean(ss[, 3])
}

sil_df <- data.frame(K = 2:10, Silhouette = sil_scores, stringsAsFactors = FALSE)
write.csv(sil_df, file.path(TBL_DIR, "clustering_silhouette.csv"), row.names = FALSE)

pdf(file.path(FIG_DIR, "silhouette_scores.pdf"), width = 7, height = 5)
ggplot(sil_df, aes(x = K, y = Silhouette)) +
  geom_line(color = "#377EB8", linewidth = 1) +
  geom_point(size = 3, color = "#E41A1C") +
  labs(title = "BRCA K-means Silhouette Scores",
       x = "Number of Clusters (K)", y = "Average Silhouette Width") +
  theme_bw(base_size = 14) +
  scale_x_continuous(breaks = 2:10)
dev.off()
cat(sprintf("  Best K by silhouette: %d (score=%.4f)\n",
            which.max(sil_scores) + 1, max(sil_scores)))

# ---- 8. Save K-means clusters for best K ----
best_k <- which.max(sil_scores) + 1
set.seed(42)
km_best <- kmeans(t(log_tpm_top), centers = best_k, nstart = 25)

cluster_df <- data.frame(
  sample_id = colnames(log_tpm_top),
  kmeans_cluster = km_best$cluster,
  stringsAsFactors = FALSE
)
write.csv(cluster_df, file.path(TBL_DIR, "clustering_kmeans.csv"), row.names = FALSE)

# Compare clusters with subtypes
if ("subtype" %in% colnames(pca_df)) {
  cluster_comp <- table(cluster_df$kmeans_cluster, pca_df$subtype)
  write.csv(cluster_comp, file.path(TBL_DIR, "clustering_vs_subtype.csv"))
  cat("  Cluster vs Subtype contingency table saved.\n")
}

cat("\n========== US-008 Complete: Clustering Analysis ==========\n")
