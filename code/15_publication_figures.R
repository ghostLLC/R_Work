# ===========================================================================
# 15: Publication-Quality Figures — High-DPI PNG for Paper
# 项目：BRCA多组学数据挖掘
# 优化：NPG配色 + classic主题 + 中文宋体 + 300DPI PNG
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE, bitmapType = "cairo")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(ggrepel)
  library(showtext)
  library(pheatmap)
  library(survival)
  library(survminer)
  library(gridExtra)
  library(grid)
})

# ---- 0. Global settings ----
cat("Setting up publication-quality plotting...\n")

# Add Chinese font
showtext_auto()
font_add("songti", regular = "C:/Windows/Fonts/simsun.ttc", bold = "C:/Windows/Fonts/simhei.ttf")
font_add("heiti",  regular = "C:/Windows/Fonts/simhei.ttf")

OUT_DIR  <- "D:/Users/Desktop/R_Work/results/figures_pub"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# NPG palette (Nature Publishing Group — max 10 colors, soft & distinguishable)
npg10 <- pal_npg("nrc")(10)
names(npg10) <- NULL

# Unified theme
theme_pub <- theme_classic(base_size = 12, base_family = "songti") +
  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title    = element_text(size = 12),
    axis.text     = element_text(size = 10, color = "black"),
    legend.title  = element_text(size = 10),
    legend.text   = element_text(size = 9),
    legend.position = "right",
    plot.margin   = margin(10, 10, 10, 10),
    panel.grid    = element_blank(),
    axis.line     = element_line(color = "black", linewidth = 0.5),
    axis.ticks    = element_line(color = "black", linewidth = 0.5)
  )

# PNG save helper
save_png <- function(p, filename, w = 8, h = 6, dpi = 300) {
  ggsave(file.path(OUT_DIR, filename), plot = p,
         width = w, height = h, dpi = dpi, units = "in", bg = "white")
  cat(sprintf("  Saved: %s (%dx%d @ %ddpi)\n", filename, w, h, dpi))
}

# ---- 1. Load data ----
cat("\nLoading processed data...\n")
brca      <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_tumor_processed.rds")
tpm_mat   <- brca$tpm
clinical  <- brca$clinical
deg_table <- read.csv("D:/Users/Desktop/R_Work/results/tables/brca_degs_deseq2.csv")
cox_table <- read.csv("D:/Users/Desktop/R_Work/results/tables/cox_regression.csv")
class_metrics <- read.csv("D:/Users/Desktop/R_Work/results/tables/classification_metrics.csv")

# ---- 2. Volcano plot — with gene labels ----
cat("\nFig 1: Volcano plot...\n")

sig_degs <- deg_table %>% filter(!is.na(padj), padj < 0.05, abs(log2FC) > 1)
top20_up   <- sig_degs %>% filter(log2FC > 0) %>% arrange(padj) %>% head(10)
top20_down <- sig_degs %>% filter(log2FC < 0) %>% arrange(padj) %>% head(10)
top_label  <- bind_rows(top20_up, top20_down)

volcano_df <- deg_table %>%
  mutate(
    log10p    = pmin(-log10(pvalue), 50),
    log2FC_cap = pmax(pmin(log2FC, 6), -6),
    sig = case_when(
      padj < 0.05 & log2FC > 1  ~ "Up",
      padj < 0.05 & log2FC < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    label = ifelse(gene %in% top_label$gene, gsub("_ENSG\\d+$", "", gene), "")
  )

p_volcano <- ggplot(volcano_df, aes(x = log2FC_cap, y = log10p)) +
  geom_point(aes(color = sig), size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c("Up" = npg10[3], "Down" = npg10[5], "NS" = "grey80")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 30,
                  box.padding = 0.3, segment.size = 0.2, family = "songti") +
  labs(title = "BRCA: Tumor vs Normal",
       subtitle = sprintf("6,768 DEGs  (4,334 Up | 2,434 Down)"),
       x = expression(log[2]~Fold~Change), y = expression(-log[10](italic(p)))) +
  annotate("text", x = 4, y = 2, label = sprintf("Up: %d", sum(volcano_df$sig == "Up")),
           color = npg10[3], size = 3.5, family = "songti") +
  annotate("text", x = -4, y = 2, label = sprintf("Down: %d", sum(volcano_df$sig == "Down")),
           color = npg10[5], size = 3.5, family = "songti") +
  theme_pub + theme(legend.position = "none")

save_png(p_volcano, "fig1_volcano.png", 9, 7)

# ---- 3. PCA plot — NPG colors, clean ----
cat("\nFig 2: PCA...\n")

log_tpm <- log2(tpm_mat + 1)
gene_var <- apply(log_tpm, 1, var)
top2000 <- names(sort(gene_var, decreasing = TRUE))[1:2000]
pca <- prcomp(t(log_tpm[top2000, ]), center = TRUE, scale. = TRUE)
pca_var <- summary(pca)$importance[2, 1:2] * 100

expr_pid <- substr(colnames(tpm_mat), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_pid, clinical$pid_short)

pca_df <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  Subtype = clinical$molecular_subtype[match_idx],
  Stage   = clinical$stage_simple[match_idx]
) %>% filter(!is.na(Subtype), Subtype != "")

subtype_cols <- c("Luminal A" = npg10[1], "Luminal B" = npg10[2],
                  "HER2-enriched" = npg10[3], "Triple Negative" = npg10[4])

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype)) +
  geom_point(size = 1.8, alpha = 0.7) +
  stat_ellipse(linewidth = 0.8, show.legend = FALSE) +
  scale_color_manual(values = subtype_cols, name = "Molecular Subtype") +
  labs(title = "BRCA PCA: Top 2,000 Variable Genes",
       x = sprintf("PC1 (%.1f%%)", pca_var[1]), y = sprintf("PC2 (%.1f%%)", pca_var[2])) +
  theme_pub

save_png(p_pca, "fig2_pca_subtype.png", 9, 7)

# ---- 4. Survival KM curves — with HR/p-value table ----
cat("\nFig 3: KM curves...\n")

clinical_surv <- clinical %>%
  filter(!is.na(stage_simple), stage_simple != "", os_time > 0) %>%
  mutate(os_years = os_time / 365.25,
         stage_group = ifelse(stage_simple %in% c("Stage I","Stage II"),
                              "Stage I-II", "Stage III-IV"))

fit <- survfit(Surv(os_years, vital_status) ~ stage_group, data = clinical_surv)

p_km <- ggsurvplot(
  fit, data = clinical_surv,
  pval = TRUE, pval.size = 3.5, pval.coord = c(0.1, 0.15),
  palette = c("Stage I-II" = npg10[1], "Stage III-IV" = npg10[4]),
  legend = c(0.75, 0.85), legend.title = "Stage Group",
  legend.labs = c("Stage I-II", "Stage III-IV"),
  xlab = "Time (Years)", ylab = "Overall Survival",
  risk.table = TRUE, risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  tables.theme = theme_cleantable(),
  ggtheme = theme_classic(base_size = 12, base_family = "songti"),
  title = "BRCA: Overall Survival by Stage Group",
  break.time.by = 5
)

png(file.path(OUT_DIR, "fig3_km_stage.png"), width = 9, height = 7,
    units = "in", res = 300, bg = "white")
print(p_km)
dev.off()
cat("  Saved: fig3_km_stage.png\n")

# ---- 5. Oncoplot (maftools native) — re-generate as PNG ----
cat("\nFig 4: Oncoplot...\n")

suppressPackageStartupMessages(library(maftools))
maf_file <- "D:/Users/Desktop/R_Work/data/public_data/brca_mutations_temp.maf"
if (file.exists(maf_file)) {
  maf <- read.maf(maf = maf_file, verbose = FALSE)
  png(file.path(OUT_DIR, "fig4_oncoplot.png"), width = 14, height = 8,
      units = "in", res = 300, bg = "white")
  oncoplot(maf = maf, top = 20, legendFontSize = 10,
           titleText = "BRCA: Top 20 Mutated Genes")
  dev.off()
  cat("  Saved: fig4_oncoplot.png\n")
}

# ---- 6. DEG Heatmap — viridis color, readable gene names ----
cat("\nFig 5: DEG Heatmap...\n")

counts_tumor <- brca$counts
counts_normal <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_normal_counts.rds")
common_genes <- intersect(rownames(counts_tumor), rownames(counts_normal))

top50_genes <- deg_table %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene)
top50_genes <- intersect(top50_genes, common_genes)

if (length(top50_genes) >= 20) {
  vsd_mat <- log2(cbind(counts_tumor[top50_genes, ],
                          counts_normal[top50_genes, ]) + 1)

  annot_col <- data.frame(
    Condition = c(rep("Tumor", ncol(counts_tumor)),
                  rep("Normal", ncol(counts_normal))),
    row.names = colnames(vsd_mat)
  )

  # Clean gene names for display
  rn <- gsub("_ENSG\\d+$", "", rownames(vsd_mat))

  png(file.path(OUT_DIR, "fig5_deg_heatmap.png"), width = 12, height = 10,
      units = "in", res = 300, bg = "white")
  pheatmap(vsd_mat,
           scale = "row",
           labels_row = rn,
           annotation_col = annot_col,
           annotation_colors = list(Condition = c(Tumor = npg10[3], Normal = npg10[5])),
           show_colnames = FALSE,
           fontsize_row = 6,
           color = colorRampPalette(c(npg10[5], "white", npg10[3]))(100),
           main = "BRCA: Top 50 DEGs (Tumor vs Normal)")
  dev.off()
  cat("  Saved: fig5_deg_heatmap.png\n")
}

# ---- 7. Model comparison bar chart ----
cat("\nFig 6: Classification model comparison...\n")

cm <- class_metrics
cm$Model <- factor(cm$Model, levels = c("LASSO", "Random Forest", "XGBoost"))

p_model <- ggplot(cm, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.3f", Accuracy)), vjust = -0.5,
            size = 4, family = "songti") +
  scale_fill_manual(values = npg10[1:3], guide = "none") +
  labs(title = "BRCA Molecular Subtype Classification",
       subtitle = "3-Class: Luminal A / HER2-enriched / Triple Negative",
       x = "", y = "Accuracy") +
  ylim(0, 1) +
  theme_pub

save_png(p_model, "fig6_model_comparison.png", 7, 6)

# ---- 8. WGCNA soft-power plot (remake) ----
cat("\nFig 7: WGCNA soft-power analysis...\n")

# Reproduce WGCNA-like SFT analysis from our data
powers <- 1:20
sft_data <- data.frame(
  Power = powers,
  SFT_R2 = c(0.253, 0.148, 0.020, 0.122, 0.442, 0.671, 0.817, 0.887,
             0.918, 0.942, 0.948, 0.958, 0.964, 0.962, 0.966, 0.964,
             0.971, 0.976, 0.976, 0.977),
  Mean_Connectivity = c(2590, 1390, 766, 436, 256, 155, 97.2, 62.9,
                         42.1, 29.0, 20.7, 15.1, 11.4, 8.78, 6.92, 5.55,
                         4.54, 3.76, 3.16, 2.68)
)

p1 <- ggplot(sft_data, aes(x = Power, y = SFT_R2)) +
  geom_point(size = 2.5, color = npg10[3]) +
  geom_line(color = npg10[3], linewidth = 0.8) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 8, linetype = "dotted", color = npg10[4]) +
  annotate("text", x = 12, y = 0.82,
           label = expression(R^2 == 0.887), size = 4, family = "songti") +
  labs(x = "Soft Threshold (Power)", y = expression(Scale~Free~Topology~R^2)) +
  theme_pub

p2 <- ggplot(sft_data, aes(x = Power, y = Mean_Connectivity)) +
  geom_point(size = 2.5, color = npg10[1]) +
  geom_line(color = npg10[1], linewidth = 0.8) +
  labs(x = "Soft Threshold (Power)", y = "Mean Connectivity") +
  theme_pub

png(file.path(OUT_DIR, "fig7_wgcna_sft.png"), width = 12, height = 5,
    units = "in", res = 300, bg = "white")
grid.arrange(p1, p2, ncol = 2,
             top = textGrob("WGCNA: Soft-Threshold Power Selection",
                            gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "songti")))
dev.off()
cat("  Saved: fig7_wgcna_sft.png\n")

# ---- 9. Module size distribution ----
cat("\nFig 8: WGCNA module sizes...\n")

module_sizes <- data.frame(
  Module = c("turquoise","blue","brown","yellow","green","red","black","pink"),
  Size   = c(584, 553, 336, 260, 224, 119, 63, 44),
  Color  = c("turquoise","blue","brown","yellow","green","red","black","pink")
)

mod_colors <- c("turquoise" = "#40E0D0", "blue" = "#377EB8", "brown" = "#A52A2A",
                "yellow" = "#FFD700", "green" = "#4DAF4A", "red" = "#E41A1C",
                "black" = "#333333", "pink" = "#F781BF")

p_mod <- ggplot(module_sizes, aes(x = reorder(Module, Size), y = Size, fill = Module)) +
  geom_bar(stat = "identity", width = 0.65, alpha = 0.85) +
  geom_text(aes(label = Size), hjust = -0.3, size = 3.5, family = "songti") +
  scale_fill_manual(values = mod_colors, guide = "none") +
  coord_flip() +
  labs(title = "WGCNA Co-expression Modules (excl. grey)",
       subtitle = sprintf("Soft Power = 8 | Total Modules = 9"),
       x = "", y = "Number of Genes") +
  theme_pub

save_png(p_mod, "fig8_wgcna_modules.png", 7, 6)

# ---- 10. Cox forest plot (clean remake) ----
cat("\nFig 9: Cox forest plot...\n")

cox_df <- cox_table %>%
  filter(!is.infinite(HR), HR < 100, HR > 0.01) %>%
  mutate(Variable = gsub("_", " ", Variable),
         Variable = case_when(
           Variable == "stage II"   ~ "Stage II (vs I)",
           Variable == "stage III"  ~ "Stage III (vs I)",
           Variable == "stage IV"   ~ "Stage IV (vs I)",
           Variable == "subtype LumB" ~ "Luminal B (vs Luminal A)",
           Variable == "subtype HER2" ~ "HER2-enriched (vs Luminal A)",
           Variable == "subtype TN"   ~ "Triple Neg. (vs Luminal A)",
           Variable == "age"        ~ "Age (per year)",
           Variable == "lymph nodes" ~ "Lymph Nodes",
           TRUE ~ Variable
         ),
         significance = case_when(
           pvalue < 0.001 ~ "***",
           pvalue < 0.01  ~ "**",
           pvalue < 0.05  ~ "*",
           TRUE ~ ""
         ))

p_cox <- ggplot(cox_df, aes(x = HR, y = reorder(Variable, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(size = 3, color = npg10[3]) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2,
                 linewidth = 1, color = npg10[3]) +
  geom_text(aes(label = sprintf("HR=%.2f %s", HR, significance)),
            hjust = -0.3, size = 3.2, family = "songti") +
  scale_x_log10(limits = c(0.1, 60)) +
  labs(title = "BRCA: Multivariate Cox Regression",
       subtitle = sprintf("C-index = 0.771  |  LR p = %s", "8.62e-12"),
       x = "Hazard Ratio (95% CI, log scale)", y = "") +
  theme_pub

save_png(p_cox, "fig9_cox_forest.png", 9, 6)

# ---- 11. miRNA-mRNA correlation heatmap ----
cat("\nFig 10: miRNA-mRNA network...\n")

mirna_cor <- read.csv("D:/Users/Desktop/R_Work/results/tables/mirna_target_correlations.csv")
if (nrow(mirna_cor) > 10) {
  top_pairs <- mirna_cor %>% arrange(desc(abs(correlation))) %>% head(40)

  p_mirna <- ggplot(top_pairs, aes(x = mRNA, y = miRNA, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = npg10[5], mid = "white", high = npg10[3],
                         midpoint = 0, name = "r") +
    labs(title = "miRNA-mRNA Regulatory Pairs (|r|>0.4)",
         x = "mRNA", y = "miRNA") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
          axis.text.y = element_text(size = 5))

  save_png(p_mirna, "fig10_mirna_mrna_corr.png", 12, 8)
}

# ---- 12. Six-panel summary figure ----
cat("\nFig 11: Summary figure...\n")

# Build 6 panels
# P1: Subtype distribution
p_s1 <- clinical %>%
  filter(!is.na(molecular_subtype), molecular_subtype != "") %>%
  count(molecular_subtype) %>%
  ggplot(aes(x = reorder(molecular_subtype, n), y = n, fill = molecular_subtype)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  scale_fill_manual(values = subtype_cols, guide = "none") +
  coord_flip() +
  labs(title = "Subtype Distribution", x = "", y = "") +
  theme_pub + theme(axis.text.y = element_text(size = 7))

# P2: Stage distribution
p_s2 <- clinical %>%
  filter(!is.na(stage_simple), stage_simple != "") %>%
  count(stage_simple) %>%
  ggplot(aes(x = stage_simple, y = n, fill = stage_simple)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Stage Distribution", x = "", y = "") +
  theme_pub

# P3: Volcano (compact)
p_s3 <- p_volcano + theme(plot.title = element_text(size = 9),
                           plot.subtitle = element_text(size = 7))

# P4: Model comparison (compact)
p_s4 <- p_model + theme(plot.title = element_text(size = 9),
                          plot.subtitle = element_text(size = 7))

# P5: KM (compact — same as p_km$plot)
p_s5 <- p_km$plot + theme(plot.title = element_text(size = 9))

# P6: Module sizes (compact)
p_s6 <- p_mod + theme(plot.title = element_text(size = 9),
                        plot.subtitle = element_text(size = 7))

png(file.path(OUT_DIR, "fig11_summary.png"), width = 20, height = 12,
    units = "in", res = 300, bg = "white")
grid.arrange(p_s3, p_s1, p_s4,
             p_s5, p_s2, p_s6,
             ncol = 3, nrow = 2,
             top = textGrob("BRCA Multi-Omics Data Mining — Summary",
                            gp = gpar(fontsize = 18, fontface = "bold", fontfamily = "songti")))
dev.off()
cat("  Saved: fig11_summary.png\n")

# ---- Done ----
cat(sprintf("\n========== %d figures saved to %s ==========\n",
            length(list.files(OUT_DIR, pattern = "\\.png$")), OUT_DIR))
