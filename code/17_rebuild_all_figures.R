# ===========================================================================
# 17: Complete Figure Rebuild — All English, Large Fonts, 300 DPI
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE, bitmapType = "cairo")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggsci)
  library(ggrepel)
  library(pheatmap)
  library(survival)
  library(survminer)
  library(gridExtra)
  library(grid)
  library(maftools)
})

npg10 <- pal_npg("nrc")(10)
OUT   <- "D:/Users/Desktop/R_Work/results/figures_pub"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# === GLOBAL THEME: BIG fonts, English only, clean ===
theme_big <- theme_classic(base_size = 20) +
  theme(
    plot.title       = element_text(size = 24, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 16, hjust = 0.5, color = "grey30"),
    axis.title       = element_text(size = 20),
    axis.text        = element_text(size = 18, color = "black"),
    legend.title     = element_text(size = 18),
    legend.text      = element_text(size = 16),
    legend.position  = "right",
    plot.margin      = margin(15, 15, 15, 15),
    panel.grid       = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.7),
    axis.ticks       = element_line(color = "black", linewidth = 0.7)
  )

theme_big_nolegend <- theme_big + theme(legend.position = "none")

# === Load data ===
brca   <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_tumor_processed.rds")
tpm    <- brca$tpm
clin   <- brca$clinical
ct     <- brca$counts
cn     <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_normal_counts.rds")
deg    <- read.csv("D:/Users/Desktop/R_Work/results/tables/brca_degs_deseq2.csv")
cox_df <- read.csv("D:/Users/Desktop/R_Work/results/tables/cox_regression.csv")
cm     <- read.csv("D:/Users/Desktop/R_Work/results/tables/classification_metrics.csv")

pid <- substr(colnames(tpm), 9, 12)
clin$pid_s <- gsub("^TCGA-\\w{2}-", "", clin$patient_id)
midx <- match(pid, clin$pid_s)

# ===================================================================
# Fig1: Volcano
# ===================================================================
cat("Fig1: Volcano\n")

sig <- deg %>% filter(!is.na(padj), padj < 0.05, abs(log2FC) > 1)
top_lbl <- bind_rows(
  sig %>% filter(log2FC > 0) %>% arrange(padj) %>% head(10),
  sig %>% filter(log2FC < 0) %>% arrange(padj) %>% head(10)
)

vdf <- deg %>%
  mutate(
    lp  = pmin(-log10(pvalue), 50),
    lfc = pmax(pmin(log2FC, 6), -6),
    sig = case_when(padj < 0.05 & log2FC > 1  ~ "Up",
                    padj < 0.05 & log2FC < -1 ~ "Down", TRUE ~ "NS"),
    lbl = ifelse(gene %in% top_lbl$gene, gsub("_ENSG\\d+$", "", gene), "")
  )

p1 <- ggplot(vdf, aes(lfc, lp)) +
  geom_point(aes(color = sig), size = 1.2, alpha = 0.45) +
  scale_color_manual(values = c("Up" = npg10[3], "Down" = npg10[5], "NS" = "grey75")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_text_repel(aes(label = lbl), size = 4.5, max.overlaps = 30,
                  box.padding = 0.4, segment.size = 0.3) +
  labs(title = "BRCA: Tumor vs Normal",
       subtitle = sprintf("6,768 DEGs (%d Up | %d Down)", sum(vdf$sig == "Up"), sum(vdf$sig == "Down")),
       x = expression(log[2]~Fold~Change), y = expression(-log[10](italic(p)))) +
  annotate("text", x = 4.5, y = 3, label = sprintf("Up: %d", sum(vdf$sig == "Up")),
           color = npg10[3], size = 6, fontface = "bold") +
  annotate("text", x = -4.5, y = 3, label = sprintf("Down: %d", sum(vdf$sig == "Down")),
           color = npg10[5], size = 6, fontface = "bold") +
  theme_big_nolegend

ggsave(file.path(OUT, "fig1_volcano.png"), p1, width = 10, height = 8, dpi = 300, bg = "white")

# ===================================================================
# Fig2: PCA
# ===================================================================
cat("Fig2: PCA\n")

lt <- log2(tpm + 1)
gv <- apply(lt, 1, var)
t2k <- names(sort(gv, decreasing = TRUE))[1:2000]
pca <- prcomp(t(lt[t2k, ]), center = TRUE, scale. = TRUE)
pv  <- summary(pca)$importance[2, 1:2] * 100

pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                     Subtype = clin$molecular_subtype[midx]) %>%
  filter(!is.na(Subtype), Subtype != "")

cols_sub <- c("Luminal A" = npg10[1], "Luminal B" = npg10[2],
              "HER2-enriched" = npg10[3], "Triple Negative" = npg10[4])

p2 <- ggplot(pca_df, aes(PC1, PC2, color = Subtype)) +
  geom_point(size = 2.5, alpha = 0.6) +
  stat_ellipse(linewidth = 1, show.legend = FALSE) +
  scale_color_manual(values = cols_sub, name = "Molecular Subtype") +
  labs(title = "BRCA PCA: Top 2,000 Variable Genes",
       x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2])) +
  theme_big

ggsave(file.path(OUT, "fig2_pca_subtype.png"), p2, width = 10, height = 8, dpi = 300, bg = "white")

# ===================================================================
# Fig3: KM — pure ggplot2, English, BIG
# ===================================================================
cat("Fig3: KM\n")

cs <- clin %>%
  filter(!is.na(stage_simple), stage_simple != "", os_time > 0) %>%
  mutate(os_years = os_time / 365.25)

fit <- survfit(Surv(os_years, vital_status) ~ stage_simple, data = cs)
lr  <- survdiff(Surv(os_years, vital_status) ~ stage_simple, data = cs)
pv_lr <- round(1 - pchisq(lr$chisq, df = length(lr$n) - 1), 6)

km_df <- data.frame(
  time  = fit$time, surv  = fit$surv,
  upper = fit$upper, lower = fit$lower,
  strata = rep(names(fit$strata), fit$strata)
) %>% mutate(strata = gsub("stage_simple=", "", strata))

stages_uniq <- unique(km_df$strata)
km_df <- bind_rows(
  data.frame(time = 0, surv = 1, upper = 1, lower = 1, strata = stages_uniq),
  km_df
)

sc <- c("Stage I" = npg10[1], "Stage II" = npg10[2],
        "Stage III" = npg10[3], "Stage IV" = npg10[4])

p3a <- ggplot(km_df, aes(time, surv, color = strata, fill = strata)) +
  geom_step(linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.06, show.legend = FALSE) +
  scale_color_manual(values = sc, name = "Stage") +
  scale_fill_manual(values = sc, guide = "none") +
  scale_x_continuous(breaks = seq(0, 25, 5), limits = c(0, 25), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.02), expand = c(0, 0),
                     labels = scales::percent_format()) +
  labs(title = "BRCA: Overall Survival by Stage",
       subtitle = sprintf("Log-rank p = %.4f  |  87 events / 1,105 patients", pv_lr),
       x = "Time (Years)", y = "Overall Survival") +
  theme_classic(base_size = 22) +
  theme(
    plot.title       = element_text(size = 26, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 16, hjust = 0.5, color = "grey30"),
    axis.title       = element_text(size = 22),
    axis.text        = element_text(size = 20, color = "black"),
    legend.title     = element_text(size = 20),
    legend.text      = element_text(size = 18),
    legend.position  = c(0.75, 0.22),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
    axis.line        = element_line(linewidth = 0.7),
    axis.ticks       = element_line(linewidth = 0.7)
  )

# Number at risk table
tcut <- c(0, 5, 10, 15, 20)
nr <- data.frame()
for (s in stages_uniq) {
  sub <- cs %>% filter(stage_simple == s) %>% mutate(os_years = os_time / 365.25)
  for (t in tcut) {
    n_at <- sum(sub$os_years >= t)
    nr <- rbind(nr, data.frame(Time = t, Stage = s, N = n_at))
  }
}
nr$Time <- factor(paste0(nr$Time, "yr"), levels = paste0(tcut, "yr"))

p3b <- ggplot(nr, aes(Time, Stage, label = N)) +
  geom_text(size = 6, fontface = "bold") +
  labs(title = "Number at Risk") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title  = element_text(size = 14, hjust = 0, face = "bold"),
    panel.grid  = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 16, face = "bold"),
    axis.title  = element_blank()
  )

png(file.path(OUT, "fig3_km_stage.png"), width = 11, height = 10, units = "in", res = 300, bg = "white")
grid.arrange(p3a, p3b + theme(plot.margin = margin(t = -5, l = 50)),
             ncol = 1, heights = c(4.2, 1))
dev.off()
cat("  Saved\n")

# ===================================================================
# Fig4: Oncoplot — top 12, big legend
# ===================================================================
cat("Fig4: Oncoplot\n")

mf <- "D:/Users/Desktop/R_Work/data/public_data/brca_mutations_temp.maf"
if (file.exists(mf)) {
  maf <- read.maf(maf = mf, verbose = FALSE)
  png(file.path(OUT, "fig4_oncoplot.png"), width = 12, height = 6,
      units = "in", res = 300, bg = "white")
  oncoplot(maf = maf, top = 12, legendFontSize = 1.5,
           titleText = "BRCA: Top 12 Mutated Genes",
           drawRowBar = TRUE, drawColBar = FALSE,
           showTumorSampleBarcodes = FALSE, gene_mar = 8)
  dev.off()
  cat("  Saved\n")
}

# ===================================================================
# Fig5: DEG Heatmap — top 25, thin dendro
# ===================================================================
cat("Fig5: Heatmap\n")

cg <- intersect(rownames(ct), rownames(cn))
t25 <- deg %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(25) %>% pull(gene)
t25 <- intersect(t25, cg)

if (length(t25) >= 10) {
  set.seed(42)
  ti <- sample(ncol(ct), min(80, ncol(ct)))
  ni <- 1:ncol(cn)

  mat <- log2(cbind(ct[t25, ti], cn[t25, ni]) + 1)
  ann <- data.frame(
    Condition = c(rep("Tumor", length(ti)), rep("Normal", length(ni))),
    row.names = colnames(mat)
  )
  rn <- gsub("_ENSG\\d+$", "", rownames(mat))

  png(file.path(OUT, "fig5_deg_heatmap.png"), width = 10, height = 8,
      units = "in", res = 300, bg = "white")
  pheatmap(mat, scale = "row", labels_row = rn,
           annotation_col = ann,
           annotation_colors = list(Condition = c(Tumor = npg10[3], Normal = npg10[5])),
           show_colnames = FALSE, fontsize_row = 10, fontsize = 14,
           treeheight_row = 15, treeheight_col = 20, border_color = NA,
           color = colorRampPalette(c(npg10[5], "gray96", npg10[3]))(80),
           main = "BRCA: Top 25 DEGs")
  dev.off()
  cat("  Saved\n")
}

# ===================================================================
# Fig6: Model comparison
# ===================================================================
cat("Fig6: Model comparison\n")

md <- data.frame(
  Model  = rep(c("LASSO", "Random Forest", "XGBoost"), each = 3),
  Metric = rep(c("Accuracy", "Kappa", "F1 (Macro)"), 3),
  Value  = c(0.894, 0.528, 0.612, 0.867, 0.365, 0.706, 0.358, 0.002, 0.232)
)
md$Model <- factor(md$Model, levels = c("LASSO", "Random Forest", "XGBoost"))

p6 <- ggplot(md, aes(Model, Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.65, alpha = 0.88) +
  geom_text(aes(label = sprintf("%.3f", Value)),
            position = position_dodge(width = 0.65),
            vjust = -0.6, size = 7, fontface = "bold") +
  scale_fill_manual(values = npg10[1:3], name = "") +
  labs(title = "BRCA Subtype Classification",
       subtitle = "3-Class: Luminal A | HER2-enriched | Triple Negative",
       x = "", y = "Score") +
  ylim(0, 1.1) +
  theme_classic(base_size = 22) +
  theme(
    plot.title       = element_text(size = 26, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 16, hjust = 0.5, color = "grey30"),
    axis.text        = element_text(size = 20, color = "black"),
    axis.title       = element_text(size = 20),
    legend.text      = element_text(size = 18),
    legend.position  = "bottom",
    axis.line        = element_line(linewidth = 0.7),
    axis.ticks       = element_line(linewidth = 0.7)
  )

ggsave(file.path(OUT, "fig6_model_comparison.png"), p6, width = 10, height = 8, dpi = 300, bg = "white")

# ===================================================================
# Fig7: WGCNA SFT
# ===================================================================
cat("Fig7: WGCNA SFT\n")

sft <- data.frame(
  Power = 1:20,
  SFT_R2 = c(0.253, 0.148, 0.020, 0.122, 0.442, 0.671, 0.817, 0.887,
             0.918, 0.942, 0.948, 0.958, 0.964, 0.962, 0.966, 0.964,
             0.971, 0.976, 0.976, 0.977),
  Mean_Conn = c(2590, 1390, 766, 436, 256, 155, 97.2, 62.9,
                42.1, 29.0, 20.7, 15.1, 11.4, 8.78, 6.92, 5.55,
                4.54, 3.76, 3.16, 2.68)
)

th7 <- theme_classic(base_size = 22) +
  theme(
    plot.title  = element_text(size = 22, face = "bold", hjust = 0.5),
    axis.text   = element_text(size = 20, color = "black"),
    axis.title  = element_text(size = 20),
    axis.line   = element_line(linewidth = 0.7),
    axis.ticks  = element_line(linewidth = 0.7)
  )

p7a <- ggplot(sft, aes(Power, SFT_R2)) +
  geom_point(size = 4, color = npg10[3]) +
  geom_line(color = npg10[3], linewidth = 1.2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_vline(xintercept = 8, linetype = "dotted", color = npg10[4], linewidth = 0.8) +
  annotate("text", x = 13, y = 0.83, label = expression(R^2 == 0.887), size = 8, parse = TRUE) +
  labs(x = "Soft Threshold (Power)", y = expression(Scale~Free~Topology~R^2)) + th7

p7b <- ggplot(sft, aes(Power, Mean_Conn)) +
  geom_point(size = 4, color = npg10[1]) +
  geom_line(color = npg10[1], linewidth = 1.2) +
  labs(x = "Soft Threshold (Power)", y = "Mean Connectivity") + th7

png(file.path(OUT, "fig7_wgcna_sft.png"), width = 16, height = 6.5,
    units = "in", res = 300, bg = "white")
grid.arrange(p7a, p7b, ncol = 2,
             top = textGrob("WGCNA: Soft-Threshold Power Selection",
                            gp = gpar(fontsize = 26, fontface = "bold")))
dev.off()
cat("  Saved\n")

# ===================================================================
# Fig8: WGCNA modules
# ===================================================================
cat("Fig8: WGCNA modules\n")

ms <- data.frame(
  Module = c("turquoise","blue","brown","yellow","green","red","black","pink"),
  Size   = c(584, 553, 336, 260, 224, 119, 63, 44)
)
mc <- c("turquoise"="#40E0D0","blue"="#377EB8","brown"="#A52A2A",
        "yellow"="#FFD700","green"="#4DAF4A","red"="#E41A1C",
        "black"="#333333","pink"="#F781BF")

p8 <- ggplot(ms, aes(x = reorder(Module, Size), y = Size, fill = Module)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.9) +
  geom_text(aes(label = Size), hjust = -0.3, size = 10, fontface = "bold") +
  scale_fill_manual(values = mc, guide = "none") +
  coord_flip(ylim = c(0, max(ms$Size) * 1.18)) +
  labs(title = "WGCNA Co-expression Modules",
       subtitle = "Soft Power = 8  |  excl. grey module",
       x = "", y = "Number of Genes") +
  theme_classic(base_size = 24) +
  theme(
    plot.title    = element_text(size = 28, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 18, hjust = 0.5, color = "grey30"),
    axis.text.y   = element_text(size = 22, color = "black"),
    axis.text.x   = element_text(size = 20, color = "black"),
    axis.title    = element_text(size = 22),
    axis.line     = element_line(linewidth = 0.7),
    axis.ticks    = element_line(linewidth = 0.7)
  )

ggsave(file.path(OUT, "fig8_wgcna_modules.png"), p8, width = 10, height = 8, dpi = 300, bg = "white")

# ===================================================================
# Fig9: Cox forest
# ===================================================================
cat("Fig9: Cox forest\n")

cdf <- cox_df %>%
  filter(!is.infinite(HR), HR < 100, HR > 0.01) %>%
  mutate(
    Variable = case_when(
      Variable == "stage_II"   ~ "Stage II (vs I)",
      Variable == "stage_III"  ~ "Stage III (vs I)",
      Variable == "stage_IV"   ~ "Stage IV (vs I)",
      Variable == "subtype_LumB"  ~ "Luminal B",
      Variable == "subtype_HER2"  ~ "HER2-enriched",
      Variable == "subtype_TN"    ~ "Triple Negative",
      Variable == "age"           ~ "Age (per 10yr)",
      Variable == "lymph_nodes"   ~ "Lymph Nodes",
      TRUE ~ Variable
    ),
    sig = case_when(pvalue < 0.001 ~ "***", pvalue < 0.01 ~ "**",
                    pvalue < 0.05 ~ "*", TRUE ~ ""),
    label_x = pmax(CI_upper * 1.6, HR * 2.5)
  )

p9 <- ggplot(cdf, aes(HR, reorder(Variable, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.8) +
  geom_point(size = 6, color = npg10[3]) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.35, linewidth = 2.5, color = npg10[3]) +
  geom_text(aes(x = label_x, label = sprintf("HR=%.1f %s", HR, sig)),
            hjust = 0, size = 8, fontface = "bold") +
  scale_x_log10(limits = c(0.02, 80)) +
  labs(title = "BRCA: Multivariate Cox Regression",
       subtitle = expression(C-index == 0.771~~~~LR~italic(p) == 8.62 %*% 10^-12),
       x = "Hazard Ratio (95% CI)", y = "") +
  theme_classic(base_size = 24) +
  theme(
    plot.title    = element_text(size = 28, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 18, hjust = 0.5, color = "grey30"),
    axis.text.y   = element_text(size = 22, color = "black", face = "bold"),
    axis.text.x   = element_text(size = 20, color = "black"),
    axis.title    = element_text(size = 22),
    axis.line     = element_line(linewidth = 0.7),
    axis.ticks    = element_line(linewidth = 0.7)
  )

ggsave(file.path(OUT, "fig9_cox_forest.png"), p9, width = 12, height = 9, dpi = 300, bg = "white")

# ===================================================================
# Fig10: miRNA-mRNA correlation
# ===================================================================
cat("Fig10: miRNA-mRNA\n")

mcor <- read.csv("D:/Users/Desktop/R_Work/results/tables/mirna_target_correlations.csv")
if (nrow(mcor) > 10) {
  tp <- mcor %>%
    arrange(desc(abs(correlation))) %>%
    head(20) %>%
    mutate(miRNA_s = gsub("^read_count_|^hsa-", "", miRNA),
           mRNA_s  = gsub("_ENSG\\d+$", "", mRNA))

  p10 <- ggplot(tp, aes(mRNA_s, miRNA_s)) +
    geom_tile(aes(fill = correlation), color = "white", linewidth = 1) +
    scale_fill_gradient2(low = npg10[5], mid = "white", high = npg10[3],
                         midpoint = 0, name = "r") +
    geom_text(aes(label = sprintf("%.2f", correlation)),
              size = 7, fontface = "bold") +
    labs(title = "miRNA-mRNA Regulatory Pairs (Top 20)",
         subtitle = paste0("|r| > 0.4  |  Total: ",
                           nrow(mcor %>% filter(abs(correlation) > 0.4)), " pairs"),
         x = "mRNA", y = "miRNA") +
    theme_minimal(base_size = 22) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 16, color = "black"),
      axis.text.y = element_text(size = 16, color = "black"),
      plot.title    = element_text(size = 28, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 18, hjust = 0.5, color = "grey30"),
      legend.title  = element_text(size = 18),
      legend.text   = element_text(size = 16),
      panel.grid    = element_blank()
    )

  ggsave(file.path(OUT, "fig10_mirna_mrna_corr.png"), p10,
         width = 13, height = 10, dpi = 300, bg = "white")
  cat("  Saved\n")
}

# ===================================================================
# Fig11: Summary 6-panel
# ===================================================================
cat("Fig11: Summary\n")

p_s1 <- clin %>%
  filter(!is.na(molecular_subtype), molecular_subtype != "") %>%
  count(molecular_subtype) %>%
  ggplot(aes(reorder(molecular_subtype, n), n, fill = molecular_subtype)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  scale_fill_manual(values = cols_sub, guide = "none") +
  coord_flip() +
  labs(title = "Subtype Distribution", x = "", y = "") +
  theme_big + theme(axis.text.y = element_text(size = 12))

p_s2 <- clin %>%
  filter(!is.na(stage_simple), stage_simple != "") %>%
  count(stage_simple) %>%
  ggplot(aes(stage_simple, n, fill = stage_simple)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Stage Distribution", x = "", y = "") +
  theme_big

p_s3 <- p1 + theme(plot.title = element_text(size = 14),
                    plot.subtitle = element_text(size = 10))

p_s4 <- p6 + theme(plot.title = element_text(size = 14),
                    plot.subtitle = element_text(size = 10))

p_s5 <- p3a + theme(plot.title = element_text(size = 14),
                     plot.subtitle = element_text(size = 10))

p_s6 <- p8 + theme(plot.title = element_text(size = 14),
                    plot.subtitle = element_text(size = 10))

png(file.path(OUT, "fig11_summary.png"), width = 22, height = 14,
    units = "in", res = 300, bg = "white")
grid.arrange(p_s3, p_s1, p_s4, p_s5, p_s2, p_s6,
             ncol = 3, nrow = 2,
             top = textGrob("BRCA Multi-Omics Data Mining — Summary",
                            gp = gpar(fontsize = 28, fontface = "bold")))
dev.off()
cat("  Saved\n")

cat("\n========== All 11 figures regenerated ==========\n")
