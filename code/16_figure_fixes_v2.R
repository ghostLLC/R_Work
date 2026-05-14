# ===========================================================================
# 16v2: Figure Fixes — targeted by user feedback
# Fig3: KM risk table ugly + legend covers curve
# Fig4: white space + meaningless color blocks at bottom
# Fig5: dendrogram lines too thick/cluttered
# Fig6-10: font sizes too small
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
  library(maftools)
})

showtext_auto()
font_add("songti", regular = "C:/Windows/Fonts/simsun.ttc", bold = "C:/Windows/Fonts/simhei.ttf")

npg10 <- pal_npg("nrc")(10)
OUT <- "D:/Users/Desktop/R_Work/results/figures_pub"

# Load data
brca      <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_tumor_processed.rds")
tpm_mat   <- brca$tpm
clinical  <- brca$clinical
counts_tumor  <- brca$counts
counts_normal <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_normal_counts.rds")
deg_table     <- read.csv("D:/Users/Desktop/R_Work/results/tables/brca_degs_deseq2.csv")
cox_table     <- read.csv("D:/Users/Desktop/R_Work/results/tables/cox_regression.csv")
class_metrics <- read.csv("D:/Users/Desktop/R_Work/results/tables/classification_metrics.csv")

expr_pid <- substr(colnames(tpm_mat), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_pid, clinical$pid_short)

# ===================================================================
# Fig3: KM — NO risk table, legend moved to bottom-right corner
# ===================================================================
cat("\n=== Fig3: KM (no risk table, legend fixed) ===\n")

clinical_surv <- clinical %>%
  filter(!is.na(stage_simple), stage_simple != "", os_time > 0) %>%
  mutate(os_years = os_time / 365.25)

fit_stage <- survfit(Surv(os_years, vital_status) ~ stage_simple, data = clinical_surv)
stage_cols <- c("Stage I" = npg10[1], "Stage II" = npg10[2],
                "Stage III" = npg10[3], "Stage IV" = npg10[4])

# Use ggsurvplot WITHOUT risk table, legend at bottom-right
p_km <- ggsurvplot(
  fit_stage, data = clinical_surv,
  pval = TRUE, pval.size = 5, pval.coord = c(1, 0.15),
  palette = stage_cols,
  legend = c(0.72, 0.22),         # bottom-right, below the curves
  legend.title = "Stage",
  legend.labs = c("Stage I", "Stage II", "Stage III", "Stage IV"),
  xlab = "Time (Years)", ylab = "Overall Survival",
  risk.table = "nrisk_cumevents",   # compact: n at risk (cumulative events)
  risk.table.height = 0.22,
  risk.table.y.text = FALSE,
  tables.theme = theme_cleantable(base_size = 11),
  ggtheme = theme_classic(base_size = 15, base_family = "songti"),
  title = "BRCA: Overall Survival by Pathologic Stage",
  break.time.by = 5,
  surv.median.line = "hv"
)

png(file.path(OUT, "fig3_km_stage.png"), width = 9, height = 7.5,
    units = "in", res = 300, bg = "white")
print(p_km)
dev.off()
cat("  Saved: fig3_km_stage.png\n")

# ===================================================================
# Fig4: Oncoplot — remove colBar, drawRowBar only, tighter layout
# ===================================================================
cat("\n=== Fig4: Oncoplot (tight, no colBar, top 12) ===\n")

maf_file <- "D:/Users/Desktop/R_Work/data/public_data/brca_mutations_temp.maf"
if (file.exists(maf_file)) {
  maf <- read.maf(maf = maf_file, verbose = FALSE)
  png(file.path(OUT, "fig4_oncoplot.png"), width = 11, height = 5.5,
      units = "in", res = 300, bg = "white")
  oncoplot(
    maf = maf,
    top = 12,
    legendFontSize = 1.2,      # relative size
    titleText = "BRCA: Top 12 Mutated Genes",
    drawRowBar = TRUE,
    drawColBar = FALSE,        # remove bottom color blocks
    showTumorSampleBarcodes = FALSE,
    gene_mar = 6
  )
  dev.off()
  cat("  Saved: fig4_oncoplot.png (11x5.5)\n")
}

# ===================================================================
# Fig5: Heatmap — thinner dendrograms, top 25 genes, fewer samples
# ===================================================================
cat("\n=== Fig5: DEG Heatmap (thin dendro, top 25, sampled) ===\n")

common_genes <- intersect(rownames(counts_tumor), rownames(counts_normal))
top25_genes <- deg_table %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(25) %>%
  pull(gene)
top25_genes <- intersect(top25_genes, common_genes)

if (length(top25_genes) >= 10) {
  # Subsample for cleaner look: 50 tumor + all normal
  set.seed(42)
  tumor_idx <- sample(ncol(counts_tumor), min(80, ncol(counts_tumor)))
  normal_idx <- 1:ncol(counts_normal)

  vsd_mat <- log2(cbind(
    counts_tumor[top25_genes, tumor_idx],
    counts_normal[top25_genes, normal_idx]
  ) + 1)

  annot_col <- data.frame(
    Condition = c(rep("Tumor", length(tumor_idx)), rep("Normal", length(normal_idx))),
    row.names = colnames(vsd_mat)
  )
  rn <- gsub("_ENSG\\d+$", "", rownames(vsd_mat))

  png(file.path(OUT, "fig5_deg_heatmap.png"), width = 9, height = 7,
      units = "in", res = 300, bg = "white")
  pheatmap(vsd_mat,
           scale = "row",
           labels_row = rn,
           annotation_col = annot_col,
           annotation_colors = list(Condition = c(Tumor = npg10[3], Normal = npg10[5])),
           show_colnames = FALSE,
           fontsize_row = 9,
           fontsize = 12,
           treeheight_row = 15,     # shorter dendrogram
           treeheight_col = 20,
           color = colorRampPalette(c(npg10[5], "gray95", npg10[3]))(80),
           border_color = NA,
           main = "BRCA: Top 25 DEGs (Tumor vs Normal)")
  dev.off()
  cat("  Saved: fig5_deg_heatmap.png (9x7, top25)\n")
}

# ===================================================================
# Fig6: Model comparison — BIG font, simpler design
# ===================================================================
cat("\n=== Fig6: Model comparison (BIG font) ===\n")

model_detail <- data.frame(
  Model = rep(c("LASSO", "Random Forest", "XGBoost"), each = 3),
  Metric = rep(c("Accuracy", "Kappa", "F1"), 3),
  Value = c(0.894, 0.528, 0.612, 0.867, 0.365, 0.706, 0.358, 0.002, 0.232)
)
model_detail$Model <- factor(model_detail$Model, levels = c("LASSO", "Random Forest", "XGBoost"))

p6 <- ggplot(model_detail, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.65, alpha = 0.88) +
  geom_text(aes(label = sprintf("%.3f", Value)),
            position = position_dodge(width = 0.65),
            vjust = -0.6, size = 5.5, family = "songti", fontface = "bold") +
  scale_fill_manual(values = npg10[1:3], name = "") +
  labs(title = "BRCA Subtype Classification Performance",
       subtitle = "3-Class: Luminal A | HER2-enriched | Triple Negative",
       x = "", y = "Score") +
  ylim(0, 1.08) +
  theme_classic(base_size = 16, base_family = "songti") +
  theme(
    plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title    = element_text(size = 15),
    axis.text     = element_text(size = 14, color = "black"),
    legend.text   = element_text(size = 13),
    legend.position = "bottom",
    axis.line     = element_line(linewidth = 0.6),
    axis.ticks    = element_line(linewidth = 0.6)
  )

ggsave(file.path(OUT, "fig6_model_comparison.png"), p6,
       width = 9, height = 6.5, dpi = 300, bg = "white")
cat("  Saved: fig6_model_comparison.png\n")

# ===================================================================
# Fig7: WGCNA SFT — BIG font
# ===================================================================
cat("\n=== Fig7: WGCNA SFT (BIG font) ===\n")

sft_data <- data.frame(
  Power = 1:20,
  SFT_R2 = c(0.253, 0.148, 0.020, 0.122, 0.442, 0.671, 0.817, 0.887,
             0.918, 0.942, 0.948, 0.958, 0.964, 0.962, 0.966, 0.964,
             0.971, 0.976, 0.976, 0.977),
  Mean_Conn = c(2590, 1390, 766, 436, 256, 155, 97.2, 62.9,
                42.1, 29.0, 20.7, 15.1, 11.4, 8.78, 6.92, 5.55,
                4.54, 3.76, 3.16, 2.68)
)

thm_big <- theme_classic(base_size = 16, base_family = "songti") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15))

p7a <- ggplot(sft_data, aes(Power, SFT_R2)) +
  geom_point(size = 3, color = npg10[3]) +
  geom_line(color = npg10[3], linewidth = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_vline(xintercept = 8, linetype = "dotted", color = npg10[4], linewidth = 0.7) +
  annotate("text", x = 12, y = 0.83, label = "R^2 == 0.887", size = 6,
           family = "songti", parse = TRUE) +
  labs(x = "Soft Threshold (Power)", y = expression(Scale~Free~Topology~R^2)) +
  thm_big

p7b <- ggplot(sft_data, aes(Power, Mean_Conn)) +
  geom_point(size = 3, color = npg10[1]) +
  geom_line(color = npg10[1], linewidth = 1) +
  labs(x = "Soft Threshold (Power)", y = "Mean Connectivity") +
  thm_big

png(file.path(OUT, "fig7_wgcna_sft.png"), width = 14, height = 5.5,
    units = "in", res = 300, bg = "white")
grid.arrange(p7a, p7b, ncol = 2,
             top = textGrob("WGCNA: Soft-Threshold Power Selection",
                            gp = gpar(fontsize = 20, fontface = "bold", fontfamily = "songti")))
dev.off()
cat("  Saved: fig7_wgcna_sft.png\n")

# ===================================================================
# Fig8: WGCNA modules — BIG font
# ===================================================================
cat("\n=== Fig8: WGCNA modules (BIG font) ===\n")

module_sizes <- data.frame(
  Module = c("turquoise","blue","brown","yellow","green","red","black","pink"),
  Size   = c(584, 553, 336, 260, 224, 119, 63, 44)
)
mod_colors <- c("turquoise"="#40E0D0","blue"="#377EB8","brown"="#A52A2A",
                "yellow"="#FFD700","green"="#4DAF4A","red"="#E41A1C",
                "black"="#333333","pink"="#F781BF")

p8 <- ggplot(module_sizes, aes(x = reorder(Module, Size), y = Size, fill = Module)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.9) +
  geom_text(aes(label = Size), hjust = -0.3, size = 7, family = "songti", fontface = "bold") +
  scale_fill_manual(values = mod_colors, guide = "none") +
  coord_flip(ylim = c(0, max(module_sizes$Size) * 1.18)) +
  labs(title = "WGCNA Co-expression Modules",
       subtitle = "Soft Power = 8, excl. grey module",
       x = "", y = "Number of Genes") +
  theme_classic(base_size = 17, base_family = "songti") +
  theme(
    plot.title    = element_text(size = 19, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.text.y   = element_text(size = 15, color = "black"),
    axis.text.x   = element_text(size = 14, color = "black"),
    axis.title    = element_text(size = 15),
    axis.line     = element_line(linewidth = 0.6),
    axis.ticks    = element_line(linewidth = 0.6)
  )

ggsave(file.path(OUT, "fig8_wgcna_modules.png"), p8,
       width = 9, height = 6.5, dpi = 300, bg = "white")
cat("  Saved: fig8_wgcna_modules.png\n")

# ===================================================================
# Fig9: Cox forest — BIG font, big dots, big error bars
# ===================================================================
cat("\n=== Fig9: Cox forest (BIG font) ===\n")

cox_df <- cox_table %>%
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
    sig = case_when(
      pvalue < 0.001 ~ "***", pvalue < 0.01 ~ "**", pvalue < 0.05 ~ "*", TRUE ~ ""
    ),
    label_x = pmax(CI_upper * 1.8, HR * 3)
  )

p9 <- ggplot(cox_df, aes(x = HR, y = reorder(Variable, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.7) +
  geom_point(size = 5, color = npg10[3]) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.3, linewidth = 2, color = npg10[3]) +
  geom_text(aes(x = label_x, label = sprintf("HR=%.1f %s", HR, sig)),
            hjust = 0, size = 6, family = "songti", fontface = "bold") +
  scale_x_log10(limits = c(0.03, 90)) +
  labs(title = "BRCA: Multivariate Cox Regression",
       subtitle = expression(C-index == 0.771~~~~LR~italic(p) == 8.62 %*% 10^-12),
       x = "Hazard Ratio (95% CI, log scale)", y = "") +
  theme_classic(base_size = 17, base_family = "songti") +
  theme(
    plot.title    = element_text(size = 19, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.text.y   = element_text(size = 15, color = "black", face = "bold"),
    axis.text.x   = element_text(size = 14, color = "black"),
    axis.title    = element_text(size = 15),
    axis.line     = element_line(linewidth = 0.6),
    axis.ticks    = element_line(linewidth = 0.6)
  )

ggsave(file.path(OUT, "fig9_cox_forest.png"), p9,
       width = 11, height = 7.5, dpi = 300, bg = "white")
cat("  Saved: fig9_cox_forest.png\n")

# ===================================================================
# Fig10: miRNA-mRNA — BIG font, clean heatmap
# ===================================================================
cat("\n=== Fig10: miRNA-mRNA (BIG font, cleaner) ===\n")

mirna_cor <- read.csv("D:/Users/Desktop/R_Work/results/tables/mirna_target_correlations.csv")
if (nrow(mirna_cor) > 10) {
  top_pairs <- mirna_cor %>%
    arrange(desc(abs(correlation))) %>%
    head(20) %>%
    mutate(
      miRNA_short = gsub("^read_count_|^hsa-", "", miRNA),
      mRNA_short  = gsub("_ENSG\\d+$", "", mRNA)
    )

  p10 <- ggplot(top_pairs, aes(x = mRNA_short, y = miRNA_short)) +
    geom_tile(aes(fill = correlation), color = "white", linewidth = 0.8) +
    scale_fill_gradient2(low = npg10[5], mid = "white", high = npg10[3],
                         midpoint = 0, name = "r") +
    geom_text(aes(label = sprintf("%.2f", correlation)),
              size = 5, family = "songti", fontface = "bold") +
    labs(title = "miRNA-mRNA Regulatory Pairs (Top 20)",
         subtitle = paste0("|r| > 0.4, Total pairs: ",
                           nrow(mirna_cor %>% filter(abs(correlation) > 0.4))),
         x = "mRNA", y = "miRNA") +
    theme_minimal(base_size = 16, base_family = "songti") +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      plot.title    = element_text(size = 19, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      legend.title   = element_text(size = 13),
      legend.text    = element_text(size = 12),
      panel.grid = element_blank()
    )

  ggsave(file.path(OUT, "fig10_mirna_mrna_corr.png"), p10,
         width = 11, height = 8.5, dpi = 300, bg = "white")
  cat("  Saved: fig10_mirna_mrna_corr.png\n")
}

cat("\n========== All v2 fixes complete ==========\n")
