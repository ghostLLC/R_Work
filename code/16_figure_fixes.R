# ===========================================================================
# 16: Figure Quality Fixes — targeted improvements
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
font_add("heiti",  regular = "C:/Windows/Fonts/simhei.ttf")

npg10 <- pal_npg("nrc")(10)
OUT    <- "D:/Users/Desktop/R_Work/results/figures_pub"

theme_pub <- theme_classic(base_size = 14, base_family = "songti") +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title    = element_text(size = 13),
    axis.text     = element_text(size = 11, color = "black"),
    legend.title  = element_text(size = 11),
    legend.text   = element_text(size = 10),
    plot.margin   = margin(12, 12, 12, 12),
    panel.grid    = element_blank(),
    axis.line     = element_line(color = "black", linewidth = 0.5),
    axis.ticks    = element_line(color = "black", linewidth = 0.5)
  )

save_png <- function(p, f, w = 8, h = 6) {
  ggsave(file.path(OUT, f), plot = p, width = w, height = h, dpi = 300, units = "in", bg = "white")
  cat(sprintf("  Saved: %s (%dx%d)\n", f, w, h))
}

# ---- Load data ----
brca      <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_tumor_processed.rds")
tpm_mat   <- brca$tpm
clinical  <- brca$clinical
cox_table <- read.csv("D:/Users/Desktop/R_Work/results/tables/cox_regression.csv")
class_metrics <- read.csv("D:/Users/Desktop/R_Work/results/tables/classification_metrics.csv")
deg_table <- read.csv("D:/Users/Desktop/R_Work/results/tables/brca_degs_deseq2.csv")
counts_tumor <- brca$counts
counts_normal <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_normal_counts.rds")

expr_pid <- substr(colnames(tpm_mat), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_pid, clinical$pid_short)

# ===================================================================
# Fig3 FIX: KM curve with cleaner risk table
# ===================================================================
cat("\n=== Fig3: KM curve (fixed risk table) ===\n")

clinical_surv <- clinical %>%
  filter(!is.na(stage_simple), stage_simple != "", os_time > 0) %>%
  mutate(os_years = os_time / 365.25)

# Use all 4 stages directly (more informative than 2-stage grouping)
fit_stage <- survfit(Surv(os_years, vital_status) ~ stage_simple, data = clinical_surv)

stage_cols <- c("Stage I" = npg10[1], "Stage II" = npg10[2],
                "Stage III" = npg10[3], "Stage IV" = npg10[4])

# Use ggsurvplot with custom risk table theme
p_km <- ggsurvplot(
  fit_stage,
  data = clinical_surv,
  pval = TRUE, pval.size = 4.5, pval.coord = c(0.5, 0.12),
  palette = stage_cols,
  legend = c(0.7, 0.85),
  legend.title = "Stage",
  legend.labs = c("Stage I", "Stage II", "Stage III", "Stage IV"),
  xlab = "Time (Years)", ylab = "Overall Survival",
  risk.table = TRUE,
  risk.table.height = 0.28,
  risk.table.y.text = TRUE,
  risk.table.y.text.col = FALSE,
  tables.theme = theme_cleantable(base_size = 11),
  ggtheme = theme_classic(base_size = 14, base_family = "songti"),
  title = "BRCA: Overall Survival by Pathologic Stage (I-IV)",
  break.time.by = 5,
  surv.median.line = "hv"
)

png(file.path(OUT, "fig3_km_stage.png"), width = 10, height = 8,
    units = "in", res = 300, bg = "white")
print(p_km)
dev.off()
cat("  Saved: fig3_km_stage.png (10x8)\n")

# ===================================================================
# Fig4 FIX: Oncoplot — top 15 genes, smaller sample subset
# ===================================================================
cat("\n=== Fig4: Oncoplot (cleaner, top 15) ===\n")

maf_file <- "D:/Users/Desktop/R_Work/data/public_data/brca_mutations_temp.maf"
if (file.exists(maf_file)) {
  maf <- read.maf(maf = maf_file, verbose = FALSE)

  png(file.path(OUT, "fig4_oncoplot.png"), width = 12, height = 7,
      units = "in", res = 300, bg = "white")
  oncoplot(
    maf = maf,
    top = 15,
    legendFontSize = 12,
    titleText = "BRCA: Top 15 Mutated Genes",
    drawRowBar = TRUE,
    drawColBar = FALSE
  )
  dev.off()
  cat("  Saved: fig4_oncoplot.png (12x7, top15)\n")
}

# ===================================================================
# Fig5 FIX: DEG Heatmap — top 30 genes, larger font
# ===================================================================
cat("\n=== Fig5: DEG Heatmap (cleaner, top 30) ===\n")

common_genes <- intersect(rownames(counts_tumor), rownames(counts_normal))
top30_genes <- deg_table %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(30) %>%
  pull(gene)
top30_genes <- intersect(top30_genes, common_genes)

if (length(top30_genes) >= 10) {
  vsd_mat <- log2(cbind(counts_tumor[top30_genes, ], counts_normal[top30_genes, ]) + 1)

  annot_col <- data.frame(
    Condition = c(rep("Tumor", ncol(counts_tumor)), rep("Normal", ncol(counts_normal))),
    row.names = colnames(vsd_mat)
  )

  rn <- gsub("_ENSG\\d+$", "", rownames(vsd_mat))

  png(file.path(OUT, "fig5_deg_heatmap.png"), width = 10, height = 8,
      units = "in", res = 300, bg = "white")
  pheatmap(vsd_mat,
           scale = "row",
           labels_row = rn,
           annotation_col = annot_col,
           annotation_colors = list(Condition = c(Tumor = npg10[3], Normal = npg10[5])),
           show_colnames = FALSE,
           fontsize_row = 8,
           fontsize = 11,
           color = colorRampPalette(c(npg10[5], "white", npg10[3]))(100),
           main = "BRCA: Top 30 DEGs (Tumor vs Normal)")
  dev.off()
  cat("  Saved: fig5_deg_heatmap.png (10x8, top30)\n")
}

# ===================================================================
# Fig6 FIX: Model comparison — richer (add F1, Kappa per class)
# ===================================================================
cat("\n=== Fig6: Model comparison (enriched) ===\n")

# Build richer comparison data
model_detail <- data.frame(
  Model = rep(c("LASSO", "Random Forest", "XGBoost"), each = 3),
  Metric = rep(c("Accuracy", "Kappa", "F1 (Macro)"), 3),
  Value = c(0.894, 0.528, 0.612,   # LASSO
            0.867, 0.365, 0.706,   # RF
            0.358, 0.002, 0.232),  # XGBoost
  stringsAsFactors = FALSE
)
model_detail$Model <- factor(model_detail$Model, levels = c("LASSO", "Random Forest", "XGBoost"))

p_model2 <- ggplot(model_detail, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.3f", Value)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.8, family = "songti") +
  scale_fill_manual(values = npg10[1:3], name = "Metric") +
  labs(title = "BRCA Molecular Subtype Classification",
       subtitle = "3 models compared on 3 metrics (3 classes)",
       x = "", y = "Score") +
  ylim(0, 1.05) +
  theme_pub + theme(legend.position = "bottom")

save_png(p_model2, "fig6_model_comparison.png", 9, 6)

# ===================================================================
# Fig8 FIX: WGCNA modules — larger font, better proportions
# ===================================================================
cat("\n=== Fig8: WGCNA modules (larger font) ===\n")

module_sizes <- data.frame(
  Module = c("turquoise","blue","brown","yellow","green","red","black","pink"),
  Size   = c(584, 553, 336, 260, 224, 119, 63, 44),
  Color  = c("turquoise","blue","brown","yellow","green","red","black","pink")
)
mod_colors <- c("turquoise"="#40E0D0","blue"="#377EB8","brown"="#A52A2A",
                "yellow"="#FFD700","green"="#4DAF4A","red"="#E41A1C",
                "black"="#333333","pink"="#F781BF")

p_mod2 <- ggplot(module_sizes, aes(x = reorder(Module, Size), y = Size, fill = Module)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.9) +
  geom_text(aes(label = Size), hjust = -0.3, size = 5, family = "songti") +
  scale_fill_manual(values = mod_colors, guide = "none") +
  coord_flip(ylim = c(0, max(module_sizes$Size) * 1.15)) +
  labs(title = "WGCNA Co-expression Modules (excl. grey)",
       subtitle = "Soft Power = 8 | 9 modules total",
       x = "", y = "Number of Genes") +
  theme_pub + theme(axis.text.y = element_text(size = 13))

save_png(p_mod2, "fig8_wgcna_modules.png", 8, 6)

# ===================================================================
# Fig9 FIX: Cox forest — bigger font, cleaner labels
# ===================================================================
cat("\n=== Fig9: Cox forest plot (larger font) ===\n")

cox_df <- cox_table %>%
  filter(!is.infinite(HR), HR < 100, HR > 0.01) %>%
  mutate(
    Variable = case_when(
      Variable == "stage_II"   ~ "Stage II (vs I)",
      Variable == "stage_III"  ~ "Stage III (vs I)",
      Variable == "stage_IV"   ~ "Stage IV (vs I)",
      Variable == "subtype_LumB" ~ "Luminal B",
      Variable == "subtype_HER2" ~ "HER2-enriched",
      Variable == "subtype_TN"   ~ "Triple Negative",
      Variable == "age"        ~ "Age (10yr)",
      Variable == "lymph_nodes" ~ "Lymph Nodes (+)",
      TRUE ~ Variable
    ),
    sig = case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01  ~ "**",
      pvalue < 0.05  ~ "*",
      TRUE ~ ""
    ),
    label_pos = pmax(HR, CI_upper) * 1.5
  )

p_cox2 <- ggplot(cox_df, aes(x = HR, y = reorder(Variable, HR))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_point(size = 4, color = npg10[3]) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.25, linewidth = 1.5, color = npg10[3]) +
  geom_text(aes(x = label_pos,
                label = sprintf("HR=%.2f %s", HR, sig)),
            hjust = 0, size = 4.5, family = "songti") +
  scale_x_log10(limits = c(0.05, 80)) +
  labs(title = "BRCA: Multivariate Cox Regression",
       subtitle = "C-index = 0.771 | LR p = 8.62e-12",
       x = "Hazard Ratio (95% CI, log scale)", y = "") +
  theme_pub + theme(axis.text.y = element_text(size = 13, face = "bold"))

save_png(p_cox2, "fig9_cox_forest.png", 10, 7)

# ===================================================================
# Fig10 FIX: miRNA-mRNA — simpler, larger text, fewer pairs
# ===================================================================
cat("\n=== Fig10: miRNA-mRNA network (simplified) ===\n")

mirna_cor <- read.csv("D:/Users/Desktop/R_Work/results/tables/mirna_target_correlations.csv")
if (nrow(mirna_cor) > 10) {
  # Take top 25 pairs for cleaner visualization
  top_pairs <- mirna_cor %>%
    arrange(desc(abs(correlation))) %>%
    head(25) %>%
    mutate(
      miRNA_short = gsub("^read_count_|^hsa-", "", miRNA),
      mRNA_short  = gsub("_ENSG\\d+$", "", mRNA)
    )

  p_mirna2 <- ggplot(top_pairs,
                     aes(x = mRNA_short, y = miRNA_short, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = npg10[5], mid = "white", high = npg10[3],
                         midpoint = 0, name = "r") +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3.5,
              family = "songti", color = "grey20") +
    labs(title = "miRNA-mRNA Regulatory Pairs (Top 25, |r|>0.4)",
         subtitle = paste0("Total pairs: ", nrow(mirna_cor %>% filter(abs(correlation) > 0.4))),
         x = "mRNA", y = "miRNA") +
    theme_minimal(base_size = 13, base_family = "songti") +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "right"
    )

  save_png(p_mirna2, "fig10_mirna_mrna_corr.png", 10, 8)
}

cat("\n========== All fig fixes complete ==========\n")
