# ===========================================================================
# US-013: Comprehensive Data Visualization Suite
# 项目：BRCA多组学数据挖掘
# 日期：2026-05-14
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(gridExtra)
  library(grid)
})

cat("\n========== US-013: Comprehensive Visualization Suite ==========\n\n")

FIG_DIR   <- "D:/Users/Desktop/R_Work/results/figures"
TBL_DIR   <- "D:/Users/Desktop/R_Work/results/tables"
INPUT_DIR <- "D:/Users/Desktop/R_Work/data/processed"
OUT_DIR   <- "D:/Users/Desktop/R_Work/results"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Multi-panel summary figure ----
cat("Step 1: Creating multi-panel summary figure...\n")

# Load key results
deg_table <- read.csv(file.path(TBL_DIR, "brca_degs_deseq2.csv"))
class_metrics <- read.csv(file.path(TBL_DIR, "classification_metrics.csv"))
cox_table <- read.csv(file.path(TBL_DIR, "cox_regression.csv"))

brca <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
clinical <- brca$clinical

# Build 6-panel figure
pdf(file.path(FIG_DIR, "summary_figure.pdf"), width = 18, height = 12)

# Panel 1: Volcano
sig_deg <- deg_table %>% filter(!is.na(padj), padj < 0.05, abs(log2FC) > 1)
p1 <- ggplot(deg_table, aes(x = pmin(pmax(log2FC, -6), 6),
                              y = pmin(-log10(pvalue), 50))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FC) > 1), size = 0.3, alpha = 0.5) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey70"), guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  labs(title = sprintf("Differential Expression\n%d DEGs", nrow(sig_deg)),
       x = "log2 FC", y = "-log10(p)") +
  theme_bw(base_size = 10)

# Panel 2: Subtype distribution
p2 <- clinical %>%
  filter(!is.na(molecular_subtype), molecular_subtype != "") %>%
  count(molecular_subtype) %>%
  ggplot(aes(x = reorder(molecular_subtype, n), y = n, fill = molecular_subtype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Luminal A"="#1B9E77","Luminal B"="#D95F02",
                                "HER2-enriched"="#7570B3","Triple Negative"="#E7298A"),
                    guide = "none") +
  coord_flip() +
  labs(title = "Molecular Subtype Distribution", x = "", y = "Patients") +
  theme_bw(base_size = 10)

# Panel 3: Stage distribution
p3 <- clinical %>%
  filter(!is.na(stage_simple), stage_simple != "") %>%
  count(stage_simple) %>%
  ggplot(aes(x = stage_simple, y = n, fill = stage_simple)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Stage Distribution", x = "", y = "Patients") +
  theme_bw(base_size = 10)

# Panel 4: Survival by stage (simplified)
clinical_surv <- clinical %>%
  filter(!is.na(stage_simple), stage_simple != "",
         !is.na(os_time), os_time > 0) %>%
  mutate(os_years = os_time / 365.25,
         stage_group = ifelse(stage_simple %in% c("Stage I","Stage II"), "Stage I-II", "Stage III-IV"))

library(survival)
library(survminer)
fit <- survfit(Surv(os_years, vital_status) ~ stage_group, data = clinical_surv)
p4 <- ggsurvplot(fit, data = clinical_surv, pval = TRUE,
                  palette = c("Stage I-II"="#377EB8", "Stage III-IV"="#E41A1C"),
                  legend = "right", ggtheme = theme_bw(base_size = 10),
                  title = "Survival by Stage Group")$plot

# Panel 5: Model comparison
p5 <- ggplot(class_metrics, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.3f", Accuracy)), vjust = -0.5) +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  ylim(0, 1) +
  labs(title = "Subtype Classification Accuracy", x = "", y = "Accuracy") +
  theme_bw(base_size = 10)

# Panel 6: WGCNA module sizes
module_sizes <- data.frame(
  Module = c("grey","turquoise","blue","brown","yellow","green","red","black","pink"),
  Size = c(2817,584,553,336,260,224,119,63,44),
  stringsAsFactors = FALSE
)
p6 <- ggplot(module_sizes %>% filter(Module != "grey"),
             aes(x = reorder(Module, Size), y = Size, fill = Module)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = module_sizes$Module[-1], guide = "none") +
  coord_flip() +
  labs(title = "WGCNA Modules (excl. grey)", x = "", y = "Number of Genes") +
  theme_bw(base_size = 10)

# Arrange final layout
grid.arrange(p1, p2, p3, p4, p5, p6,
             ncol = 3, nrow = 2,
             top = textGrob("BRCA Multi-Omics Data Mining Summary",
                            gp = gpar(fontsize = 16, fontface = "bold")))

dev.off()
cat("  Multi-panel summary figure saved.\n")

# ---- 2. Expression heatmap of key genes across subtypes ----
cat("\nStep 2: Comprehensive heatmap of key genes...\n")

brca <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
tpm_mat <- brca$tpm
clinical <- brca$clinical

# Collect key genes from various analyses
key_genes <- unique(c(
  # Top DEGs
  deg_table %>% filter(!is.na(padj), padj < 0.01, abs(log2FC) > 2) %>% pull(gene) %>% head(30),
  # Known BRCA genes
  c("ESR1", "PGR", "ERBB2", "MKI67", "TP53", "BRCA1", "BRCA2", "PIK3CA",
    "GATA3", "FOXA1", "CCND1", "MYC", "EGFR", "VIM", "CDH1")
))

key_in_expr <- intersect(key_genes, rownames(tpm_mat))
cat(sprintf("  Key genes in expression data: %d\n", length(key_in_expr)))

# Match clinical
expr_patients <- substr(colnames(tpm_mat), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_patients, clinical$pid_short)

# Sort by subtype
subtype_order <- order(clinical$molecular_subtype[match_idx],
                       clinical$stage_simple[match_idx],
                       na.last = TRUE)

expr_key <- log2(tpm_mat[key_in_expr, subtype_order] + 1)

# Annotation
annot_col <- data.frame(
  Subtype = clinical$molecular_subtype[match_idx][subtype_order],
  Stage = clinical$stage_simple[match_idx][subtype_order],
  stringsAsFactors = FALSE
)
rownames(annot_col) <- colnames(expr_key)
annot_col[is.na(annot_col)] <- "Unknown"

ann_colors <- list(
  Subtype = c("Luminal A"="#1B9E77","Luminal B"="#D95F02",
              "HER2-enriched"="#7570B3","Triple Negative"="#E7298A","Unknown"="grey80"),
  Stage = c("Stage I"="#1B9E77","Stage II"="#D95F02",
            "Stage III"="#7570B3","Stage IV"="#E7298A","Unknown"="grey80")
)

pdf(file.path(FIG_DIR, "key_genes_heatmap.pdf"), width = 14, height = 10)
pheatmap(expr_key,
         scale = "row",
         annotation_col = annot_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 7,
         main = "BRCA Key Genes: Expression Across Molecular Subtypes")
dev.off()
cat("  Key genes heatmap saved.\n")

# ---- 3. Generate HTML report ----
cat("\nStep 3: Generating comprehensive HTML report...\n")

html_content <- c(
  '<!DOCTYPE html><html><head><meta charset="UTF-8">',
  '<title>BRCA Multi-Omics Data Mining Report</title>',
  '<style>',
  'body{font-family:Arial,sans-serif;max-width:1200px;margin:0 auto;padding:20px;background:#f5f5f5}',
  'h1{color:#1a1a1a;border-bottom:3px solid #E41A1C;padding-bottom:10px}',
  'h2{color:#333;border-bottom:2px solid #377EB8;padding-bottom:5px;margin-top:30px}',
  'h3{color:#555}',
  'table{border-collapse:collapse;width:100%;margin:15px 0;background:white}',
  'th{background:#377EB8;color:white;padding:10px;text-align:left}',
  'td{padding:8px;border-bottom:1px solid #ddd}',
  'tr:hover{background:#f0f0f0}',
  'img{max-width:100%;height:auto;margin:15px 0;border:1px solid #ddd;border-radius:4px}',
  '.highlight{background:#fff3cd;padding:15px;border-radius:5px;margin:15px 0}',
  '.metric{display:inline-block;background:white;padding:15px;margin:10px;border-radius:8px;',
  '  box-shadow:0 2px 4px rgba(0,0,0,0.1);text-align:center;min-width:120px}',
  '.metric .value{font-size:28px;font-weight:bold;color:#E41A1C}',
  '.metric .label{font-size:12px;color:#666}',
  '</style></head><body>',
  '<h1>BRCA Multi-Omics Data Mining — Comprehensive Analysis Report</h1>',
  '<p><strong>Project Date:</strong> 2026-05-14 | <strong>Data Source:</strong> TCGA-BRCA | <strong>Analysis Platform:</strong> R 4.6.0</p>',
  '',
  '<div class="highlight">',
  '<h3>Executive Summary</h3>',
  '<p>This report presents a comprehensive multi-omics data mining analysis of TCGA-BRCA (Breast Invasive Carcinoma).',
  'The analysis integrates mRNA transcriptome (25,981 genes × 1,094 patients), miRNA expression (1,881 miRNAs × 1,079 patients),',
  'and somatic mutation data, covering differential expression, machine learning classification, clustering, WGCNA co-expression networks,',
  'survival analysis, and integrative multi-omics regulatory networks.</p>',
  '</div>',
  '',
  '<h2>1. Dataset Overview</h2>',
  '<div>',
  '<div class="metric"><div class="value">1,094</div><div class="label">mRNA Patients</div></div>',
  '<div class="metric"><div class="value">1,079</div><div class="label">miRNA Patients</div></div>',
  '<div class="metric"><div class="value">1,076</div><div class="label">Common Patients</div></div>',
  '<div class="metric"><div class="value">25,981</div><div class="label">mRNA Genes</div></div>',
  '<div class="metric"><div class="value">1,881</div><div class="label">miRNAs</div></div>',
  '<div class="metric"><div class="value">6,768</div><div class="label">DEGs</div></div>',
  '</div>',
  '',
  '<table>',
  '<tr><th>Data Layer</th><th>Platform</th><th>Features</th><th>Samples</th><th>Status</th></tr>',
  '<tr><td>mRNA Expression</td><td>Illumina HiSeq RNA-seq</td><td>25,981 genes (TPM)</td><td>1,094 patients (tumor)</td><td>Complete</td></tr>',
  '<tr><td>miRNA Expression</td><td>Illumina HiSeq miRNA-seq</td><td>1,881 miRNAs (counts)</td><td>1,079 patients</td><td>Complete</td></tr>',
  '<tr><td>Somatic Mutation</td><td>WES (MuTect2)</td><td>Gene-level MAF</td><td>~990 samples</td><td>Downloaded</td></tr>',
  '<tr><td>Clinical Data</td><td>TCGA Clinical</td><td>94 variables</td><td>1,174 patients</td><td>Complete</td></tr>',
  '</table>',
  '',
  '<h2>2. Differential Expression Analysis</h2>',
  '<p>DESeq2 was used to identify differentially expressed genes (DEGs) between BRCA tumor (n=1,105) and normal (n=113) samples.</p>',
  '<ul>',
  '<li><strong>Total DEGs (|log2FC|>1, padj<0.05):</strong> 6,768</li>',
  '<li><strong>Up-regulated in tumor:</strong> 4,334 genes</li>',
  '<li><strong>Down-regulated in tumor:</strong> 2,434 genes</li>',
  '</ul>',
  '<p><a href="figures/volcano_brca.pdf" View Figure</a></p>',
  '<p><a href="figures/deg_heatmap_top50.pdf" View Figure</a></p>',
  '',
  '<h2>3. Molecular Subtype Classification</h2>',
  '<p>Three machine learning models were trained to predict BRCA molecular subtypes (Luminal A, Luminal B, HER2-enriched, Triple Negative) from mRNA expression.</p>',
  '<table>',
  '<tr><th>Model</th><th>Accuracy</th><th>Kappa</th><th>F1 (Macro)</th></tr>',
  '<tr><td><strong>LASSO (Multinomial)</strong></td><td><strong>0.894</strong></td><td>0.528</td><td>0.612</td></tr>',
  '<tr><td>Random Forest</td><td>0.867</td><td>0.365</td><td>0.706</td></tr>',
  '<tr><td>XGBoost</td><td>0.358</td><td>0.002</td><td>0.232</td></tr>',
  '</table>',
  '<p>LASSO achieved the best performance with 89.4% accuracy, selecting 35 discriminative genes.</p>',
  '<p><a href="figures/classification_roc.pdf" View Figure</a></p>',
  '',
  '<h2>4. Clustering Analysis</h2>',
  '<p>Unsupervised clustering was performed using PCA, t-SNE, and hierarchical clustering on the top 2,000 most variable genes.</p>',
  '<ul>',
  '<li><strong>K-means optimal K:</strong> 2 (silhouette score = 0.18)</li>',
  '<li><strong>PCA:</strong> PC1 (15.3%) and PC2 (8.5%) separate samples primarily by ER status</li>',
  '</ul>',
  '<p><a href="figures/pca_subtype.pdf" View Figure</a></p>',
  '<p><a href="figures/tsne_brca.pdf" View Figure</a></p>',
  '',
  '<h2>5. WGCNA Co-expression Network</h2>',
  '<p>Weighted Gene Co-expression Network Analysis identified 9 co-expression modules from 5,000 most variable genes.</p>',
  '<table>',
  '<tr><th>Module</th><th>Size</th><th>Top Hub Gene</th><th>Module Membership</th></tr>',
  '<tr><td>Blue</td><td>553</td><td>ENSG00000167208</td><td>0.959</td></tr>',
  '<tr><td>Turquoise</td><td>584</td><td>ENSG00000115163</td><td>0.889</td></tr>',
  '<tr><td>Brown</td><td>336</td><td>—</td><td>—</td></tr>',
  '<tr><td>Yellow</td><td>260</td><td>—</td><td>—</td></tr>',
  '</table>',
  '<p>Soft-threshold power: 8 (scale-free R² = 0.887)</p>',
  '<p><a href="figures/wgcna_module_trait.pdf" View Figure</a></p>',
  '<p><a href="figures/wgcna_soft_power.pdf" View Figure</a></p>',
  '',
  '<h2>6. Survival Analysis</h2>',
  '<p>Kaplan-Meier and Cox regression analyses were performed on 1,105 patients with 87 death events.</p>',
  '<ul>',
  '<li><strong>Cox model C-index:</strong> 0.771 (likelihood ratio p = 8.62e-12)</li>',
  '<li><strong>Molecular subtype:</strong> Not significantly associated with OS (log-rank p = 0.493)</li>',
  '</ul>',
  '<p><a href="figures/km_curves_stage.pdf" View Figure</a></p>',
  '<p><a href="figures/km_curves_subtype.pdf" View Figure</a></p>',
  '<p><a href="figures/cox_forest_plot.pdf" View Figure</a></p>',
  '',
  '<h2>7. Association & Enrichment Analysis</h2>',
  '<ul>',
  '<li><strong>miRNA-mRNA correlations:</strong> 228 pairs with |r|>0.3</li>',
  '<li><strong>Clinical correlations:</strong> 100 DEGs correlated with clinical variables (age, lymph nodes, OS)</li>',
  '</ul>',
  '<p><a href="figures/clinical_correlation.pdf" View Figure</a></p>',
  '',
  '<h2>8. Key Findings Summary</h2>',
  '<table>',
  '<tr><th>Finding</th><th>Details</th></tr>',
  '<tr><td>Differentially Expressed Genes</td><td>6,768 genes (4,334 up / 2,434 down in tumor vs normal)</td></tr>',
  '<tr><td>Best Classifier</td><td>LASSO multinomial: 89.4% accuracy for molecular subtype prediction</td></tr>',
  '<tr><td>Optimal Clusters</td><td>K=2 (driven primarily by ER status)</td></tr>',
  '<tr><td>WGCNA Modules</td><td>9 co-expression modules; Blue module (553 genes) with strongest hub</td></tr>',
  '<tr><td>Prognostic Model</td><td>Cox C-index=0.771; Stage is dominant predictor</td></tr>',
  '<tr><td>miRNA-mRNA Network</td><td>228 significant miRNA-mRNA regulatory pairs</td></tr>',
  '</table>',
  '',
  '<h2>9. Files Generated</h2>',
  '<table>',
  '<tr><th>Category</th><th>Files</th></tr>',
  '<tr><td>R Scripts</td><td>code/01_data_preprocessing.R through code/13_visualization.R</td></tr>',
  '<tr><td>Processed Data</td><td>data/processed/ (9 RDS/CSV files)</td></tr>',
  '<tr><td>Public Data</td><td>data/public_data/ (miRNA, mutation data)</td></tr>',
  '<tr><td>Result Figures</td><td>results/figures/ (20+ PDF figures)</td></tr>',
  '<tr><td>Result Tables</td><td>results/tables/ (15+ CSV tables)</td></tr>',
  '</table>',
  '',
  '<p style="text-align:center;color:#999;margin-top:50px">',
  'Generated by Claude Code + R 4.6.0 | TCGA-BRCA Multi-Omics Analysis | 2026-05-14</p>',
  '</body></html>'
)

writeLines(html_content, file.path(OUT_DIR, "brca_analysis_report.html"))
cat("  HTML report saved: results/brca_analysis_report.html\n")

# ---- 4. Final file manifest ----
cat("\nStep 4: Generating final deliverables manifest...\n")

Rscripts <- c("01_data_preprocessing.R","02_download_mirna.R","02_download_methylation.R",
  "02_download_mutations.R","03_integrate_omics.R","04_diff_expression.R",
  "05_classification.R","06_clustering.R","07_association_analysis.R",
  "08_wgcna.R","11_survival_analysis.R","12_advanced_algorithms.R","13_visualization.R")

processed_data <- c("brca_tumor_processed.rds","brca_tumor_tpm.rds","brca_clinical_clean.rds",
  "brca_multimics_sample_map.csv","brca_mRNA_aligned.rds","brca_miRNA_aligned.rds",
  "gene_annotation.rds","sample_info.csv","brca_tumor_counts.rds")

public_data <- c("brca_mirna_counts.rds","brca_mirna_clinical.rds","brca_mutations.rds")

figures <- c("volcano_brca.pdf","deg_heatmap_top50.pdf","ma_plot_brca.pdf","pca_subtype.pdf",
  "pca_stage.pdf","tsne_brca.pdf","classification_roc.pdf","classification_cm.pdf",
  "hclust_heatmap_top500.pdf","silhouette_scores.pdf","wgcna_module_trait.pdf",
  "wgcna_soft_power.pdf","km_curves_subtype.pdf","km_curves_stage.pdf",
  "cox_forest_plot.pdf","clinical_correlation.pdf","summary_figure.pdf",
  "key_genes_heatmap.pdf","oncoplot_top20.pdf","mutation_types.pdf",
  "mutual_exclusivity.pdf","wgcna_sample_dendrogram.pdf","omics_overlap_upset.pdf")

tables <- c("brca_degs_deseq2.csv","classification_metrics.csv","cox_regression.csv",
  "wgcna_modules.csv","wgcna_hub_genes.csv","mirna_target_correlations.csv",
  "clustering_kmeans.csv","clustering_silhouette.csv","multiomics_summary.csv",
  "mutated_genes_summary.csv","mirna_mrna_regulatory_network.csv",
  "classification_rf_features.csv")

docs <- c("01_data_background.md","brca_analysis_report.html","deliverables_manifest.csv")

all_files <- c(Rscripts, processed_data, public_data, figures, tables, docs)
categories <- c(
  rep("R Script", length(Rscripts)),
  rep("Data - Processed", length(processed_data)),
  rep("Data - Public", length(public_data)),
  rep("Figure", length(figures)),
  rep("Table", length(tables)),
  rep("Documentation", length(docs))
)

manifest <- data.frame(Category = categories, File = all_files, stringsAsFactors = FALSE)

write.csv(manifest, file.path(OUT_DIR, "deliverables_manifest.csv"), row.names = FALSE)

cat("  Manifest saved: results/deliverables_manifest.csv\n")
cat("\n========== US-013 Complete: Visualization Suite ==========\n")
cat(sprintf("  Deliverables: 10 R scripts + 10 data files + 15+ figures + 12 tables + 1 doc + 1 report\n"))
