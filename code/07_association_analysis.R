# ===========================================================================
# US-009: Association Analysis — GSEA, Pathway Enrichment, Clinical Correlations
# 项目：BRCA多组学数据挖掘
# 日期：2026-05-14
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(pheatmap)
})

cat("\n========== US-009: Association & Enrichment Analysis ==========\n\n")

INPUT_DIR  <- "D:/Users/Desktop/R_Work/data/processed"
FIG_DIR    <- "D:/Users/Desktop/R_Work/results/figures"
TBL_DIR    <- "D:/Users/Desktop/R_Work/results/tables"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load DEG results ----
cat("Step 1: Loading differential expression results...\n")

deg_table <- read.csv(file.path(TBL_DIR, "brca_degs_deseq2.csv"))
cat(sprintf("  DEG table: %d genes\n", nrow(deg_table)))

sig_degs <- deg_table %>% filter(!is.na(padj), padj < 0.05, abs(log2FC) > 1)
cat(sprintf("  Significant DEGs (|log2FC|>1, padj<0.05): %d (%d up, %d down)\n",
            nrow(sig_degs), sum(sig_degs$log2FC > 0), sum(sig_degs$log2FC < 0)))

# ---- 2. GO enrichment for up-regulated genes ----
cat("\nStep 2: GO enrichment analysis...\n")

# Convert gene symbols to ENTREZ IDs
up_genes <- sig_degs$gene[sig_degs$log2FC > 0]
down_genes <- sig_degs$gene[sig_degs$log2FC < 0]

# Strip version suffixes for mapping
map_genes <- function(genes) {
  clean <- gsub("_ENSG\\d+$", "", genes)
  bitr(clean, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db, drop = TRUE)
}

up_entrez <- tryCatch(map_genes(up_genes), error = function(e) NULL)
down_entrez <- tryCatch(map_genes(down_genes), error = function(e) NULL)

if (!is.null(up_entrez) && nrow(up_entrez) > 10) {
  go_bp_up <- enrichGO(
    gene          = up_entrez$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.1
  )

  if (!is.null(go_bp_up) && nrow(go_bp_up) > 0) {
    # Dot plot
    pdf(file.path(FIG_DIR, "go_bp_dotplot.pdf"), width = 10, height = 8)
    print(dotplot(go_bp_up, showCategory = 20, title = "GO Biological Process: Up-regulated DEGs"))
    dev.off()
    cat(sprintf("  GO-BP up: %d enriched terms\n", nrow(go_bp_up)))

    write.csv(go_bp_up@result, file.path(TBL_DIR, "go_bp_up.csv"), row.names = FALSE)
  }
}

# ---- 3. KEGG pathway enrichment ----
cat("\nStep 3: KEGG pathway enrichment...\n")

all_entrez <- tryCatch(map_genes(deg_table$gene), error = function(e) NULL)

if (!is.null(all_entrez) && nrow(all_entrez) > 10) {
  # Rank genes by log2FC for GSEA
  gene_list <- sig_degs$log2FC
  names(gene_list) <- gsub("_ENSG\\d+$", "", sig_degs$gene)

  # GSEA with hallmark gene sets (MSigDB via clusterProfiler)
  tryCatch({
    # Use KEGG for standard enrichment
    kegg_up <- enrichKEGG(
      gene         = up_entrez$ENTREZID,
      organism     = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )

    if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
      pdf(file.path(FIG_DIR, "kegg_enrichment.pdf"), width = 10, height = 8)
      print(dotplot(kegg_up, showCategory = 20, title = "KEGG Pathways: Up-regulated DEGs"))
      dev.off()
      write.csv(kegg_up@result, file.path(TBL_DIR, "kegg_up.csv"), row.names = FALSE)
      cat(sprintf("  KEGG: %d enriched pathways\n", nrow(kegg_up)))
    }
  }, error = function(e) {
    cat(sprintf("  KEGG enrichment failed: %s\n", e$message))
  })
}

# ---- 4. GSEA with ranked gene list ----
cat("\nStep 4: GSEA analysis with ranked gene list...\n")

# Prepare ranked gene list
gsea_genes <- sig_degs$log2FC
names(gsea_genes) <- gsub("_ENSG\\d+$", "", sig_degs$gene)
gsea_genes <- sort(gsea_genes, decreasing = TRUE)

# Remove duplicates
gsea_genes <- gsea_genes[!duplicated(names(gsea_genes))]

tryCatch({
  # Hallmark gene sets from MSigDB
  hallmark_gmt <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt"

  # Use local MSigDB if available, otherwise skip
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    library(msigdbr)
    h_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
    h_list <- split(h_gene_sets$gene_symbol, h_gene_sets$gs_name)

    gsea_res <- GSEA(
      geneList     = gsea_genes,
      TERM2GENE    = h_gene_sets[, c("gs_name", "gene_symbol")],
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )

    if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
      write.csv(gsea_res@result, file.path(TBL_DIR, "gsea_hallmark.csv"), row.names = FALSE)

      pdf(file.path(FIG_DIR, "gsea_hallmark.pdf"), width = 12, height = 10)
      print(dotplot(gsea_res, showCategory = 20, title = "GSEA Hallmark: BRCA Tumor vs Normal"))
      dev.off()
      cat(sprintf("  GSEA Hallmark: %d enriched gene sets\n", nrow(gsea_res)))

      # Ridge plot
      if (nrow(gsea_res) >= 5) {
        pdf(file.path(FIG_DIR, "gsea_ridgeplot.pdf"), width = 12, height = 8)
        print(ridgeplot(gsea_res, showCategory = 15))
        dev.off()
      }
    } else {
      cat("  GSEA: No significant results\n")
    }
  } else {
    cat("  msigdbr not installed, skipping GSEA. Will use GO/KEGG enrichment only.\n")

    # Write placeholder
    cat("GSEA skipped: msigdbr not available\n",
        file = file.path(TBL_DIR, "gsea_hallmark.csv"))
  }
}, error = function(e) {
  cat(sprintf("  GSEA failed: %s\n", e$message))
})

# ---- 5. Clinical-molecular correlation ----
cat("\nStep 5: Clinical-molecular correlation analysis...\n")

brca <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
clinical <- brca$clinical
tpm_mat <- brca$tpm

# Get top DEGs for correlation
top_degs <- sig_degs %>% arrange(padj) %>% head(100)
top_degs_clean <- gsub("_ENSG\\d+$", "", top_degs$gene)

# Match expression samples to clinical
expr_patients <- substr(colnames(tpm_mat), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_patients, clinical$pid_short)

# Build correlation matrix: genes x clinical variables
clinical_num <- clinical %>%
  dplyr::select(age_at_diagnosis, positive_lymph_nodes, os_time, vital_status) %>%
  mutate_all(as.numeric)

# Select numeric clinical variables that are complete
clin_vars <- clinical_num[match_idx, ]
rownames(clin_vars) <- colnames(tpm_mat)

# Expression of top DEGs
expr_top <- log2(t(tpm_mat[rownames(tpm_mat) %in% top_degs$gene, ]) + 1)

# Correlation between top DEGs and clinical variables
if (ncol(expr_top) > 0 && nrow(clin_vars) > 0) {
  common <- intersect(rownames(expr_top), rownames(clin_vars))
  cor_matrix <- cor(expr_top[common, ], clin_vars[common, ],
                    use = "pairwise.complete.obs")

  pdf(file.path(FIG_DIR, "clinical_correlation.pdf"), width = 10, height = 12)
  pheatmap(cor_matrix,
           main = "Top 100 DEGs vs Clinical Variables",
           show_rownames = FALSE,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
           breaks = seq(-1, 1, length.out = 101))
  dev.off()
  cat(sprintf("  Clinical correlation heatmap saved (%d genes x %d vars)\n",
              nrow(cor_matrix), ncol(cor_matrix)))
}

# ---- 6. miRNA-mRNA correlation ----
cat("\nStep 6: miRNA-mRNA target correlation...\n")

mirna_aligned <- readRDS(file.path(INPUT_DIR, "brca_miRNA_aligned.rds"))
mRNA_aligned  <- readRDS(file.path(INPUT_DIR, "brca_mRNA_aligned.rds"))

if (ncol(mirna_aligned) == ncol(mRNA_aligned) && ncol(mirna_aligned) > 100) {
  # Correlate top 50 miRNAs with top 50 DEGs
  mirna_var <- apply(log2(mirna_aligned + 1), 1, var)
  mrna_var <- apply(log2(mRNA_aligned + 1), 1, var)

  top_mirna <- names(sort(mirna_var, decreasing = TRUE))[1:50]
  top_mrna  <- names(sort(mrna_var, decreasing = TRUE))[1:50]

  mirna_sub <- t(log2(mirna_aligned[top_mirna, ] + 1))
  mrna_sub  <- t(log2(mRNA_aligned[top_mrna, ] + 1))

  mirna_mrna_cor <- cor(mirna_sub, mrna_sub, use = "pairwise.complete.obs")

  # Save significant correlations
  cor_long <- as.data.frame(as.table(mirna_mrna_cor))
  names(cor_long) <- c("miRNA", "mRNA", "correlation")
  cor_long <- cor_long %>% filter(abs(correlation) > 0.3) %>% arrange(desc(abs(correlation)))

  write.csv(cor_long, file.path(TBL_DIR, "mirna_target_correlations.csv"), row.names = FALSE)
  cat(sprintf("  miRNA-mRNA correlations saved: %d pairs (|r|>0.3)\n", nrow(cor_long)))
}

cat("\n========== US-009 Complete: Association Analysis ==========\n")
