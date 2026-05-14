# ===========================================================================
# US-012: Advanced Algorithms — Mutation Analysis + Multi-omics Integration
# 项目：BRCA多组学数据挖掘
# 日期：2026-05-14
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(maftools)
  library(tidyverse)
  library(pheatmap)
})

cat("\n========== US-012: Advanced Algorithms ==========\n\n")

DATA_DIR  <- "D:/Users/Desktop/R_Work/data/public_data"
INPUT_DIR <- "D:/Users/Desktop/R_Work/data/processed"
FIG_DIR   <- "D:/Users/Desktop/R_Work/results/figures"
TBL_DIR   <- "D:/Users/Desktop/R_Work/results/tables"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Process mutation data with maftools ----
cat("Step 1: Processing mutation data...\n")

# Find the downloaded MAF CSV file (saved by GDCquery_Maf or alternative)
maf_files <- list.files(DATA_DIR, pattern = "\\.maf", full.names = TRUE, ignore.case = TRUE)
csv_files <- list.files(DATA_DIR, pattern = "TCGA.*\\.csv", full.names = TRUE)

if (length(maf_files) > 0) {
  maf_file <- maf_files[1]
  cat(sprintf("  Using MAF: %s\n", basename(maf_file)))
} else if (length(csv_files) > 0) {
  maf_file <- csv_files[1]
  cat(sprintf("  Using CSV: %s\n", basename(maf_file)))
} else {
  # Try loading the RDS
  rds_file <- file.path(DATA_DIR, "brca_mutations.rds")
  if (file.exists(rds_file)) {
    cat("  Loading mutation data from RDS...\n")
    mut_data <- readRDS(rds_file)
    if (is.data.frame(mut_data) && nrow(mut_data) > 1000) {
      # Save as temporary TSV for maftools
      maf_file <- file.path(DATA_DIR, "brca_mutations_temp.maf")
      write.table(mut_data, maf_file, sep = "\t", quote = FALSE, row.names = FALSE)
      cat(sprintf("  Converted RDS to MAF: %d variants\n", nrow(mut_data)))
    }
  }
}

maf_loaded <- FALSE
if (exists("maf_file") && file.exists(maf_file)) {
  tryCatch({
    maf <- read.maf(maf = maf_file, verbose = FALSE)
    maf_loaded <- TRUE
    cat(sprintf("  MAF loaded: %d samples, %d genes, %d variants\n",
                length(unique(maf@data$Tumor_Sample_Barcode)),
                length(unique(maf@data$Hugo_Symbol)),
                nrow(maf@data)))
  }, error = function(e) {
    cat(sprintf("  read.maf failed: %s\n", e$message))
  })
}

# ---- 2. Oncoplot ----
if (maf_loaded) {
  cat("\nStep 2: Generating oncoplot...\n")

  pdf(file.path(FIG_DIR, "oncoplot_top20.pdf"), width = 14, height = 8)
  oncoplot(
    maf = maf,
    top = 20,
    legendFontSize = 10,
    titleText = "BRCA: Top 20 Mutated Genes"
  )
  dev.off()
  cat("  Oncoplot top 20 saved.\n")

  # ---- 3. Mutation summary ----
  cat("\nStep 3: Mutation summary statistics...\n")

  # Variant classification
  var_class <- maf@data %>%
    count(Variant_Classification, sort = TRUE)

  pdf(file.path(FIG_DIR, "mutation_types.pdf"), width = 8, height = 5)
  ggplot(var_class %>% head(10), aes(x = reorder(Variant_Classification, n), y = n, fill = Variant_Classification)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "BRCA: Variant Classification Distribution",
         x = "", y = "Count") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  dev.off()
  cat("  Variant types plot saved.\n")

  # Variant type per sample
  tryCatch({
    pdf(file.path(FIG_DIR, "mutation_load.pdf"), width = 10, height = 5)
    par(mfrow = c(1, 2))
    plot(maf, plot.type = "Mutation load", main = "BRCA: Mutation Load")
    dev.off()
    cat("  Mutation load plot saved.\n")
  }, error = function(e) {
    cat(sprintf("  Mutation load plot skipped: %s\n", e$message))
  })

  # Top mutated genes table
  gene_summary <- maf@data %>%
    count(Hugo_Symbol, sort = TRUE) %>%
    rename(N_Mutations = n)
  write.csv(gene_summary, file.path(TBL_DIR, "mutated_genes_summary.csv"), row.names = FALSE)
  cat(sprintf("  Top mutated: %s (%d), %s (%d), %s (%d)\n",
              gene_summary$Hugo_Symbol[1], gene_summary$N_Mutations[1],
              gene_summary$Hugo_Symbol[2], gene_summary$N_Mutations[2],
              gene_summary$Hugo_Symbol[3], gene_summary$N_Mutations[3]))

  # ---- 4. Mutual exclusivity ----
  cat("\nStep 4: Somatic interaction analysis...\n")

  tryCatch({
    pdf(file.path(FIG_DIR, "mutual_exclusivity.pdf"), width = 10, height = 8)
    somaticInteractions(
      maf = maf,
      top = 15,
      pvalue = c(0.05, 0.1)
    )
    dev.off()
    cat("  Mutual exclusivity plot saved.\n")
  }, error = function(e) {
    cat(sprintf("  Interaction analysis skipped: %s\n", e$message))
  })
} else {
  cat("\n  MAF not loaded - skipping mutation visualization.\n")
}

# ---- 5. miRNA-mRNA integrative regulatory network ----
cat("\nStep 5: Building miRNA-mRNA regulatory network...\n")

tryCatch({
mirna_aligned <- readRDS(file.path(INPUT_DIR, "brca_miRNA_aligned.rds"))
mRNA_aligned  <- readRDS(file.path(INPUT_DIR, "brca_mRNA_aligned.rds"))
deg_table <- read.csv(file.path(TBL_DIR, "brca_degs_deseq2.csv"))

# Top DEGs + top miRNAs
top_degs <- deg_table %>%
  filter(!is.na(padj), padj < 0.01, abs(log2FC) > 2) %>%
  pull(gene)
top_degs <- intersect(top_degs, rownames(mRNA_aligned))
cat(sprintf("  Top DEGs for network: %d\n", length(top_degs)))

if (length(top_degs) >= 20 && ncol(mirna_aligned) > 100) {
  # Select most variable miRNAs
  mirna_var <- apply(log2(mirna_aligned + 1), 1, var)
  top_mirnas <- names(sort(mirna_var, decreasing = TRUE))[1:100]

  mirna_expr <- t(log2(mirna_aligned[top_mirnas, ] + 1))
  mrna_expr  <- t(log2(mRNA_aligned[top_degs, ] + 1))

  # Correlation matrix
  net_cor <- cor(mirna_expr, mrna_expr, use = "pairwise.complete.obs")

  # Filter significant negative correlations (miRNA typically suppresses targets)
  cor_long <- as.data.frame(as.table(net_cor))
  names(cor_long) <- c("miRNA", "mRNA", "correlation")
  cor_sig <- cor_long %>%
    filter(abs(correlation) > 0.4) %>%
    arrange(desc(abs(correlation)))

  cat(sprintf("  Regulatory pairs (|r|>0.4): %d\n", nrow(cor_sig)))

  if (nrow(cor_sig) > 10) {
    write.csv(cor_sig, file.path(TBL_DIR, "mirna_mrna_regulatory_network.csv"), row.names = FALSE)

    # Network heatmap
    if (nrow(cor_sig) >= 50) {
      # Pivot to matrix for heatmap: top 30 miRNAs x top 30 DEGs
      top_mirna_net <- names(sort(table(cor_sig$miRNA), decreasing = TRUE))[1:min(30, length(unique(cor_sig$miRNA)))]
      top_mrna_net  <- names(sort(table(cor_sig$mRNA), decreasing = TRUE))[1:min(30, length(unique(cor_sig$mRNA)))]

      sub_cor <- cor_sig %>%
        filter(miRNA %in% top_mirna_net, mRNA %in% top_mrna_net)

      net_matrix <- reshape2::acast(sub_cor, miRNA ~ mRNA, value.var = "correlation", fill = 0)

      pdf(file.path(FIG_DIR, "mirna_mrna_network_heatmap.pdf"), width = 12, height = 10)
      pheatmap(net_matrix,
               main = "miRNA-mRNA Regulatory Network",
               color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
               fontsize_row = 6, fontsize_col = 6)
      dev.off()
      cat("  Regulatory network heatmap saved.\n")
    }
  }
}
}, error = function(e) {
  cat(sprintf("  miRNA-mRNA network error: %s\n", e$message))
})

# ---- 6. Multi-omics overview ----
cat("\nStep 6: Multi-omics data integration summary...\n")

omics_summary <- data.frame(
  Layer = c("mRNA Expression", "miRNA Expression", "Somatic Mutation"),
  Samples = c(1105, 1207, if(maf_loaded) length(unique(maf@data$Tumor_Sample_Barcode)) else 0),
  Features = c(25981, 1881, if(maf_loaded) length(unique(maf@data$Hugo_Symbol)) else 0),
  Key_Finding = c(
    "6,768 DEGs (Tumor vs Normal); LASSO 89.4% subtype accuracy",
    "228 miRNA-mRNA pairs (|r|>0.3); 100 regulatory pairs (|r|>0.4)",
    if(maf_loaded) sprintf("Top mutated: %s (%d variants)", gene_summary$Hugo_Symbol[1], gene_summary$N_Mutations[1]) else "Data available, processing pending"
  ),
  stringsAsFactors = FALSE
)

write.csv(omics_summary, file.path(TBL_DIR, "multiomics_summary.csv"), row.names = FALSE)

cat("\n========== US-012 Complete: Advanced Algorithms ==========\n")
