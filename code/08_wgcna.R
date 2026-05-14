# ===========================================================================
# US-010: WGCNA Co-expression Network Analysis
# 项目：BRCA多组学数据挖掘
# 日期：2026-05-14
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(pheatmap)
})

# WGCNA settings
enableWGCNAThreads(4)
options(stringsAsFactors = FALSE)

cat("\n========== US-010: WGCNA Co-expression Network ==========\n\n")

INPUT_DIR  <- "D:/Users/Desktop/R_Work/data/processed"
FIG_DIR    <- "D:/Users/Desktop/R_Work/results/figures"
TBL_DIR    <- "D:/Users/Desktop/R_Work/results/tables"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load and prepare data ----
cat("Step 1: Loading and preparing expression data...\n")

brca <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
tpm_mat <- brca$tpm
clinical <- brca$clinical

cat(sprintf("  TPM: %d genes x %d samples\n", nrow(tpm_mat), ncol(tpm_mat)))

# Select top 5000 most variable genes for WGCNA
log_tpm <- log2(tpm_mat + 1)
gene_var <- apply(log_tpm, 1, var)
top5000 <- names(sort(gene_var, decreasing = TRUE))[1:min(5000, nrow(log_tpm))]

datExpr <- t(log_tpm[top5000, ])
cat(sprintf("  WGCNA input: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))

# ---- 2. Sample clustering to detect outliers ----
cat("\nStep 2: Sample clustering for outlier detection...\n")

gsg <- goodSamplesGenes(datExpr, verbose = 3)
cat(sprintf("  All genes OK: %s | All samples OK: %s\n", gsg$allOK, gsg$allOK))

if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat(sprintf("  After filtering: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))
}

# Sample dendrogram
sampleTree <- hclust(dist(datExpr), method = "average")

pdf(file.path(FIG_DIR, "wgcna_sample_dendrogram.pdf"), width = 12, height = 6)
par(cex = 0.5)
plot(sampleTree, main = "BRCA Sample Clustering", sub = "", xlab = "")
dev.off()

# ---- 3. Soft-threshold power selection ----
cat("\nStep 3: Selecting soft-threshold power...\n")

powers <- c(1:20)
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 5,
  networkType = "signed"
)

pdf(file.path(FIG_DIR, "wgcna_soft_power.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.8, col = "blue", lty = 2)

# Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

par(mfrow = c(1, 1))
dev.off()

# Select power
soft_power <- sft$powerEstimate
if (is.na(soft_power)) {
  soft_power <- 6
  cat(sprintf("  No optimal power found, using default: %d\n", soft_power))
} else {
  cat(sprintf("  Optimal soft power: %d (R^2=%.3f)\n",
              soft_power, sft$fitIndices[soft_power, 2]))
}

# ---- 4. Blockwise module construction ----
cat("\nStep 4: Blockwise module construction (this may take several minutes)...\n")

net <- blockwiseModules(
  datExpr,
  power = soft_power,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 2,
  maxBlockSize = 5000,
  nThreads = 4
)

module_colors <- labels2colors(net$colors)
cat(sprintf("  Modules identified: %d\n", length(unique(module_colors))))

# Save module assignments
module_df <- data.frame(
  Gene = colnames(datExpr),
  Module = net$colors,
  ModuleColor = module_colors,
  stringsAsFactors = FALSE
)
write.csv(module_df, file.path(TBL_DIR, "wgcna_modules.csv"), row.names = FALSE)

# Module size distribution
mod_sizes <- table(module_colors)
cat(sprintf("  Module sizes: %s\n",
            paste(names(sort(mod_sizes, decreasing = TRUE))[1:min(10, length(mod_sizes))],
                  sort(mod_sizes, decreasing = TRUE)[1:min(10, length(mod_sizes))],
                  sep = "=", collapse = ", ")))

# ---- 5. Module-trait correlation ----
cat("\nStep 5: Module-trait correlation...\n")

# Prepare clinical traits
expr_patients <- substr(rownames(datExpr), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_patients, clinical$pid_short)

traits <- data.frame(
  Age         = as.numeric(clinical$age_at_diagnosis[match_idx]),
  LymphNodes  = as.numeric(clinical$positive_lymph_nodes[match_idx]),
  OS_Time     = as.numeric(clinical$os_time[match_idx]),
  OS_Status   = as.numeric(clinical$vital_status[match_idx]),
  stringsAsFactors = FALSE
)

# Encode subtypes as binary
for (st in unique(clinical$molecular_subtype[match_idx])) {
  if (!is.na(st) && st != "") {
    traits[[paste0("Subtype_", gsub("[ -]", "_", st))]] <-
      as.numeric(clinical$molecular_subtype[match_idx] == st)
  }
}

# Encode stages
for (stg in unique(clinical$stage_simple[match_idx])) {
  if (!is.na(stg) && stg != "") {
    traits[[paste0("Stage_", stg)]] <-
      as.numeric(clinical$stage_simple[match_idx] == stg)
  }
}

# Replace NAs with median
for (i in seq_along(traits)) {
  nas <- is.na(traits[[i]])
  if (any(nas)) {
    traits[[i]][nas] <- median(traits[[i]], na.rm = TRUE)
  }
}

# Module eigengenes
MEs <- net$MEs
colnames(MEs) <- gsub("^ME", "M", colnames(MEs))

# Correlate modules with traits
module_trait_cor <- cor(MEs, traits, use = "p")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(datExpr))

# Heatmap
pdf(file.path(FIG_DIR, "wgcna_module_trait.pdf"), width = 14, height = 10)

text_matrix <- matrix(sprintf("%.3f", module_trait_cor),
                      nrow = nrow(module_trait_cor),
                      dimnames = dimnames(module_trait_cor))

pheatmap(module_trait_cor,
         main = "WGCNA Module-Trait Correlations",
         display_numbers = text_matrix, fontsize_number = 6,
         fontsize = 8,
         color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = TRUE,
         cluster_cols = TRUE)

dev.off()
cat("  Module-trait heatmap saved.\n")

# ---- 6. Hub gene identification ----
cat("\nStep 6: Identifying hub genes in key modules...\n")

# Find modules most correlated with survival
surv_cor <- abs(module_trait_cor[, "OS_Status"])
sig_modules <- names(sort(surv_cor[surv_cor > 0.1], decreasing = TRUE))

if (length(sig_modules) == 0) {
  # Fallback: take top 3 modules by size
  sig_modules <- names(sort(mod_sizes, decreasing = TRUE))[1:3]
}

cat(sprintf("  Key modules: %s\n", paste(sig_modules, collapse = ", ")))

hub_genes_list <- list()

# Build color-to-number mapping
color_to_num <- setNames(unique(net$colors), labels2colors(unique(net$colors)))
color_to_num <- color_to_num[!duplicated(names(color_to_num))]

# Calculate gene-module membership (kME) for key modules
for (mod_color in sig_modules) {
  mod_num <- color_to_num[mod_color]
  if (is.na(mod_num)) {
    # Try extracting number from M-prefixed name
    mod_num <- as.numeric(gsub("M", "", mod_color))
  }
  mod_genes <- colnames(datExpr)[net$colors == mod_num]

  # Build ME column name (MEs columns use color names)
  me_col <- paste0("M", mod_num)

  if (length(mod_genes) >= 20) {
    mod_expr <- datExpr[, mod_genes, drop = FALSE]

    # Module membership
    kME <- cor(mod_expr, MEs[, me_col, drop = FALSE], use = "p")
    colnames(kME) <- "MM"

    hub_df <- data.frame(
      Gene = rownames(kME),
      MM = kME[, 1],
      stringsAsFactors = FALSE
    ) %>% arrange(desc(abs(MM))) %>% head(20)

    hub_genes_list[[mod_color]] <- hub_df
    cat(sprintf("    %s: %d genes, top hub=%s (MM=%.3f)\n",
                mod_color, length(mod_genes), hub_df$Gene[1], abs(hub_df$MM[1])))
  }
}

# Combine hub genes
hub_all <- bind_rows(hub_genes_list, .id = "Module")
write.csv(hub_all, file.path(TBL_DIR, "wgcna_hub_genes.csv"), row.names = FALSE)

# ---- 7. Network visualization of top module ----
cat("\nStep 7: Exporting network of key module for visualization...\n")

if (length(sig_modules) > 0) {
  top_mod_color <- sig_modules[1]
  top_mod_num <- color_to_num[top_mod_color]
  if (is.na(top_mod_num)) top_mod_num <- as.numeric(gsub("M", "", top_mod_color))
  top_genes <- colnames(datExpr)[net$colors == top_mod_num]

  me_idx <- which(module_colors == top_mod_color)[1]
  if (is.na(me_idx)) me_idx <- top_mod_num + 1
  me_col <- colnames(MEs)[me_idx]

  if (length(top_genes) >= 30) {
    # For large modules, take top 50 hub genes
    kME_all <- cor(datExpr[, top_genes], MEs[, me_col, drop = FALSE], use = "p")
    top_hub_idx <- order(abs(kME_all[, 1]), decreasing = TRUE)[1:min(50, length(top_genes))]
    top_hub_genes <- top_genes[top_hub_idx]

    # Export edges for Cytoscape
    top_expr <- datExpr[, top_hub_genes]
    adj <- abs(cor(top_expr, use = "p"))^soft_power
    diag(adj) <- 0

    # Threshold edges
    adj[adj < 0.3] <- 0

    # Save as edge list
    edges <- which(adj > 0, arr.ind = TRUE)
    if (nrow(edges) > 0) {
      edge_df <- data.frame(
        Source = rownames(adj)[edges[, 1]],
        Target = colnames(adj)[edges[, 2]],
        Weight = adj[edges],
        stringsAsFactors = FALSE
      )
      edge_df <- edge_df %>% filter(Source < Target) %>% arrange(desc(Weight))
      write.csv(edge_df, file.path(TBL_DIR, "wgcna_network_edges.csv"), row.names = FALSE)
      cat(sprintf("  Network edges exported: %d edges for %d genes\n",
                  nrow(edge_df), length(top_hub_genes)))
    }

    # Simple network heatmap
    if (sum(adj > 0) > 10) {
      pdf(file.path(FIG_DIR, "wgcna_network.pdf"), width = 10, height = 10)
      tryCatch({
        pheatmap(adj,
                 main = sprintf("WGCNA %s Module: Gene Co-expression Network", top_mod_color),
                 show_rownames = (ncol(adj) <= 30),
                 show_colnames = (ncol(adj) <= 30),
                 color = colorRampPalette(c("white", "#FFA500", "#E41A1C"))(50))
      }, error = function(e) {
        cat(sprintf("  Network heatmap failed: %s\n", e$message))
      })
      dev.off()
      cat("  Network heatmap saved.\n")
    }
  }
}

cat("\n========== US-010 Complete: WGCNA Analysis ==========\n")
cat(sprintf("  Total modules: %d | Soft power: %d\n", length(unique(module_colors)), soft_power))
