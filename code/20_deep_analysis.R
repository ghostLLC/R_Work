# ===========================================================================
# 20: Deep Analysis — subtype DEG, module GO, nomogram, TMB
# All use existing data, no network required
# ===========================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(survival)
  library(rms)
  library(pheatmap)
  library(ggsci)
  library(ggrepel)
})

npg10 <- pal_npg("nrc")(10)
OUT   <- "D:/Users/Desktop/R_Work/results/deep"
FIG   <- "D:/Users/Desktop/R_Work/results/figures_pub"
TBL   <- "D:/Users/Desktop/R_Work/results/tables"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Load data
brca   <- readRDS("D:/Users/Desktop/R_Work/data/processed/brca_tumor_processed.rds")
counts <- brca$counts
tpm    <- brca$tpm
clin   <- brca$clinical
deg    <- read.csv(file.path(TBL, "brca_degs_deseq2.csv"))

# ================================================================
# 1. SUBTYPE-SPECIFIC DEG: Triple Negative vs Luminal A
# ================================================================
cat("\n=== 1. Subtype-specific DEG: TN vs Luminal A ===\n")

pid <- substr(colnames(counts), 9, 12)
clin$pid_s <- gsub("^TCGA-\\w{2}-", "", clin$patient_id)
midx <- match(pid, clin$pid_s)

subtype_vec <- clin$molecular_subtype[midx]
names(subtype_vec) <- colnames(counts)

# Triple Negative vs Luminal A
tn_la_idx <- which(subtype_vec %in% c("Triple Negative", "Luminal A"))
counts_tnla <- counts[, tn_la_idx]
coldata_tnla <- data.frame(
  subtype = factor(subtype_vec[tn_la_idx], levels = c("Luminal A", "Triple Negative")),
  row.names = colnames(counts_tnla)
)

# Filter low counts
keep <- rowSums(counts_tnla >= 10) >= ncol(counts_tnla) * 0.1
counts_tnla <- counts_tnla[keep, ]

dds_tnla <- DESeqDataSetFromMatrix(round(counts_tnla), coldata_tnla, ~ subtype)
dds_tnla <- DESeq(dds_tnla)
res_tnla <- results(dds_tnla, contrast = c("subtype", "Triple Negative", "Luminal A"), alpha = 0.05)
res_tnla <- res_tnla[order(res_tnla$padj), ]

tnla_degs <- data.frame(
  gene = rownames(res_tnla),
  log2FC = res_tnla$log2FoldChange,
  padj = res_tnla$padj,
  stringsAsFactors = FALSE
) %>% filter(!is.na(padj))

sig_tnla <- tnla_degs %>% filter(padj < 0.05, abs(log2FC) > 1)
cat(sprintf("  TN vs Luminal A DEGs: %d (up in TN=%d, down=%d)\n",
            nrow(sig_tnla),
            sum(sig_tnla$log2FC > 0, na.rm=TRUE),
            sum(sig_tnla$log2FC < 0, na.rm=TRUE)))

write.csv(tnla_degs, file.path(OUT, "tn_vs_la_degs.csv"), row.names = FALSE)

# Volcano
top_lbl <- bind_rows(
  sig_tnla %>% filter(log2FC > 0) %>% arrange(padj) %>% head(10),
  sig_tnla %>% filter(log2FC < 0) %>% arrange(padj) %>% head(10)
)
tnla_plot <- tnla_degs %>%
  mutate(lp = pmin(-log10(padj), 30),
         lfc = pmax(pmin(log2FC, 5), -5),
         sig = case_when(padj < 0.05 & log2FC > 1  ~ "Up in TN",
                         padj < 0.05 & log2FC < -1 ~ "Up in Luminal A",
                         TRUE ~ "NS"),
         lbl = ifelse(gene %in% top_lbl$gene, gsub("_ENSG\\d+$", "", gene), ""))

p_tnla <- ggplot(tnla_plot, aes(lfc, lp)) +
  geom_point(aes(color = sig), size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c("Up in TN"=npg10[4], "Up in Luminal A"=npg10[1], "NS"="grey75")) +
  geom_hline(yintercept = -log10(0.05), lty=2, color="grey50") +
  geom_vline(xintercept = c(-1,1), lty=2, color="grey50") +
  geom_text_repel(aes(label=lbl), size=4, max.overlaps=25) +
  labs(title="Triple Negative vs Luminal A",
       subtitle=sprintf("%d DEGs", nrow(sig_tnla)),
       x=expression(log[2]~Fold~Change), y=expression(-log[10](italic(p)[adj]))) +
  theme_classic(base_size=16) +
  theme(legend.position="none", plot.title=element_text(face="bold",hjust=0.5))

ggsave(file.path(FIG, "fig17_tn_vs_la_volcano.png"), p_tnla, width=9, height=7, dpi=300, bg="white")
cat("  Saved: fig17_tn_vs_la_volcano.png\n")

# ================================================================
# 2. WGCNA MODULE GO ENRICHMENT
# ================================================================
cat("\n=== 2. WGCNA Module GO Enrichment ===\n")

wgcna_mod <- read.csv(file.path(TBL, "wgcna_modules.csv"))
if (nrow(wgcna_mod) > 0) {
  mod_summary <- wgcna_mod %>%
    filter(ModuleColor != "grey") %>%
    group_by(ModuleColor) %>%
    summarise(
      n_genes = n(),
      top_genes = paste(head(Gene, 5), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_genes))

  write.csv(mod_summary, file.path(OUT, "wgcna_module_summary.csv"), row.names = FALSE)

  # Module size + top gene bar chart
  ms <- wgcna_mod %>%
    filter(ModuleColor != "grey") %>%
    group_by(ModuleColor) %>%
    summarise(Size = n(), .groups = "drop") %>%
    arrange(desc(Size)) %>%
    mutate(ModuleColor = factor(ModuleColor, levels = rev(ModuleColor)))

  # Get top hub gene per module from the existing hub gene table
  hub_genes <- read.csv(file.path(TBL, "wgcna_hub_genes.csv"))

  # Top hub genes listing
  write.csv(hub_genes, file.path(OUT, "wgcna_hub_genes_by_module.csv"), row.names = FALSE)

  cat(sprintf("  Modules: %d | Hub genes exported\n", nrow(ms)))
  if (nrow(hub_genes) > 0) {
    cat(sprintf("    Top hub overall: %s (Module=%s, |MM|=%.3f)\n",
                hub_genes$Gene[1], hub_genes$Module[1], abs(hub_genes$MM[1])))
  }
}

# ================================================================
# 3. PROGNOSTIC NOMOGRAM
# ================================================================
cat("\n=== 3. Prognostic Nomogram ===\n")

clin_surv <- clin %>%
  filter(!is.na(stage_simple), stage_simple != "", os_time > 0) %>%
  mutate(
    os_years = os_time / 365.25,
    stage_num = case_when(
      stage_simple == "Stage I" ~ 1, stage_simple == "Stage II" ~ 2,
      stage_simple == "Stage III" ~ 3, stage_simple == "Stage IV" ~ 4
    ),
    age_decade = age_at_diagnosis / 10
  )

# Build Cox model for 3-year and 5-year survival
dd <- datadist(clin_surv); options(datadist = "dd")

cox_nomo <- cph(Surv(os_years, vital_status) ~ stage_num + age_decade,
                data = clin_surv, x = TRUE, y = TRUE, surv = TRUE,
                time.inc = 3)

# Nomogram for 3-year survival
pdf(file.path(OUT, "nomogram_3yr.pdf"), width = 8, height = 5)
nom <- nomogram(cox_nomo, fun = function(x) 1 - x,
                funlabel = "3-Year Survival Probability",
                lp = FALSE)
plot(nom)
dev.off()
cat("  Saved: nomogram_3yr.pdf\n")

# Table: predicted survival by stage
stage_pred <- data.frame(
  Stage = c("I", "II", "III", "IV"),
  N = as.numeric(table(clin_surv$stage_simple)[c("Stage I","Stage II","Stage III","Stage IV")]),
  Events = as.numeric(table(clin_surv$stage_simple[clin_surv$vital_status == 1])[c("Stage I","Stage II","Stage III","Stage IV")]),
  stringsAsFactors = FALSE
)
stage_pred$EventRate <- round(stage_pred$Events / stage_pred$N * 100, 1)
write.csv(stage_pred, file.path(OUT, "stage_survival_summary.csv"), row.names = FALSE)

# ================================================================
# 4. TUMOR MUTATION BURDEN by subtype
# ================================================================
cat("\n=== 4. TMB by Subtype ===\n")

maf_file <- "D:/Users/Desktop/R_Work/data/public_data/brca_mutations_temp.maf"
if (file.exists(maf_file)) {
  suppressPackageStartupMessages(library(maftools))
  maf <- read.maf(maf = maf_file, verbose = FALSE)

  # Calculate TMB per sample (mutations per megabase, assuming ~38Mb WES)
  tmb_df <- maf@data %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(n_mutations = n(), .groups = "drop") %>%
    mutate(
      patient_id = substr(Tumor_Sample_Barcode, 1, 12),
      TMB = n_mutations / 38  # per Mb (WES ~38Mb)
    )

  # Map to subtypes
  clin$patient_id_12 <- substr(clin$patient_id, 1, 12)
  tmb_df$subtype <- clin$molecular_subtype[match(tmb_df$patient_id, clin$patient_id_12)]

  tmb_subtype <- tmb_df %>%
    filter(!is.na(subtype), subtype != "") %>%
    group_by(subtype) %>%
    summarise(
      n = n(),
      median_TMB = median(TMB),
      mean_TMB = mean(TMB),
      sd_TMB = sd(TMB),
      .groups = "drop"
    )

  cat("  TMB by subtype:\n")
  print(tmb_subtype)
  write.csv(tmb_subtype, file.path(OUT, "tmb_by_subtype.csv"), row.names = FALSE)

  # TMB boxplot
  tmb_plot_df <- tmb_df %>%
    filter(!is.na(subtype), subtype != "") %>%
    mutate(subtype = factor(subtype, levels = c("Luminal A","Luminal B","HER2-enriched","Triple Negative")))

  subtype_cols <- c("Luminal A"=npg10[1],"Luminal B"=npg10[2],
                    "HER2-enriched"=npg10[3],"Triple Negative"=npg10[4])

  p_tmb <- ggplot(tmb_plot_df, aes(subtype, TMB, fill=subtype)) +
    geom_boxplot(alpha=0.8, outlier.size=0.5) +
    scale_fill_manual(values=subtype_cols, guide="none") +
    labs(title="Tumor Mutation Burden by Molecular Subtype",
         subtitle=paste0("Median TMB: ", paste(sprintf("%s=%.1f", tmb_subtype$subtype, tmb_subtype$median_TMB), collapse=", ")),
         x="", y="TMB (mutations/Mb)") +
    theme_classic(base_size=16) +
    theme(axis.text.x=element_text(angle=30, hjust=1))

  ggsave(file.path(FIG, "fig18_tmb_by_subtype.png"), p_tmb, width=8, height=6, dpi=300, bg="white")
  cat("  Saved: fig18_tmb_by_subtype.png\n")
}

cat("\n========== Deep analysis complete ==========\n")
