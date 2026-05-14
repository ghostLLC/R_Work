# ===========================================================================
# 模块1：BRCA数据预处理
# 项目：BRCA深度基因组分析
# 日期：2026-05-14
# 说明：加载TCGA-BRCA表达和临床数据，执行清洗、标准化、缺失值处理
# ===========================================================================

# 0. 环境设置 ===============================================================
Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE, timeout = 600)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(mice)
  library(VIM)
})

cat("\n========== 模块1：BRCA数据预处理 ==========\n\n")

# 路径定义
DATA_DIR    <- "D:/Users/Desktop/RProject/data"
PRIVATE_DIR <- file.path(DATA_DIR, "private_data")
OUTPUT_DIR  <- file.path(DATA_DIR, "processed")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. 加载私有数据 ----
cat("Step 1: 加载TCGA-BRCA私有数据...\n")

e <- new.env()
load(file.path(PRIVATE_DIR, "brca_exp.Rdata"), envir = e)
counts_raw <- e$brca
rm(e); gc()

e <- new.env()
load(file.path(PRIVATE_DIR, "brca_clinical.Rdata"), envir = e)
clinical_raw <- e$brca_clinical
rm(e); gc()

cat(sprintf("  表达矩阵: %d 基因 x %d 样本\n", nrow(counts_raw), ncol(counts_raw)))
cat(sprintf("  临床数据: %d 患者 x %d 临床变量\n", nrow(clinical_raw), ncol(clinical_raw)))
cat(sprintf("  数据类型: %s (范围 %.0f - %.0f)\n", typeof(counts_raw), min(counts_raw), max(counts_raw)))

# ---- 2. TCGA Barcode解析 ----
cat("\nStep 2: 解析TCGA barcode，区分肿瘤/正常样本...\n")

barcodes    <- colnames(counts_raw)
patient_ids <- substr(barcodes, 9, 12)
sample_code <- substr(barcodes, 14, 15)

sample_info <- data.frame(
  barcode     = barcodes,
  patient_id  = patient_ids,
  sample_code = sample_code,
  sample_type = case_when(
    sample_code == "01" ~ "Tumor",
    sample_code == "06" ~ "Metastatic",
    sample_code == "11" ~ "Normal",
    TRUE                 ~ "Other"
  ),
  stringsAsFactors = FALSE
)

cat(sprintf("  样本分布: %s\n", paste(names(table(sample_info$sample_type)), 
                                   table(sample_info$sample_type), sep="=", collapse=", ")))
write.csv(sample_info, file.path(OUTPUT_DIR, "sample_info.csv"), row.names = FALSE)

# ---- 3. 拆分肿瘤/正常样本矩阵 ----
cat("\nStep 3: 拆分肿瘤/正常表达矩阵...\n")

tumor_idx   <- which(sample_info$sample_type == "Tumor")
normal_idx  <- which(sample_info$sample_type == "Normal")

counts_tumor  <- counts_raw[, tumor_idx, drop = FALSE]
counts_normal <- counts_raw[, normal_idx, drop = FALSE]

cat(sprintf("  肿瘤样本: %d, 正常样本: %d\n", ncol(counts_tumor), ncol(counts_normal)))

# ---- 4. 低表达基因过滤 ----
cat("\nStep 4: 过滤低表达基因...\n")

min_samples <- max(2, ceiling(ncol(counts_tumor) * 0.1))
keep_genes  <- rowSums(counts_tumor >= 10) >= min_samples

counts_tumor  <- counts_tumor[keep_genes, , drop = FALSE]
counts_normal <- counts_normal[keep_genes, , drop = FALSE]
counts_raw    <- counts_raw[keep_genes, , drop = FALSE]

cat(sprintf("  过滤前基因: %d, 过滤后: %d (保留 %.1f%%)\n", 
            length(keep_genes), sum(keep_genes), 100 * mean(keep_genes)))

# ---- 5. ENSEMBL注释下载 (一次性获取基因symbol+位置) ----
cat("\nStep 5: 下载Ensembl基因注释 (symbol + 基因长度)...\n")

ensembl_ids   <- rownames(counts_tumor)
ensembl_clean <- gsub("\\.\\d+$", "", ensembl_ids)
cat(sprintf("  共 %d 个ENSEMBL基因ID\n", length(ensembl_clean)))

gene_map_file <- file.path(OUTPUT_DIR, "ensembl_gene_map.rds")

if (file.exists(gene_map_file)) {
  cat("  从本地缓存加载基因注释...\n")
  gene_map <- readRDS(gene_map_file)
} else {
  ftp_url <- "https://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.113.gene.txt.gz"
  cat(sprintf("  下载Ensembl基因注释: %s\n", ftp_url))
  
  tmp_gz <- tempfile(fileext = ".gz")
  success <- tryCatch({
    download.file(ftp_url, tmp_gz, mode = "wb", timeout = 600)
    TRUE
  }, error = function(e) { FALSE })
  
  if (success) {
    gene_data <- read.table(gzfile(tmp_gz), header = TRUE, sep = "\t",
                            quote = "", comment.char = "", stringsAsFactors = FALSE,
                            fill = TRUE, na.strings = "")
    unlink(tmp_gz)
    
    # 提取: gene_stable_id, gene_name, gene_biotype, description, start, end
    avail_cols <- intersect(c("gene_stable_id", "gene_name", "gene_biotype", 
                               "description", "start_position", "end_position"),
                            colnames(gene_data))
    gene_map <- gene_data[, avail_cols]
    
    # 计算基因长度
    if (all(c("start_position", "end_position") %in% avail_cols)) {
      gene_map$gene_length <- gene_map$end_position - gene_map$start_position + 1
    }
    
    saveRDS(gene_map, gene_map_file)
    cat(sprintf("  下载完成: %d个基因注释\n", nrow(gene_map)))
  } else {
    cat("  Ensembl下载失败，使用备用方案\n")
    gene_map <- data.frame(
      gene_stable_id = ensembl_clean,
      gene_name = ensembl_clean,
      gene_biotype = "unknown",
      description = "unknown",
      gene_length = 1000,
      stringsAsFactors = FALSE
    )
  }
}

# ---- 6. 基因名称和长度赋值 ----
cat("\nStep 6: 基因名称和长度赋值...\n")

# 构建映射：ensembl_clean → symbol
map_idx <- match(ensembl_clean, gene_map$gene_stable_id)

gene_symbol <- ifelse(!is.na(map_idx) & gene_map$gene_name[map_idx] != "",
                      gene_map$gene_name[map_idx], ensembl_ids)
# 处理缺失
gene_symbol[is.na(gene_symbol) | gene_symbol == ""] <- ensembl_ids[is.na(gene_symbol) | gene_symbol == ""]
# 处理重复symbol
dup_idx <- which(duplicated(gene_symbol) | duplicated(gene_symbol, fromLast = TRUE))
gene_symbol[dup_idx] <- paste0(gene_symbol[dup_idx], "_", ensembl_clean[dup_idx])

gene_length_vec <- if ("gene_length" %in% colnames(gene_map) && !is.null(map_idx)) {
  len <- gene_map$gene_length[map_idx]
  len[is.na(len) | len <= 0] <- median(len[!is.na(len) & len > 0], na.rm = TRUE)
  if (all(is.na(len))) len <- rep(1000, length(len))
  len
} else {
  rep(1000, nrow(counts_tumor))
}

cat(sprintf("  有效基因长度: %d (中位数: %.0f bp)\n",
            sum(gene_length_vec > 0), median(gene_length_vec, na.rm = TRUE)))
cat(sprintf("  唯一基因symbol: %d\n", length(unique(gene_symbol))))

# 为矩阵设置行名
rownames(counts_tumor)  <- gene_symbol
rownames(counts_normal) <- gene_symbol
rownames(counts_raw)    <- gene_symbol

# 构建基因注释表
gene_anno <- data.frame(
  ensembl_id      = ensembl_ids,
  ensembl_clean   = ensembl_clean,
  gene_symbol     = gene_symbol,
  gene_length     = gene_length_vec,
  gene_biotype    = ifelse(!is.na(map_idx), gene_map$gene_biotype[map_idx], "unknown"),
  description     = ifelse(!is.na(map_idx), gene_map$description[map_idx], "unknown"),
  stringsAsFactors = FALSE
)

# ---- 7. Counts → TPM 标准化 ----
cat("\nStep 7: Counts → TPM标准化...\n")

counts_to_tpm <- function(counts, gene_lengths) {
  rpk <- counts / (gene_lengths / 1000)
  sweep(rpk, 2, colSums(rpk) / 1e6, "/")
}

tpm_tumor  <- counts_to_tpm(counts_tumor, gene_length_vec)
tpm_normal <- counts_to_tpm(counts_normal, gene_length_vec)

cat(sprintf("  TPM完成 (log2范围: 肿瘤 %.2f-%.2f, 正常 %.2f-%.2f)\n",
            min(log2(tpm_tumor + 1)), max(log2(tpm_tumor + 1)),
            min(log2(tpm_normal + 1)), max(log2(tpm_normal + 1))))

# ---- 8. 临床数据清洗 ----
cat("\nStep 8: 临床数据清洗与变量提取...\n")

# 提取关键临床变量
key_vars <- c("bcr_patient_barcode", "vital_status", "days_to_death",
              "days_to_last_followup", "age_at_initial_pathologic_diagnosis",
              "gender", "race_list", "stage_event_pathologic_stage",
              "stage_event_tnm_categories",
              "breast_carcinoma_estrogen_receptor_status",
              "breast_carcinoma_progesterone_receptor_status",
              "lab_proc_her2_neu_immunohistochemistry_receptor_status",
              "histological_type", "menopause_status",
              "person_neoplasm_cancer_status",
              "number_of_lymphnodes_positive_by_he",
              "radiation_therapy", "postoperative_rx_tx")

clinical <- clinical_raw[, intersect(key_vars, colnames(clinical_raw)), drop = FALSE]

# 重命名
name_map <- c(
  "bcr_patient_barcode" = "patient_id",
  "vital_status" = "vital_status",
  "days_to_death" = "days_to_death",
  "days_to_last_followup" = "days_to_last_followup",
  "age_at_initial_pathologic_diagnosis" = "age_at_diagnosis",
  "gender" = "gender",
  "race_list" = "race",
  "stage_event_pathologic_stage" = "pathologic_stage",
  "stage_event_tnm_categories" = "tnm_categories",
  "breast_carcinoma_estrogen_receptor_status" = "er_status",
  "breast_carcinoma_progesterone_receptor_status" = "pr_status",
  "lab_proc_her2_neu_immunohistochemistry_receptor_status" = "her2_status",
  "histological_type" = "histological_type",
  "menopause_status" = "menopause_status",
  "person_neoplasm_cancer_status" = "cancer_status",
  "number_of_lymphnodes_positive_by_he" = "positive_lymph_nodes",
  "radiation_therapy" = "radiation_therapy",
  "postoperative_rx_tx" = "postoperative_rx"
)
for (nm in intersect(names(name_map), colnames(clinical))) {
  colnames(clinical)[colnames(clinical) == nm] <- name_map[nm]
}

# 计算OS和生存事件
clinical <- clinical %>%
  mutate(
    vital_status = ifelse(vital_status == "Dead", 1, 0),
    os_time = ifelse(vital_status == 1,
                     as.numeric(days_to_death),
                     as.numeric(days_to_last_followup)),
    os_time = ifelse(os_time <= 0 | is.na(os_time), 1, os_time)
  )

# 简化病理分期
clinical <- clinical %>%
  mutate(
    stage_simple = case_when(
      grepl("Stage I[A-C]?$|Stage I$", pathologic_stage, ignore.case = TRUE) ~ "Stage I",
      grepl("Stage II[A-C]?$|Stage II$", pathologic_stage, ignore.case = TRUE) ~ "Stage II",
      grepl("Stage III[A-C]?$|Stage III$", pathologic_stage, ignore.case = TRUE) ~ "Stage III",
      grepl("Stage IV[A-C]?$|Stage IV$", pathologic_stage, ignore.case = TRUE) ~ "Stage IV",
      TRUE ~ NA_character_
    )
  )

# 分子亚型分类
clinical <- clinical %>%
  mutate(
    er_status_clean = case_when(
      grepl("Positive", er_status, ignore.case = TRUE) ~ "Positive",
      grepl("Negative", er_status, ignore.case = TRUE) ~ "Negative",
      TRUE ~ NA_character_
    ),
    pr_status_clean = case_when(
      grepl("Positive", pr_status, ignore.case = TRUE) ~ "Positive",
      grepl("Negative", pr_status, ignore.case = TRUE) ~ "Negative",
      TRUE ~ NA_character_
    ),
    her2_status_clean = case_when(
      grepl("Positive", her2_status, ignore.case = TRUE) ~ "Positive",
      grepl("Negative", her2_status, ignore.case = TRUE) ~ "Negative",
      grepl("Equivocal|Indeterminate", her2_status, ignore.case = TRUE) ~ "Equivocal",
      TRUE ~ NA_character_
    ),
    molecular_subtype = case_when(
      (er_status_clean == "Positive" | pr_status_clean == "Positive") & her2_status_clean == "Negative" ~ "Luminal A",
      (er_status_clean == "Positive" | pr_status_clean == "Positive") & her2_status_clean == "Positive" ~ "Luminal B",
      er_status_clean == "Negative" & pr_status_clean == "Negative" & her2_status_clean == "Positive" ~ "HER2-enriched",
      er_status_clean == "Negative" & pr_status_clean == "Negative" & her2_status_clean == "Negative" ~ "Triple Negative",
      TRUE ~ NA_character_
    )
  )

cat(sprintf("  提取 %d 个关键临床变量\n", ncol(clinical)))

# ---- 9. 缺失值处理 ----
cat("\nStep 9: 缺失值处理...\n")

impute_cols <- c("age_at_diagnosis", "positive_lymph_nodes", "os_time",
                 "stage_simple", "molecular_subtype", "race")

impute_cols <- intersect(impute_cols, colnames(clinical))

for (v in names(clinical)) {
  miss_n <- sum(is.na(clinical[[v]]))
  if (miss_n > 0) cat(sprintf("    %-25s: %d (%.1f%%)\n", v, miss_n, 100*miss_n/nrow(clinical)))
}

# 简单中位数/众数填充
for (v in intersect(impute_cols, colnames(clinical))) {
  nas <- is.na(clinical[[v]])
  if (any(nas)) {
    if (is.numeric(clinical[[v]])) {
      clinical[[v]][nas] <- median(clinical[[v]], na.rm = TRUE)
    } else {
      tbl <- table(clinical[[v]])
      clinical[[v]][nas] <- names(tbl)[which.max(tbl)]
    }
  }
}
cat("  缺失值填充完成\n")

# ---- 10. 临床-表达对齐 ----
cat("\nStep 10: 临床-表达样本对齐...\n")

expr_patients_id <- substr(colnames(counts_tumor), 9, 12)
clinical$patient_id_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)

match_idx <- match(expr_patients_id, clinical$patient_id_short)
matched   <- !is.na(match_idx)

counts_tumor_final  <- counts_tumor[, matched, drop = FALSE]
tpm_tumor_final     <- tpm_tumor[, matched, drop = FALSE]
clinical_final      <- clinical[match_idx[matched], , drop = FALSE]

cat(sprintf("  匹配: %d / %d (%.1f%%)\n",
            ncol(counts_tumor_final), ncol(counts_tumor), 100 * mean(matched)))

# ---- 11. 保存数据 ----
cat("\nStep 11: 保存处理后数据...\n")

saveRDS(counts_tumor_final, file.path(OUTPUT_DIR, "brca_tumor_counts.rds"))
saveRDS(tpm_tumor_final, file.path(OUTPUT_DIR, "brca_tumor_tpm.rds"))
saveRDS(clinical_final, file.path(OUTPUT_DIR, "brca_clinical_clean.rds"))
saveRDS(gene_anno, file.path(OUTPUT_DIR, "gene_annotation.rds"))

if (ncol(counts_normal) > 0) {
  saveRDS(counts_normal, file.path(OUTPUT_DIR, "brca_normal_counts.rds"))
  saveRDS(tpm_normal, file.path(OUTPUT_DIR, "brca_normal_tpm.rds"))
}

# 完整打包
saveRDS(list(
  counts           = counts_tumor_final,
  tpm              = tpm_tumor_final,
  clinical         = clinical_final,
  gene_annotation  = gene_anno,
  sample_info      = sample_info,
  preprocessing_date = Sys.time()
), file.path(OUTPUT_DIR, "brca_tumor_processed.rds"))

# ---- 12. 预处理报告 ----
cat("\n========== 预处理完成报告 ==========\n")
cat(sprintf("  输出目录: %s\n", OUTPUT_DIR))
cat(sprintf("  肿瘤样本: %d (表达矩阵 %d基因 x %d样本)\n",
            ncol(counts_tumor_final), nrow(counts_tumor_final), ncol(counts_tumor_final)))
cat(sprintf("  正常样本: %d\n", ncol(counts_normal)))
cat(sprintf("  临床变量: %d\n", ncol(clinical_final)))

if ("stage_simple" %in% colnames(clinical_final)) {
  cat(sprintf("  病理分期: %s\n",
              paste(names(table(clinical_final$stage_simple)),
                    table(clinical_final$stage_simple), sep="=", collapse=", ")))
}
if ("molecular_subtype" %in% colnames(clinical_final)) {
  cat(sprintf("  分子亚型: %s\n",
              paste(names(table(clinical_final$molecular_subtype)),
                    table(clinical_final$molecular_subtype), sep="=", collapse=", ")))
}
cat(sprintf("  生存: 死亡=%d 存活=%d\n",
            sum(clinical_final$vital_status == 1, na.rm = TRUE),
            sum(clinical_final$vital_status == 0, na.rm = TRUE)))
cat("=========================================\n")

# 摘要表
summary_table <- data.frame(
  Metric = c("原始基因数","过滤后基因数","肿瘤样本数","正常样本数",
             "临床变量数","死亡事件","存活",
             "Stage I","Stage II","Stage III","Stage IV",
             "Luminal A","Luminal B","HER2-enriched","Triple Negative"),
  Value = c(
    length(keep_genes), nrow(counts_tumor_final), ncol(counts_tumor_final),
    ncol(counts_normal), ncol(clinical_final),
    sum(clinical_final$vital_status == 1, na.rm = TRUE),
    sum(clinical_final$vital_status == 0, na.rm = TRUE),
    sum(clinical_final$stage_simple == "Stage I", na.rm = TRUE),
    sum(clinical_final$stage_simple == "Stage II", na.rm = TRUE),
    sum(clinical_final$stage_simple == "Stage III", na.rm = TRUE),
    sum(clinical_final$stage_simple == "Stage IV", na.rm = TRUE),
    sum(clinical_final$molecular_subtype == "Luminal A", na.rm = TRUE),
    sum(clinical_final$molecular_subtype == "Luminal B", na.rm = TRUE),
    sum(clinical_final$molecular_subtype == "HER2-enriched", na.rm = TRUE),
    sum(clinical_final$molecular_subtype == "Triple Negative", na.rm = TRUE)
  )
)
write.csv(summary_table, file.path(OUTPUT_DIR, "preprocessing_summary.csv"), row.names = FALSE)

cat("预处理全部完成!\n")
