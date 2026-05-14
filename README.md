# BRCA Multi-Omics Data Mining

基于 TCGA-BRCA 数据的乳腺癌多组学数据挖掘完整分析项目。

## 项目概述

本项目对 TCGA-BRCA（乳腺浸润性癌）队列进行了系统性的多组学数据挖掘分析，涵盖以下数据层与分析方法：

| 数据层 | 平台 | 特征数 | 样本量 |
|--------|------|--------|--------|
| mRNA 转录组 | Illumina HiSeq RNA-seq | 25,981 基因 | 1,094 例（肿瘤）+ 113 例（正常） |
| miRNA 转录组 | Illumina HiSeq miRNA-seq | 1,881 miRNA | 1,079 例 |
| 体细胞突变 | WES (MuTect2) | 15,413 基因 | 990 例 |
| 临床数据 | TCGA Clinical | 25 核心变量 | 1,174 例 |

**核心整合队列**：1,076 例（mRNA ∩ miRNA），重叠率 98.4%。

## 分析方法

| 模块 | 方法 | 核心结果 |
|------|------|---------|
| 差异表达分析 | DESeq2 | 6,768 DEGs（4,334↑ / 2,434↓） |
| 功能富集分析 | g:Profiler (gprofiler2) | 上调富集 Cell Cycle (p=1.5×10⁻²¹)，下调富集 ECM organization |
| 分子亚型分类 | LASSO / Random Forest / XGBoost | LASSO 89.4% 准确率（35 基因） |
| 聚类分析 | PCA / t-SNE / K-means / Consensus Clustering | PC1 15.3%，最佳 K=2（ER+/ER- 二分） |
| WGCNA | Signed network, soft power=8 | 9 个共表达模块，blue 模块 hub MM=0.959 |
| 生存分析 | KM / Cox / 单变量 Cox 预后筛选 | C-index=0.771，Stage IV HR=8.74 |
| 突变分析 | maftools | PIK3CA (37.3%), TP53 (35.2%) |
| 关联规则挖掘 | arules (Apriori) | 基因→亚型/受体状态关联规则 |
| miRNA-mRNA 网络 | Pearson 相关矩阵 | 2,065 对 (|r|>0.4) |

## 项目结构

```
R_Work/
├── code/                    # 分析脚本（18个R脚本）
│   ├── 01_data_preprocessing.R        # mRNA预处理
│   ├── 02_download_mirna.R            # miRNA下载
│   ├── 02_download_mutations.R        # 突变下载
│   ├── 03_integrate_omics.R           # 多组学整合
│   ├── 04_diff_expression.R           # 差异表达
│   ├── 05_classification.R            # 分类模型
│   ├── 06_clustering.R                # 聚类分析
│   ├── 07_association_analysis.R      # 关联分析
│   ├── 08_wgcna.R                     # WGCNA
│   ├── 11_survival_analysis.R         # 生存分析
│   ├── 12_advanced_algorithms.R       # 突变+网络
│   ├── 13_visualization.R             # 可视化套件
│   ├── 14_generate_paper.R            # 生成docx论文
│   ├── 15_publication_figures.R       # 发表级图表
│   ├── 17_rebuild_all_figures_v2.R    # 图表标准版
│   └── 18_new_figures.R              # 富集/关联/聚类新图
├── data/
│   ├── private_data/                  # 原始TCGA数据(不上传)
│   ├── public_data/                   # 下载的公共数据(不上传)
│   └── processed/                     # 处理后数据(不上传)
├── results/
│   ├── figures_pub/                   # 发表级图表 (16张PNG, 300 DPI)
│   ├── tables/                        # 结果表格 (CSV)
│   ├── enrichment/                    # 富集分析结果
│   ├── association/                   # 关联规则结果
│   ├── consensus_clustering/          # 共识聚类结果
│   └── BRCA_Paper_final.docx          # 最终论文
├── docs/
│   ├── 01_data_background.md          # 数据背景文档
│   ├── BRCA_Project_Report.md         # 项目完整报告
│   └── paper_final_20260514.md        # 论文Markdown源文件
└── .omc/
    └── figure_standards.md            # 图表标准规范
```

## 环境依赖

- **R**: 4.6.0
- **Bioconductor**: 3.23
- **操作系统**: Windows 11

### 核心 R 包

```r
# 数据获取
BiocManager::install("TCGAbiolinks")
BiocManager::install("TCGAbiolinksGUI.data")

# 差异表达与富集
BiocManager::install("DESeq2")
install.packages("gprofiler2")

# 分类与聚类
install.packages(c("caret", "randomForest", "glmnet"))
install.packages("ConsensusClusterPlus")

# 网络分析
install.packages("WGCNA")

# 生存分析
install.packages(c("survival", "survminer"))

# 突变分析
BiocManager::install("maftools")

# 关联规则
install.packages(c("arules", "arulesViz"))

# 可视化
install.packages(c("ggplot2", "ggrepel", "ggsci", "pheatmap", "gridExtra"))
install.packages("showtext")  # 中文字体支持

# 文档生成
install.packages(c("officer", "flextable"))

# 基因注释
BiocManager::install(c("org.Hs.eg.db", "GO.db"))
```

---

## 数据下载指南

> 由于 TCGA 数据受 GDC 访问协议保护且文件较大，本仓库**不包含任何数据文件**。以下为完整的数据获取流程。

### 1. 私有数据（private_data）

`data/private_data/` 目录下的 `.Rdata` 文件为项目自有的 TCGA-BRCA 数据子集，包含以下 12 个文件：

```
brca_exp.Rdata          # BRCA mRNA 表达矩阵 (60,660 genes × 1,226 samples)
brca_clinical.Rdata     # BRCA 临床数据 (1,174 patients × 94 variables)
chol_exp.Rdata          # CHOL 胆管癌表达矩阵
chol_clinical.Rdata     # CHOL 临床数据
coad_exp.Rdata          # COAD 结肠癌表达矩阵
coad_clinical.Rdata     # COAD 临床数据
kirc_exp.Rdata          # KIRC 肾癌表达矩阵
kirc_clinical.Rdata     # KIRC 临床数据
luad_exp.Rdata          # LUAD 肺腺癌表达矩阵
luad_clinical.Rdata     # LUAD 临床数据
stad_exp.Rdata          # STAD 胃癌表达矩阵
stad_clinical.Rdata     # STAD 临床数据
```

这些文件应放置在 `data/private_data/` 目录下，格式为标准 R `.Rdata`，每个文件包含一个数据对象（表达矩阵为 `matrix`，临床数据为 `data.frame`）。

### 2. 公共数据下载（public_data）

以下数据通过 **TCGAbiolinks** R 包从 GDC API 在线下载。

#### 2.1 miRNA 表达数据

```r
library(TCGAbiolinks)

# 查询 BRCA miRNA-seq 数据
mirna_query <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  access        = "open"
)

# 下载
GDCdownload(mirna_query, method = "api", files.per.chunk = 20)

# 导入并提取 read_count 矩阵
mirna_data <- GDCprepare(mirna_query)
# 输出: data.frame, 1,881 miRNAs × columns (含 read_count / RPM / cross-mapped)
# 脚本 02_download_mirna.R 中已处理 data.frame 格式的提取逻辑
```

#### 2.2 体细胞突变数据

```r
library(TCGAbiolinks)

# 查询 BRCA Masked Somatic Mutation (MuTect2 pipeline)
mut_query <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type     = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access        = "open"
)

# 下载 (~32 MB, 992 files, 50 chunks)
GDCdownload(mut_query, method = "api", files.per.chunk = 20)

# 导入
maf_data <- GDCprepare(mut_query)
# 输出: data.frame, ~89,000+ variants
# 使用 maftools::read.maf() 读取为 MAF 对象进行下游分析
```

#### 2.3 DNA 甲基化数据（可选，本项目未使用）

```r
# 注意: BRCA 450K 甲基化数据约为 11.7 GB (895 files, 90 chunks)
# 下载耗时长且内存需求大，不推荐全量下载
# 如需使用，建议通过 GDCquery 添加 sample_type 过滤或使用启动子探针子集

methyl_query <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "DNA Methylation",
  platform      = "Illumina Human Methylation 450",
  data.type     = "Methylation Beta Value",
  access        = "open"
)
```

### 3. 数据预处理

下载完成后，按以下顺序运行预处理脚本：

```r
# Step 1: mRNA 预处理（必须首先运行）
source("code/01_data_preprocessing.R")
# 输出: data/processed/brca_tumor_processed.rds (含 counts, tpm, clinical)

# Step 2: miRNA 下载与提取
source("code/02_download_mirna.R")
# 输出: data/public_data/brca_mirna_counts.rds, brca_mirna_clinical.rds

# Step 3: 突变数据下载
source("code/02_download_mutations.R")
# 输出: data/public_data/brca_mutations.rds

# Step 4: 多组学患者级整合
source("code/03_integrate_omics.R")
# 输出: data/processed/brca_mRNA_aligned.rds, brca_miRNA_aligned.rds
```

### 4. 数据目录最终结构

完成下载和预处理后，`data/` 目录结构如下：

```
data/
├── private_data/           # 手动放置的私有数据（不上传）
│   ├── brca_exp.Rdata
│   ├── brca_clinical.Rdata
│   └── ... (其他癌种)
├── public_data/            # TCGAbiolinks 下载（不上传）
│   ├── brca_mirna_counts.rds       # 1,881 miRNAs × 1,207 samples
│   ├── brca_mirna_clinical.rds     # 样本元数据
│   ├── brca_mirna_se.rds           # SummarizedExperiment 原始对象
│   ├── brca_mutations.rds          # 89,568 variants
│   └── brca_mutations_temp.maf     # maftools 兼容格式
└── processed/              # 预处理后数据（不上传）
    ├── brca_tumor_processed.rds    # 完整肿瘤数据包 (counts + tpm + clinical)
    ├── brca_tumor_counts.rds       # 肿瘤 counts 矩阵
    ├── brca_tumor_tpm.rds          # 肿瘤 TPM 矩阵
    ├── brca_normal_counts.rds      # 正常组织 counts
    ├── brca_normal_tpm.rds         # 正常组织 TPM
    ├── brca_clinical_clean.rds     # 清洗后临床数据
    ├── brca_mRNA_aligned.rds       # 整合后 mRNA (1,076 patients)
    ├── brca_miRNA_aligned.rds      # 整合后 miRNA (1,076 patients)
    ├── brca_multimics_sample_map.csv  # 样本映射表
    ├── gene_annotation.rds         # 基因注释
    └── sample_info.csv             # 样本信息
```

---

## 运行分析

在数据准备完毕后，可单独运行任一分析脚本：

```r
# 差异表达分析
source("code/04_diff_expression.R")

# 分类模型
source("code/05_classification.R")

# 聚类分析
source("code/06_clustering.R")

# 关联分析
source("code/07_association_analysis.R")

# WGCNA（内存需求较高，建议单独运行）
source("code/08_wgcna.R")

# 生存分析
source("code/11_survival_analysis.R")

# 突变+多组学整合
source("code/12_advanced_algorithms.R")

# 可视化套件（含HTML报告）
source("code/13_visualization.R")

# 生成发表级图表
source("code/17_rebuild_all_figures_v2.R")
source("code/18_new_figures.R")

# 生成论文 .docx
source("code/14_generate_paper.R")
```

---

## 图表标准

本项目所有发表级图表遵循统一规范（详见 `.omc/figure_standards.md`）：

- 格式：PNG, 300 DPI, 白色背景
- 语言：英文
- 配色：NPG (Nature Publishing Group) 色板，通过 `ggsci::pal_npg("nrc")`
- 字号：`theme_classic(base_size=15-17)`，标题 18-21pt bold
- 所有图表由 `code/17_rebuild_all_figures_v2.R` 及 `code/18_new_figures.R` 一键生成

## 关键技术决策与踩坑记录

1. **分子亚型分类 Bug**：原始代码中 `ER+ | PR+ & HER2-` 因 R 运算符优先级被解析为 `ER+ | (PR+ & HER2-)`，导致 ER+/HER2+ 患者被误归入 Luminal A。修正为 `(ER+ | PR+) & HER2-`。
2. **甲基化数据放弃**：BRCA Illumina 450K 甲基化数据总量 11.7 GB (895 files)，下载和加载均超出当前环境承受范围。
3. **LASSO-Cox 零结果**：87/1,105 的事件率过低，LASSO-Cox 无法选出稳定的预后基因。改用单变量 Cox 初筛策略。
4. **Ensembl FTP 下载失败**：GRCh38.113 gene.txt.gz 返回 404，fallback 为统一基因长度 1,000 bp。
5. **g:Profiler 替代 clusterProfiler**：clusterProfiler 的 GO.db/org.Hs.eg.db 依赖链存在安装问题，改用 gprofiler2 在线 API 完成 GO/KEGG/Reactome/WikiPathway 富集分析。
6. **ggsurvplot 风险表中文渲染异常**：survminer 的 ggsurvplot 使用 base R 图形设备，showtext 对其无效。改用纯 ggplot2 手绘 KM 曲线 + grid.arrange 拼接风险表。

## 参考文献

- Sung H, et al. Global Cancer Statistics 2022. CA Cancer J Clin, 2024.
- Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. Nature, 2012.
- Perou CM, et al. Molecular portraits of human breast tumours. Nature, 2000.
- Love MI, et al. DESeq2. Genome Biology, 2014.
- Langfelder P, Horvath S. WGCNA. BMC Bioinformatics, 2008.
- Kolberg L, et al. g:Profiler (2023 update). Nucleic Acids Research, 2023.
- Mayakonda A, et al. Maftools. Genome Research, 2018.
