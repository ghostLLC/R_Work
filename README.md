# BRCA Multi-Omics Data Mining

[![R 4.6.0](https://img.shields.io/badge/R-4.6.0-blue)](https://www.r-project.org/)
[![Bioconductor 3.23](https://img.shields.io/badge/Bioconductor-3.23-green)](https://www.bioconductor.org/)
[![LaTeX](https://img.shields.io/badge/LaTeX-MiKTeX%2025.12-orange)](https://miktex.org/)
[![License](https://img.shields.io/badge/License-MIT-lightgrey)](LICENSE)

基于 TCGA-BRCA 数据的乳腺癌多组学数据挖掘完整分析项目。涵盖 mRNA 转录组、miRNA 表达谱和体细胞突变三个组学层次，综合运用差异表达、功能富集、机器学习分类、聚类、WGCNA 共表达网络、生存分析、突变分析和关联规则挖掘等方法。

---

## 核心结果速览

| 分析模块 | 方法 | 关键结果 |
|---------|------|---------|
| 差异表达 | DESeq2 (Tumor vs Normal) | **6,768 DEGs** (4,334↑ / 2,434↓) |
| 亚型间差异 | DESeq2 (TN vs Luminal A) | **4,888 DEGs** (2,270↑ in TN / 2,618↑ in LA) |
| 功能富集 | g:Profiler (GO/KEGG/Reactome) | 上调富集 Cell Cycle (p=1.5×10⁻²¹) |
| 分子亚型分类 | LASSO / Random Forest / XGBoost | **LASSO 89.4%** 准确率 (35 基因) |
| 聚类 | PCA / t-SNE / K-means / Consensus | PC1 15.3%, K=2 (ER+/ER- 二分) |
| WGCNA | Signed network (power=8) | **9 个模块**, blue 模块 hub MM=0.959 |
| 生存分析 | KM / Cox / Nomogram | **C-index=0.771**, Stage IV HR=8.74 |
| 突变分析 | maftools | **PIK3CA 37.3%**, TP53 35.2% |
| TMB 分析 | 亚型间对比 | TN 最高, Luminal A 最低 |
| 关联规则 | Apriori (arules) | 基因→亚型/受体状态规则 |
| miRNA 网络 | Pearson 相关 | 2,065 对 (|r|>0.4) |
| 预后列线图 | rms nomogram | 3 年生存预测 |

---

## 项目结构

```
R_Work/
├── code/                          # R 分析脚本
│   ├── 01_data_preprocessing.R          # mRNA 预处理
│   ├── 02_download_mirna.R              # miRNA 下载
│   ├── 02_download_mutations.R          # 突变下载
│   ├── 03_integrate_omics.R             # 多组学整合
│   ├── 04_diff_expression.R             # DESeq2 差异表达
│   ├── 05_classification.R              # 分子亚型分类
│   ├── 06_clustering.R                  # 聚类分析
│   ├── 07_association_analysis.R        # 关联分析
│   ├── 08_wgcna.R                       # WGCNA
│   ├── 11_survival_analysis.R           # 生存分析
│   ├── 12_advanced_algorithms.R         # 突变+网络
│   ├── 13_visualization.R               # 可视化套件
│   ├── 14_generate_paper.R              # 生成 .docx
│   ├── 15_publication_figures.R         # 发表级图表 (初版)
│   ├── 17_rebuild_all_figures_v2.R      # 发表级图表 (标准版)
│   ├── 18_new_figures.R                # 富集/关联新图
│   ├── 19_insert_figures_tables_v2.R    # .docx 图文混排
│   └── 20_deep_analysis.R              # 深度分析 (亚型DEG/列线图/TMB)
├── paper/                         # LaTeX 论文
│   ├── BRCA_Paper.tex                   # 主文件 (xelatex 编译)
│   ├── references.bib                    # GB/T 7714 参考文献
│   └── BRCA_Paper.pdf                   # 编译输出 (18 页)
├── results/
│   ├── figures_pub/                     # 发表级图表 (18 张 PNG, 300 DPI)
│   ├── tables/                          # 结果表格 (CSV)
│   ├── enrichment/                      # 富集分析结果
│   ├── association/                     # 关联规则结果
│   ├── consensus_clustering/            # 共识聚类结果
│   ├── deep/                            # 深度分析结果
│   ├── BRCA_Paper_Final.docx            # .docx 论文
│   └── brca_analysis_report.html        # HTML 交互报告
├── docs/                           # 文档
│   ├── paper_final_20260514.md          # 论文 Markdown 源文件
│   └── BRCA_Project_Report.md           # 项目完整报告
└── .omc/
    └── figure_standards.md              # 图表标准规范
```

---

## 环境依赖

### R 包

```r
# 核心分析
BiocManager::install(c("TCGAbiolinks", "DESeq2", "maftools", "org.Hs.eg.db"))
install.packages(c("gprofiler2", "caret", "randomForest", "glmnet", "WGCNA",
                    "survival", "survminer", "ConsensusClusterPlus",
                    "arules", "arulesViz", "rms"))

# 可视化
install.packages(c("ggplot2", "ggrepel", "ggsci", "pheatmap", "gridExtra",
                    "showtext", "officer", "flextable"))
```

### LaTeX (编译论文 PDF)

需要 MiKTeX 或 TeX Live，安装后执行：

```bash
cd paper
xelatex BRCA_Paper.tex
biber BRCA_Paper
xelatex BRCA_Paper.tex
xelatex BRCA_Paper.tex
```

或上传 `BRCA_Paper.tex` + `references.bib` + `results/figures_pub/*.png` 到 [Overleaf](https://www.overleaf.com) 在线编译。

---

## 数据下载

> 由于数据受 GDC 协议保护且文件较大，本仓库不包含数据文件。

### 私有数据

`data/private_data/` 需手动放置以下 `.Rdata` 文件（TCGA-BRCA 表达和临床数据子集）：

```
brca_exp.Rdata          # mRNA 表达矩阵
brca_clinical.Rdata     # 临床数据
```

### 公共数据下载 (TCGAbiolinks)

**miRNA 表达**:
```r
library(TCGAbiolinks)
query <- GDCquery(project="TCGA-BRCA", data.category="Transcriptome Profiling",
                   data.type="miRNA Expression Quantification",
                   workflow.type="BCGSC miRNA Profiling", access="open")
GDCdownload(query, method="api", files.per.chunk=20)
data <- GDCprepare(query)
```

**体细胞突变**:
```r
query <- GDCquery(project="TCGA-BRCA", data.category="Simple Nucleotide Variation",
                   data.type="Masked Somatic Mutation",
                   workflow.type="Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(query, method="api", files.per.chunk=20)
maf <- GDCprepare(query)
```

---

## 运行分析

预处理完成后，各模块独立可运行：

```r
source("code/04_diff_expression.R")     # 差异表达
source("code/05_classification.R")      # 分类模型
source("code/06_clustering.R")          # 聚类
source("code/07_association_analysis.R") # 关联分析
source("code/08_wgcna.R")              # WGCNA
source("code/11_survival_analysis.R")   # 生存分析
source("code/12_advanced_algorithms.R") # 突变+网络
source("code/20_deep_analysis.R")       # 深度分析
```

---

## 图表标准

所有发表级图表遵循统一规范 (`code/17_rebuild_all_figures_v2.R`)：

- 格式: PNG, 300 DPI, 白色背景
- 语言: 英文标注
- 配色: NPG (Nature Publishing Group) 色板
- 字号: `theme_classic(base_size=15-17)`, title 18-21pt bold
- 尺寸: 8-14英寸宽, 5-10英寸高

---

## 踩坑记录

1. **分子亚型 Bug**: `ER+ | PR+ & HER2-` 因运算符优先级被误解析，修正为 `(ER+ | PR+) & HER2-`
2. **甲基化数据**: 450K 芯片 11.7 GB，放弃下载，后续可用启动子探针子集
3. **LASSO-Cox 零结果**: 87 例死亡/1,105 例患者，改为单变量 Cox 初筛
4. **clusterProfiler 依赖**: GO.db/org.Hs.eg.db 安装链失败，改用 gprofiler2 在线 API
5. **ggsurvplot 中文**: base R 设备不支持 showtext，改用纯 ggplot2 手绘
6. **.docx 中文乱码**: `Sys.setlocale("LC_ALL","English")` 破坏了编码，去掉即可
7. **LaTeX 图浮动**: `[htbp]` 导致图堆末尾，改为 `[H]` (float 包) 强制就地

---

## 参考文献

- Sung H, et al. Global Cancer Statistics 2022. *CA Cancer J Clin*, 2024.
- TCGA Network. Comprehensive molecular portraits of human breast tumours. *Nature*, 2012.
- Perou CM, et al. Molecular portraits of human breast tumours. *Nature*, 2000.
- Love MI, et al. DESeq2. *Genome Biology*, 2014.
- Langfelder P, Horvath S. WGCNA. *BMC Bioinformatics*, 2008.
- Kolberg L, et al. g:Profiler 2023 update. *Nucleic Acids Research*, 2023.
- Mayakonda A, et al. Maftools. *Genome Research*, 2018.
