# BRCA Multi-Omics Data Mining

[![R 4.6.0](https://img.shields.io/badge/R-4.6.0-blue)](https://www.r-project.org/)
[![Bioconductor 3.23](https://img.shields.io/badge/Bioconductor-3.23-green)](https://www.bioconductor.org/)
[![LaTeX](https://img.shields.io/badge/LaTeX-MiKTeX%2025.12-orange)](https://miktex.org/)
[![License](https://img.shields.io/badge/License-MIT-lightgrey)](LICENSE)

基于 TCGA-BRCA 数据的乳腺癌多组学数据挖掘完整分析项目。涵盖 mRNA 转录组、miRNA 表达谱和体细胞突变三个组学层次，综合运用差异表达、功能富集、机器学习分类、聚类分析、WGCNA 共表达网络、生存分析、突变分析、TMB 分析和关联规则挖掘等方法。

---

## 核心结果速览

| 分析模块 | 方法 | 关键结果 |
|---------|------|---------|
| 差异表达 | DESeq2 (Tumor vs Normal) | **6,768 DEGs** (4,334↑ / 2,434↓) |
| 亚型间差异 | DESeq2 (TN vs Luminal A) | **4,888 DEGs** |
| 功能富集 | g:Profiler (GO/KEGG/Reactome) | 上调富集 Cell Cycle (p=1.5×10⁻²¹) |
| 分子亚型分类 | LASSO / Random Forest / XGBoost | **LASSO 89.4%** 准确率 (35 基因) |
| 聚类分析 | PCA / t-SNE / K-means / Consensus | PC1 15.3%, K=2 (ER+/ER- 二分) |
| WGCNA | Signed network (power=8) | **9 个模块**, blue 模块 hub MM=0.959 |
| 生存分析 | KM / Cox / Nomogram | **C-index=0.771**, Stage IV HR=8.74 |
| 突变分析 | maftools | **PIK3CA 37.3%**, TP53 35.2% |
| TMB 分析 | 亚型间对比 | TN 最高, Luminal A 最低 |
| 关联规则 | Apriori (arules) | 基因→亚型/受体状态规则 |
| miRNA-mRNA 网络 | Pearson 相关 | 2,065 对 (|r|>0.4) |
| 预后列线图 | rms nomogram | 3 年生存预测 |

---

## 项目结构

```
R_Work/
├── code/                              # R 分析脚本
│   ├── 01_data_preprocessing.R              # mRNA 数据预处理
│   ├── 02_download_mirna.R                  # miRNA 表达下载
│   ├── 02_download_methylation.R            # 甲基化下载 (450K, 未使用)
│   ├── 02_download_mutations.R              # 体细胞突变下载
│   ├── 03_integrate_omics.R                 # 多组学数据整合
│   ├── 04_diff_expression.R                 # DESeq2 差异表达分析
│   ├── 05_classification.R                  # ML 分子亚型分类
│   ├── 06_clustering.R                      # 聚类分析 (PCA/t-SNE/K-means/Consensus)
│   ├── 07_association_analysis.R            # Apriori 关联规则挖掘
│   ├── 08_wgcna.R                           # WGCNA 共表达网络
│   ├── 11_survival_analysis.R               # KM / Cox / 列线图
│   ├── 12_advanced_algorithms.R             # maftools 突变 + miRNA-mRNA 网络
│   ├── 13_visualization.R                   # 可视化套件
│   ├── 14_generate_paper.R                  # .docx 论文生成
│   ├── 15_publication_figures.R             # 发表级图表 (初版)
│   ├── 16_fig3_pure_ggplot.R                # Fig3 KM 曲线 (纯 ggplot)
│   ├── 16_figure_fixes.R                    # 图表修复
│   ├── 16_figure_fixes_v2.R                 # 图表修复 v2
│   ├── 17_rebuild_all_figures.R             # 重建全部图表
│   ├── 17_rebuild_all_figures_v2.R          # 重建全部图表 (标准版, 300 DPI)
│   ├── 18_new_figures.R                     # 富集/关联新图
│   ├── 19_insert_figures_tables.R           # .docx 图文混排
│   ├── 19_insert_figures_tables_v2.R        # .docx 图文混排 v2
│   ├── 19_insert_figures_tables_v3.R        # .docx 图文混排 v3
│   ├── 19_v3.R                              # v3 版论文生成
│   ├── 20_deep_analysis.R                   # 深度分析 (亚型DEG/列线图/TMB)
│   ├── 20b_tmb_fix.R                        # TMB 分析修复
│   ├── explore_data.R                       # 探索性数据分析
│   └── _fix_*.R / _audit_*.R / _verify_*.R  # 审计修复与验证脚本 (7 个)
├── paper/                             # LaTeX 论文
│   ├── references.bib                        # GB/T 7714 参考文献
│   └── v0518-*/                             # 时间戳版本快照 (20 个)
│       ├── BRCA_Paper_*.tex                  # 论文源文件
│       ├── BRCA_Paper_*.pdf                  # 编译输出
│       └── references.bib                    # 各版参考文献
├── results/                           # 分析结果
│   ├── figures_pub/                         # 发表级图表 (18 张 PNG, 300 DPI)
│   ├── tables/                              # 结果表格 (22 个 CSV)
│   ├── enrichment/                          # 富集分析结果 (g:Profiler)
│   ├── association/                         # 关联规则结果
│   ├── consensus_clustering/                # 共识聚类结果
│   ├── deep/                                # 深度分析 (列线图/TMB/亚型DEG)
│   ├── figures/                             # 初步图表
│   ├── BRCA_Paper_Final.docx                # .docx 论文 (含图表)
│   ├── BRCA_Paper.docx                      # .docx 论文
│   ├── BRCA_Paper_with_Figures.docx         # .docx 论文 (图文版)
│   ├── brca_analysis_report.html            # HTML 交互报告
│   └── deliverables_manifest.csv            # 交付清单
├── docs/                              # 文档与草稿
│   ├── paper_final_20260514.md              # 论文 Markdown 定稿
│   ├── BRCA_Project_Report.md               # 项目完整报告
│   ├── 01_data_background.md                # 数据背景说明
│   ├── paper_draft_full.md                  # 论文完整草稿
│   ├── paper_draft_polished.md              # 论文润色稿
│   ├── paper_draft_final.md                 # 论文终稿
│   ├── paper_s1~s4_humanized.md             # 章节人性化稿 (4 个)
│   ├── paper_section_*.md                   # 分章节草稿 (5 个)
│   └── template_*.txt                       # 模板内容
└── .gitignore                         # Git 忽略规则 (排除 data/ / 临时文件)
```

---

## 环境依赖

### R 包

```r
# 核心分析
BiocManager::install(c("TCGAbiolinks", "DESeq2", "maftools", "org.Hs.eg.db"))
install.packages(c("gprofiler2", "caret", "randomForest", "glmnet", "xgboost",
                    "WGCNA", "survival", "survminer", "ConsensusClusterPlus",
                    "arules", "arulesViz", "rms"))

# 可视化
install.packages(c("ggplot2", "ggrepel", "ggsci", "pheatmap", "gridExtra",
                    "showtext", "officer", "flextable"))
```

### LaTeX (编译论文 PDF)

需要 MiKTeX 或 TeX Live，安装后进入对应版本目录编译：

```bash
cd paper/v0518-1821     # 使用最新版本
xelatex BRCA_Paper_0518-1821.tex
biber BRCA_Paper_0518-1821
xelatex BRCA_Paper_0518-1821.tex
xelatex BRCA_Paper_0518-1821.tex
```

或上传 `.tex` + `references.bib` + `results/figures_pub/*.png` 到 [Overleaf](https://www.overleaf.com) 在线编译。

---

## 数据说明

> 由于数据受 GDC 协议保护且文件较大，本仓库不包含数据文件。`data/` 目录已通过 `.gitignore` 排除。

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
query <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
                  data.type = "miRNA Expression Quantification",
                  workflow.type = "BCGSC miRNA Profiling", access = "open")
GDCdownload(query, method = "api", files.per.chunk = 20)
data <- GDCprepare(query)
```

**体细胞突变**:
```r
query <- GDCquery(project = "TCGA-BRCA", data.category = "Simple Nucleotide Variation",
                  data.type = "Masked Somatic Mutation",
                  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(query, method = "api", files.per.chunk = 20)
maf <- GDCprepare(query)
```

---

## 运行分析

预处理完成后，各模块独立可运行：

```r
source("code/04_diff_expression.R")      # 差异表达
source("code/05_classification.R")       # 分类模型
source("code/06_clustering.R")           # 聚类分析
source("code/07_association_analysis.R") # 关联规则
source("code/08_wgcna.R")                # WGCNA
source("code/11_survival_analysis.R")    # 生存分析
source("code/12_advanced_algorithms.R")  # 突变 + miRNA 网络
source("code/20_deep_analysis.R")        # 深度分析 (亚型DEG/列线图/TMB)
```

生成发表级图表和论文：

```r
source("code/17_rebuild_all_figures_v2.R")  # 生成全部 18 张 300 DPI 图表
source("code/19_insert_figures_tables_v3.R") # 生成图文混排 .docx 论文
```

---

## 图表标准

所有发表级图表遵循统一规范 (`code/17_rebuild_all_figures_v2.R`)：

- 格式: PNG, 300 DPI, 白色背景
- 语言: 英文标注
- 配色: NPG (Nature Publishing Group) 色板
- 字号: `theme_classic(base_size = 15-17)`, title 18-21pt bold
- 尺寸: 8-14 英寸宽, 5-10 英寸高

### 图表清单 (18 张)

| 编号 | 文件名 | 内容 |
|------|--------|------|
| Fig1 | `fig1_volcano.png` | Tumor vs Normal 火山图 |
| Fig2 | `fig2_pca_subtype.png` | PCA 亚型分布 |
| Fig3 | `fig3_km_stage.png` | AJCC Stage KM 生存曲线 |
| Fig4 | `fig4_oncoplot.png` | 瀑布图 (top 20 突变基因) |
| Fig5 | `fig5_deg_heatmap.png` | DEG 热图 |
| Fig6 | `fig6_model_comparison.png` | 三模型分类性能对比 |
| Fig7 | `fig7_wgcna_sft.png` | WGCNA 软阈值选择 |
| Fig8 | `fig8_wgcna_modules.png` | WGCNA 模块聚类 |
| Fig9 | `fig9_cox_forest.png` | 多变量 Cox 森林图 |
| Fig10 | `fig10_mirna_mrna_corr.png` | miRNA-mRNA 相关性 |
| Fig11 | `fig11_summary.png` | 多组学总览 |
| Fig12 | `fig12_enrichment_up.png` | 上调 DEG 富集 |
| Fig13 | `fig13_enrichment_down.png` | 下调 DEG 富集 |
| Fig14 | `fig14_enrichment_comparison.png` | 富集对比 |
| Fig15 | `fig15_association_rules.png` | 关联规则网络 |
| Fig16 | `fig16_enrichment_top10.png` | Top 10 富集条目 |
| Fig17 | `fig17_tn_vs_la_volcano.png` | TN vs Luminal A 火山图 |
| Fig18 | `fig18_tmb_by_subtype.png` | TMB 亚型箱线图 |

---

## 踩坑记录

1. **分子亚型 Bug**: `ER+ | PR+ & HER2-` 因运算符优先级被误解析，修正为 `(ER+ | PR+) & HER2-`
2. **甲基化数据**: 450K 芯片 11.7 GB，放弃下载，后续可用启动子探针子集
3. **LASSO-Cox 零结果**: 87 例死亡 / 1,105 例患者，改为单变量 Cox 初筛
4. **clusterProfiler 依赖**: GO.db / org.Hs.eg.db 安装链失败，改用 gprofiler2 在线 API
5. **ggsurvplot 中文**: base R 设备不支持 showtext，改用纯 ggplot2 手绘 KM 曲线
6. **.docx 中文乱码**: `Sys.setlocale("LC_ALL", "English")` 破坏编码，去掉即可
7. **LaTeX 图浮动**: `[htbp]` 导致图堆末尾，改为 `[H]` (float 包) 强制就地
8. **TMB 亚型方向**: 初始版本 TN 最低，经修复后确认为 TN 最高 (符合预期)

---

## 参考文献

- Sung H, et al. Global Cancer Statistics 2022. *CA Cancer J Clin*, 2024.
- TCGA Network. Comprehensive molecular portraits of human breast tumours. *Nature*, 2012.
- Perou CM, et al. Molecular portraits of human breast tumours. *Nature*, 2000.
- Love MI, et al. DESeq2. *Genome Biology*, 2014.
- Langfelder P, Horvath S. WGCNA. *BMC Bioinformatics*, 2008.
- Kolberg L, et al. g:Profiler 2023 update. *Nucleic Acids Research*, 2023.
- Mayakonda A, et al. Maftools. *Genome Research*, 2018.
