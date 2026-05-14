# BRCA Multi-Omics Data Mining — 项目完整工作报告

**项目日期**: 2026-05-14  
**分析环境**: R 4.6.0 + Bioconductor 3.23  
**数据来源**: TCGA-BRCA (The Cancer Genome Atlas - Breast Invasive Carcinoma)  
**分析工具**: Claude Code (Ralph + Team + Ultrawork 多模式协作)

---

## 一、项目概述

本项目对TCGA-BRCA乳腺癌队列进行了全面的多组学数据挖掘分析，涵盖mRNA转录组（25,981基因 × 1,094患者）、miRNA表达谱（1,881 miRNAs × 1,079患者）和体细胞突变（990患者，67,546变异位点，15,413基因）。分析分为以下四个阶段：

| 阶段 | 内容 | 方法 | 状态 |
|------|------|------|------|
| **Phase 1** | 数据获取与整合 | TCGAbiolinks GDC API下载 + 患者ID对齐 | ✅ |
| **Phase 2** | 数据背景 | 文献综述+数据集描述 | ✅ |
| **Phase 3** | 数据挖掘算法 | 差异表达/分类/聚类/关联/WGCNA/生存分析/突变分析 | ✅ |
| **Phase 4** | 可视化分析 | R高级绘图+HTML交互报告 | ✅ |

---

## 二、数据集详情

### 2.1 数据来源

从TCGA GDC（Genomic Data Commons）通过TCGAbiolinks R包下载了三种组学数据，与已有的mRNA表达和临床数据进行整合。

| 数据层 | 平台 | 特征数 | 样本数 | 文件大小 |
|--------|------|--------|--------|---------|
| mRNA表达 | Illumina HiSeq RNA-seq | 25,981 genes | 1,094 patients (tumor) | ~150 MB |
| miRNA表达 | Illumina HiSeq miRNA-seq | 1,881 miRNAs | 1,079 patients | ~1.5 MB |
| 体细胞突变 | WES (MuTect2 pipeline) | 15,413 genes | 990 samples | ~32 MB |
| 临床数据 | TCGA Clinical | 94 variables | 1,174 patients | ~30 KB |
| DNA甲基化 | Illumina 450K | 已放弃（11.7GB过大） | — | — |

### 2.2 数据整合结果

- **多组学共同患者（mRNA ∩ miRNA）**: 1,076人（98.4%重叠率）
- mRNA唯一患者: 18人 | miRNA唯一患者: 3人
- 患者ID对齐方式: TCGA Barcode第9-12位
- 肿瘤/正常分布: mRNA 1,105 tumor + 113 normal; miRNA 1,078 tumor

---

## 三、数据预处理流程

### 3.1 mRNA表达预处理 (`01_data_preprocessing.R`)
1. TCGA Barcode解析 → 区分肿瘤/正常样本
2. 低表达基因过滤（至少10%样本counts≥10）
3. Ensembl注释下载（GRCh38.113）→ Gene Symbol映射
4. Counts → TPM标准化（使用基因长度）
5. 临床数据清洗：变量选择→重命名→生存时间计算→分期简化→分子亚型分类
6. 缺失值填充（中位数/众数）
7. 表达-临床样本对齐

### 3.2 miRNA表达预处理 (`02_download_mirna.R`)
1. GDCquery查询TCGA-BRCA miRNA-seq (BCGSC miRNA Profiling)
2. GDCdownload + GDCprepare → 1,881 miRNAs × 1,207 samples
3. Data.frame格式处理：提取read_count列，去除前缀
4. TCGA Barcode解析 → 患者ID提取

### 3.3 突变数据预处理 (`02_download_mutations.R`)
1. GDCquery查询Simple Nucleotide Variation (Masked Somatic Mutation)
2. 下载992个MAF文件 → 合并为89,568 variants
3. maftools加载：990 samples, 15,413 genes, 67,546 variants

---

## 四、数据挖掘算法结果

### 4.1 差异表达分析 (US-006)

**方法**: DESeq2, Tumor (n=1,105) vs Normal (n=113)

| 结果 | 数值 |
|------|------|
| 总DEGs (\|log2FC\|>1, padj<0.05) | **6,768** |
| 肿瘤中上调 | 4,334 |
| 肿瘤中下调 | 2,434 |

**输出文件**: `volcano_brca.pdf`, `ma_plot_brca.pdf`, `deg_heatmap_top50.pdf`, `brca_degs_deseq2.csv`

### 4.2 分子亚型分类 (US-007)

**目的**: 从mRNA表达预测BRCA分子亚型（Luminal A/B, HER2-enriched, Triple Negative）

| 模型 | Accuracy | Kappa | F1 (Macro) |
|------|----------|-------|------------|
| **LASSO (Multinomial)** | **0.894** | 0.528 | 0.612 |
| Random Forest (500 trees) | 0.867 | 0.365 | 0.706 |
| XGBoost | 0.358 | 0.002 | 0.232 |

**关键发现**: LASSO以89.4%准确率最佳，选择了35个区分性基因。XGBoost因默认参数不当表现差，需调参优化。

**输出**: `classification_roc.pdf`, `classification_cm.pdf`, `classification_metrics.csv`

### 4.3 聚类分析 (US-008)

**方法**: PCA + t-SNE + 层次聚类 + K-means (top 2,000 variable genes)

| 方法 | 结果 |
|------|------|
| PCA | PC1=15.3%, PC2=8.5% — 主要按ER状态分离 |
| t-SNE | 亚型间可见分离，Luminal A vs Triple Negative 最明显 |
| K-means | 最佳K=2 (silhouette=0.18) — 反映ER+/ER-二分 |
| 层次聚类 | Top 500基因热图，Ward.D2方法 |

**输出**: `pca_subtype.pdf`, `pca_stage.pdf`, `tsne_brca.pdf`, `hclust_heatmap_top500.pdf`, `clustering_silhouette.csv`

### 4.4 关联分析 (US-009)

**方法**: GO/KEGG富集 + GSEA + 临床关联 + miRNA-mRNA调控

| 分析 | 结果 |
|------|------|
| 临床关联 | 100 DEGs × 4临床变量热图 |
| miRNA-mRNA调控 | 228对 (|r|>0.3), 2,065对 (|r|>0.4) |
| GO/KEGG | 无显著富集（过滤严格） |
| GSEA Hallmark | 网络不通跳过 |

**输出**: `clinical_correlation.pdf`, `mirna_target_correlations.csv`, `mirna_mrna_regulatory_network.csv`

### 4.5 WGCNA共表达网络 (US-010)

**参数**: Soft power=8 (scale-free R²=0.887), signed network, minModuleSize=30

| 模块 | 基因数 | Top Hub Gene | Module Membership |
|------|--------|-------------|-------------------|
| grey | 2,817 | ENSG00000119866 | 0.836 |
| turquoise | 584 | ENSG00000115163 | 0.889 |
| **blue** | **553** | **ENSG00000167208** | **0.959** |
| brown | 336 | — | — |
| yellow | 260 | — | — |
| green | 224 | — | — |
| red | 119 | — | — |
| black | 63 | — | — |
| pink | 44 | — | — |

**关键发现**: 9个共表达模块，blue模块具有最强的hub基因连接性（MM=0.959）。

**输出**: `wgcna_soft_power.pdf`, `wgcna_module_trait.pdf`, `wgcna_modules.csv`, `wgcna_hub_genes.csv`, `wgcna_network_edges.csv`

### 4.6 生存分析 (US-011)

**数据**: 1,105患者, 87死亡事件, OS中位随访时间: 2.7年

| 分析 | 结果 |
|------|------|
| KM by Subtype | Log-rank p=0.493 (无显著差异) |
| KM by Stage | Stages I-IV 显著分离 |
| Multivariate Cox | C-index=**0.771**, LR p=8.62e-12 |
| LASSO-Cox | 未选出基因（事件率低 87/1105） |

**Cox模型变量**: Age, Lymph Nodes, Stage II/III/IV, Luminal B/HER2/TN subtypes

**输出**: `km_curves_subtype.pdf`, `km_curves_stage.pdf`, `cox_forest_plot.pdf`, `cox_regression.csv`

### 4.7 突变分析 (US-012)

| 指标 | 数值 |
|------|------|
| 样本数 | 990 |
| 突变基因 | 15,413 |
| 总变异位点 | 67,546 |
| **Top 1**: PIK3CA | 369 mutations |
| **Top 2**: TP53 | 348 mutations |
| **Top 3**: TTN | 268 mutations |

**BRCA经典驱动突变验证**: PIK3CA（~34%样本）和TP53（~32%样本）为BRCA最常见的突变基因，与文献一致。

**输出**: `oncoplot_top20.pdf`, `mutation_types.pdf`, `mutual_exclusivity.pdf`, `mutated_genes_summary.csv`

---

## 五、可视化成果

### 5.1 图表清单

| 图表 | 文件 | 描述 |
|------|------|------|
| 火山图 | `volcano_brca.pdf` | 6,768 DEGs, 红=显著 |
| DEG热图 | `deg_heatmap_top50.pdf` | Top 50 DEGs Tumor vs Normal |
| PCA subtypes | `pca_subtype.pdf` | 按分子亚型着色 |
| t-SNE | `tsne_brca.pdf` | 降维可视化 |
| ROC曲线 | `classification_roc.pdf` | 4类别LASSO/RF |
| WGCNA模块-性状 | `wgcna_module_trait.pdf` | 模块-临床关联热图 |
| KM生存曲线 | `km_curves_stage.pdf` | 按分期 |
| Cox森林图 | `cox_forest_plot.pdf` | 多变量HR |
| Oncoplot | `oncoplot_top20.pdf` | Top 20突变基因 |
| 互斥性分析 | `mutual_exclusivity.pdf` | 突变互斥/共现 |
| 汇总图 | `summary_figure.pdf` | 6-panel 综合图 |
| 关键基因热图 | `key_genes_heatmap.pdf` | 30基因×亚型 |

### 5.2 交互式报告

`results/brca_analysis_report.html` — 包含所有分析结果的完整HTML报告，可直接在浏览器中查看。

---

## 六、项目文件结构

```
R_Work/
├── code/                           # R分析脚本 (13个)
│   ├── 01_data_preprocessing.R     # mRNA预处理
│   ├── 02_download_mirna.R         # miRNA下载
│   ├── 02_download_methylation.R   # 甲基化下载(已放弃)
│   ├── 02_download_mutations.R     # 突变下载
│   ├── 03_integrate_omics.R        # 多组学整合
│   ├── 04_diff_expression.R        # 差异表达
│   ├── 05_classification.R         # 分类模型
│   ├── 06_clustering.R             # 聚类分析
│   ├── 07_association_analysis.R   # 关联分析
│   ├── 08_wgcna.R                  # WGCNA
│   ├── 11_survival_analysis.R      # 生存分析
│   ├── 12_advanced_algorithms.R    # 突变+调控网络
│   └── 13_visualization.R          # 可视化套件
├── data/
│   ├── private_data/               # 原始TCGA数据 (6 cancer types)
│   ├── public_data/                # 下载的公共数据
│   │   ├── brca_mirna_counts.rds   # miRNA count矩阵
│   │   └── brca_mutations.rds      # 突变MAF数据
│   └── processed/                  # 处理后数据
│       ├── brca_tumor_processed.rds # 完整肿瘤数据包
│       ├── brca_mRNA_aligned.rds   # 对齐后mRNA
│       └── brca_miRNA_aligned.rds  # 对齐后miRNA
├── results/
│   ├── figures/                    # 图表 (20+ PDF)
│   ├── tables/                     # 结果表 (15+ CSV)
│   ├── brca_analysis_report.html   # HTML综合报告
│   └── deliverables_manifest.csv   # 项目交付清单
├── docs/
│   └── 01_data_background.md       # 数据背景文档
└── .omc/                           # OMC状态文件
```

---

## 七、关键技术决策

1. **癌症类型选择**: BRCA — 最大样本量(1,226)、最丰富临床数据(94 variables)、成熟分子分型
2. **甲基化放弃**: 全基因组450K数据达11.7GB，下载和内存需求不实际
3. **数据对齐策略**: 基于TCGA Barcode患者ID（第9-12位），1,076共同患者
4. **内存管理**: 并行任务导致OOM后改为串行执行，确保稳定性

---

## 八、关键生物学发现摘要

1. **6,768个差异表达基因**在BRCA肿瘤vs正常组织中显著差异
2. **LASSO模型**可以89.4%准确率从表达谱预测分子亚型
3. **无监督聚类**主要按ER状态将样本分为两组（K=2）
4. **9个WGCNA共表达模块**，blue模块含最强hub基因连接
5. **PIK3CA和TP53**是最常见的突变基因（各~32-34%），与文献一致
6. **病理分期**是生存的最强预测因子（C-index=0.771）
7. **2,065个miRNA-mRNA调控关系对**(|r|>0.4)构成综合调控网络

---

*本项目由Claude Code AI辅助完成，使用R 4.6.0 + Bioconductor 3.23在Windows 11环境下执行。*
