# BRCA Multi-Omics Data Mining: Data Background and Context

## 1. 癌症的全球疾病负担与乳腺癌的流行病学

### 1.1 癌症的全球疾病负担

癌症是全球第二大死亡原因，据世界卫生组织国际癌症研究机构（IARC）GLOBOCAN 2022统计，全球每年新发癌症病例约2000万例，癌症相关死亡约970万例。癌症不仅造成巨大的生命损失，还带来沉重的社会经济负担——包括直接医疗成本、生产力损失和护理负担。预计到2050年，全球癌症负担将增至每年3500万新发病例，增长主要来自人口老龄化和生活方式的改变。

癌症的发病机制极为复杂，涉及基因组学、表观遗传学、转录组学、蛋白质组学和代谢组学等多个层面的异常调控。Hallmarks of Cancer（Hanahan & Weinberg, 2011）概括了癌症的十大特征：持续的增殖信号、逃避生长抑制、抵抗细胞死亡、无限复制能力、诱导血管生成、激活侵袭和转移、基因组不稳定和突变、肿瘤促进炎症、代谢重编程以及免疫逃逸。这些特征背后的分子机制构成了癌症研究的核心内容。

### 1.2 乳腺癌的流行病学与临床特征

乳腺癌（Breast Invasive Carcinoma, BRCA）是全球女性中最常见的恶性肿瘤，也是癌症相关死亡的主要原因之一。2022年全球新发乳腺癌约230万例，占所有癌症的11.7%，死亡约67万例。

乳腺癌的重要临床和分子特征包括：

- **分子亚型**（基于免疫组化标志物）：
  - **Luminal A型**（ER+/PR+/HER2-/Ki-67低）：约占40-50%，预后最好
  - **Luminal B型**（ER+/PR+/HER2±/Ki-67高）：约占20-30%
  - **HER2富集型**（ER-/PR-/HER2+）：约占10-15%，靶向治疗有效
  - **三阴性乳腺癌/TNBC**（ER-/PR-/HER2-）：约占15-20%，预后最差

- **TNM分期系统**：基于原发肿瘤大小（T）、淋巴结转移（N）和远处转移（M）
- **组织学分级**：Nottingham分级系统（1-3级）

### 1.3 乳腺癌发病的分子机制复杂性

乳腺癌的发生发展涉及多条关键信号通路的异常：

1. **激素受体通路**：雌激素受体（ER）和孕激素受体（PR）通路异常激活驱动约70%的乳腺癌
2. **HER2/ERBB2信号**：约15-20%的乳腺癌存在HER2基因扩增，激活PI3K/AKT/mTOR和MAPK通路
3. **PIK3CA突变**：约30-40%的ER+乳腺癌携带PIK3CA激活突变
4. **TP53突变**：尤其在TNBC中高频发生（~80%）
5. **BRCA1/2突变**：同源重组修复缺陷（HRD），与乳腺癌遗传易感性相关
6. **细胞周期调控异常**：CCND1扩增、CDKN2A缺失
7. **表观遗传改变**：启动子甲基化导致的基因沉默

## 2. 数据挖掘在癌症研究中的重要性

随着高通量测序技术的发展，癌症研究进入了"大数据"时代。TCGA（The Cancer Genome Atlas）、ICGC（International Cancer Genome Consortium）等大型项目产生了海量的多组学数据。数据挖掘技术在这些数据的分析和解读中扮演着不可替代的角色：

### 2.1 关键应用领域

1. **分子分型与精准医学**：通过聚类分析、分类算法识别具有临床意义的分子亚型
2. **生物标志物发现**：通过差异表达分析、生存分析筛选诊断/预后/预测性生物标志物
3. **药物靶点识别**：通过WGCNA等网络分析方法识别关键驱动基因
4. **多组学整合**：整合基因组、转录组、表观组数据揭示疾病的系统生物学机制
5. **预后预测模型**：基于机器学习构建临床结果预测模型

### 2.2 常用数据挖掘方法

- **差异表达分析**（DESeq2, edgeR, limma）：识别实验条件间的差异基因
- **聚类分析**（层次聚类、共识聚类、t-SNE/UMAP）：发现样本的自然分组
- **分类算法**（随机森林、LASSO、XGBoost）：预测样本的分子亚型或临床结果
- **生存分析**（Kaplan-Meier、Cox回归）：评估基因表达与预后的关联
- **WGCNA**（加权基因共表达网络分析）：识别共表达模块和枢纽基因
- **通路富集分析**（GSEA、GO、KEGG）：揭示生物过程层面的调控变化
- **多组学因子分析**（MOFA）：跨组学维度提取潜在因子

## 3. TCGA-BRCA数据集概述

The Cancer Genome Atlas（TCGA）是由美国国家癌症研究所（NCI）和国家人类基因组研究所（NHGRI）联合发起的里程碑式癌症基因组项目。TCGA-BRCA队列包含约1,100例乳腺癌患者的全面多组学数据。

### 3.1 本项目中使用的数据集

| 数据类型 | 来源 | 维度 | 平台 | 说明 |
|---------|------|-----|------|------|
| **mRNA表达** | TCGA-BRCA | 60,660基因 × 1,226样本 | Illumina HiSeq RNA-seq | Counts矩阵，含肿瘤和正常组织 |
| **miRNA表达** | TCGA-BRCA | ~1,881 miRNA × ~1,000样本 | Illumina HiSeq miRNA-seq | BCGSC miRNA Profiling |
| **DNA甲基化** | TCGA-BRCA | ~450,000 CpG位点 × ~800样本 | Illumina 450K甲基化芯片 | Beta值矩阵 |
| **体细胞突变** | TCGA-BRCA | 基因/样本级别 | WES (MuTect2) | MAF格式，突变注释 |
| **临床数据** | TCGA-BRCA | 1,174患者 × 94变量 | - | 生存、分期、受体状态等 |

### 3.2 TCGA Barcode结构

TCGA样本条码遵循标准格式：`TCGA-XX-YYYY-ZZ-WW-XXX-YY`

- `TCGA`：项目标识
- `XX`：组织来源地点（Tissue Source Site）
- `YYYY`：参与者编号（Participant ID）
- `ZZ`：样本类型代码（01=原发肿瘤, 06=转移, 11=正常组织）
- `WW`：样本顺序（Vial）
- `XXX`：分析物（Portion/Analyte）
- `YY`：检测平台（Plate/Center）

### 3.3 数据预处理流程

1. **基因过滤**：去除低表达基因（至少在10%样本中counts ≥ 10）
2. **标准化**：Counts → TPM（Transcripts Per Million）转换，使用Ensembl基因长度
3. **基因注释**：ENSEMBL ID → Gene Symbol映射，处理重复和缺失
4. **缺失值处理**：中位数填充（数值变量）或众数填充（分类变量）
5. **样本对齐**：表达矩阵与临床数据按患者ID匹配

## 4. 本分析的整体方案

### 4.1 分析目标

1. 识别BRCA肿瘤与正常组织间的差异表达基因
2. 构建分子亚型分类模型，筛选亚型特异性生物标志物
3. 通过共识聚类发现新的分子分层方案
4. 识别与临床特征相关的共表达模块
5. 构建和验证多基因预后风险评分
6. 多组学整合揭示BRCA的跨组学调控网络

### 4.2 技术路线

```
数据下载（miRNA + 甲基化 + 突变）→ 多组学整合对齐 → 差异表达分析
→ 分类模型（RF/LASSO/XGBoost）→ 聚类分析（共识聚类/UMAP）
→ 关联分析（GSEA/通路富集）→ WGCNA共表达网络
→ 生存分析（KM/Cox）→ 突变特征分析 → 多组学可视化
```

### 4.3 工具与软件

- **R 4.6.0**：主要分析环境
- **Bioconductor 3.22**：TCGAbiolinks, DESeq2, edgeR, limma, clusterProfiler, maftools, ELMER
- **主要R包**：survival, survminer, WGCNA, caret, randomForest, glmnet, xgboost, ggplot2, ComplexHeatmap, pheatmap

---

*参考文献*
- Hanahan D, Weinberg RA. Hallmarks of Cancer: The Next Generation. Cell. 2011;144(5):646-674.
- Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. Nature. 2012;490(7418):61-70.
- Sung H, et al. Global Cancer Statistics 2022. CA Cancer J Clin. 2024;74(3):229-263.
- Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics. 2008;9:559.
