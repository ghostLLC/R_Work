# ===========================================================================
# 19: Insert figures + tables into paper .docx
# ===========================================================================

options(stringsAsFactors = FALSE, encoding = "UTF-8")

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
  library(dplyr)
})

cat("\n========== Inserting Figures & Tables ==========\n\n")

FIG_DIR  <- "D:/Users/Desktop/R_Work/results/figures_pub"
TBL_DIR  <- "D:/Users/Desktop/R_Work/results/tables"
OUTPUT   <- "D:/Users/Desktop/R_Work/results/BRCA_Paper_with_Figures.docx"
DRAFT    <- "D:/Users/Desktop/R_Work/docs/paper_final_20260514.md"

dir.create(dirname(OUTPUT), showWarnings = FALSE, recursive = TRUE)

# =============================================================
# Helper functions
# =============================================================

add_heading <- function(doc, text, level = 1) {
  body_add_par(doc, text, style = paste0("heading ", level))
}

add_text <- function(doc, text) {
  lines <- strsplit(text, "\n")[[1]]
  for (line in lines) {
    line <- trimws(line)
    if (nchar(line) < 2 && !grepl("^$", line)) {
      doc <- body_add_par(doc, "", style = "Normal")
    } else if (nchar(line) >= 2) {
      doc <- body_add_par(doc, line, style = "Normal")
    }
  }
  doc
}

add_fig <- function(doc, path, caption, width = 5.5, height = 4) {
  if (file.exists(path)) {
    doc <- body_add_par(doc, "", style = "Normal")
    doc <- body_add_img(doc, src = path, width = width, height = height)
    doc <- body_add_par(doc, caption, style = "Normal")
    cat(sprintf("  [FIG] %s\n", basename(path)))
  } else {
    cat(sprintf("  [MISS] %s\n", basename(path)))
  }
  doc
}

add_table <- function(doc, df, caption, col_widths = NULL) {
  ft <- flextable(df) %>%
    theme_booktabs() %>%
    fontsize(size = 9, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    align(align = "center", part = "all") %>%
    autofit()
  if (!is.null(col_widths)) {
    ft <- width(ft, width = col_widths)
  }
  doc <- body_add_par(doc, "", style = "Normal")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, caption, style = "Normal")
  cat(sprintf("  [TBL] %s\n", caption))
  doc
}

# =============================================================
# Build document
# =============================================================
cat("Creating document...\n")
doc <- read_docx()

# ---- 摘要 ----
cat("\n--- Abstract ---\n")
doc <- add_heading(doc, "基于多组学数据挖掘的乳腺癌分子特征分析", 1)
doc <- add_heading(doc, "摘要", 2)

abstract_text <- "乳腺癌是全球发病率最高的恶性肿瘤，其分子异质性是精准诊疗面临的核心挑战。本研究基于TCGA-BRCA多组学数据（mRNA转录组25,981个基因、miRNA表达谱1,881个、体细胞突变15,413个基因及临床注释94项），对1,076例多组学共同患者进行了系统性的数据挖掘分析。研究综合运用DESeq2差异表达分析、LASSO多分类回归、随机森林分类、加权基因共表达网络分析（WGCNA）、Cox比例风险回归、maftools突变景观分析、g:Profiler功能富集分析、arules关联规则挖掘及ConsensusClusterPlus共识聚类等算法，构建了从数据预处理到生物学发现的全流程分析管线。主要结果包括：鉴定6,768个显著差异表达基因（|log2FC|>1, FDR<0.05）；LASSO分类器以89.4%准确率实现三分类分子亚型预测；WGCNA识别9个共表达模块，blue模块枢纽基因模块隶属度达0.959；多变量Cox模型C-index为0.771；PIK3CA（37.3%）和TP53（35.2%）为最高频突变基因。功能富集分析揭示上调基因显著富集于细胞周期和DNA复制通路。本研究为乳腺癌多组学数据挖掘提供了完整的分析框架和可复现的计算流程。"
doc <- add_text(doc, abstract_text)
doc <- add_text(doc, "关键词：乳腺癌；TCGA；多组学数据挖掘；差异表达分析；WGCNA；预后模型；功能富集")
doc <- body_add_break(doc)

# ---- 1. 引言 ----
cat("\n--- Section 1: Introduction ---\n")
doc <- add_heading(doc, "1 引言", 2)
doc <- add_heading(doc, "1.1 研究背景", 3)
doc <- add_text(doc, "据国际癌症研究机构（IARC）GLOBOCAN 2022统计，全球每年新发癌症约2,000万例，死亡约970万例。其中，女性乳腺癌以约230万例新发病例居全球恶性肿瘤发病率首位，同年导致约67万例死亡[1]。在中国，乳腺癌的发病率和死亡率持续上升，且发病年龄显著早于欧美国家，构成重大的公共卫生负担。")
doc <- add_text(doc, "乳腺癌的发生发展涉及多基因突变累积、表观遗传重塑及信号通路交互失调等复杂的分子过程[2]。传统单一维度研究方法难以全面揭示乳腺癌的系统性分子特征。高通量测序技术的快速发展及成本的持续下降，使大规模、多层次癌症基因组数据的积累成为可能。以癌症基因组图谱（The Cancer Genome Atlas, TCGA）为代表的国际大型合作项目，为研究者提供了涵盖基因组、转录组、表观基因组及蛋白质组的海量公开数据资源[3]。")
doc <- add_text(doc, "在数据驱动的研究范式下，如何从高维度、高噪声的癌症多组学数据中高效提取具有生物学意义和临床价值的信息，是生物信息学领域的核心问题。将差异表达筛选、网络建模、正则化回归、聚类分析及关联规则挖掘等方法系统应用于乳腺癌多组学数据的整合分析，不仅有助于识别新的分子标志物和预后预测因子，也可为药物靶点发现和个体化治疗方案的制定提供数据驱动的视角。")

doc <- add_heading(doc, "1.2 乳腺癌分子分型", 3)
doc <- add_text(doc, "乳腺癌是由多种分子特征各异的亚型构成的异质性疾病集合。Perou等[4]于2000年首次基于基因表达谱提出乳腺癌分子分型体系，目前国际广泛认可的固有分子亚型包括四类：Luminal A型（约占50%-60%，ER+/PR+/HER2-/Ki-67低，预后最佳）；Luminal B型（约占15%-20%，ER+但Ki-67较高）；HER2过表达型（约占10%-15%，以ERBB2基因扩增为驱动特征）；三阴性型（约占15%-20%，ER-/PR-/HER2-，侵袭性最强）。尽管四分类体系已广泛应用于临床实践，同一亚型内部的生物学异质性仍然显著，整合多组学数据进一步挖掘分子分型的内在结构具有重要的研究价值。")

doc <- add_heading(doc, "1.3 研究内容与组织结构", 3)
doc <- add_text(doc, "本研究以TCGA-BRCA乳腺癌多组学数据为研究对象，综合运用多种数据挖掘方法，系统开展乳腺癌分子特征的整合分析。论文组织结构如下：第2章介绍数据来源与预处理流程；第3章系统阐述所采用的各项数据挖掘算法及分析结果；第4章展示多维度可视化分析；最后为结论与展望。")
doc <- body_add_break(doc)

# ---- 2. 数据来源与预处理 ----
cat("\n--- Section 2: Data ---\n")
doc <- add_heading(doc, "2 数据来源与预处理", 2)

doc <- add_heading(doc, "2.1 TCGA-BRCA数据概述", 3)
doc <- add_text(doc, "TCGA由美国NCI和NHGRI于2006年联合启动，覆盖33种癌症类型、超过20,000例肿瘤样本[3]。乳腺浸润性癌（BRCA）是样本量最大的实体瘤项目。本研究从Genomic Data Commons（GDC）获取Level 3处理级别的数据，涵盖四个核心数据层：mRNA转录组（HiSeq RNA-seq）、miRNA转录组（HiSeq miRNA-seq）、体细胞突变（WES, MuTect2）及临床信息（TCGA Clinical）。TCGA样本条形码前12个字符构成患者级唯一标识符，是多组学数据集成的核心锚点。")

# Table 2.1: Data overview
tbl_data_overview <- data.frame(
  Data_Layer = c("mRNA Expression", "miRNA Expression", "Somatic Mutation", "Clinical Data"),
  Platform = c("HiSeq RNA-seq (STAR-Counts)", "HiSeq miRNA-seq", "WES (MuTect2)", "TCGA Clinical"),
  Features = c("25,981 genes", "1,881 miRNAs", "15,413 genes", "25 core variables"),
  Samples = c("1,094 tumor + 113 normal", "1,079 patients", "990 samples", "1,174 patients"),
  Normalization = c("TPM", "—", "—", "—"),
  stringsAsFactors = FALSE
)
doc <- add_table(doc, tbl_data_overview, "表2.1 多组学数据集总览")

doc <- add_heading(doc, "2.2 数据获取与整合", 3)
doc <- add_text(doc, "所有多组学数据通过TCGAbiolinks R包[11]从GDC API统一获取（GDC Data Release 42.0）。数据获取遵循查询-下载-导入三阶段流水线模式。多组学整合采用患者级对齐策略：提取各数据层样本条形码的前12个字符作为患者ID，取交集确定核心患者队列。整合结果：mRNA覆盖1,094例肿瘤（另含113例正常组织），miRNA覆盖1,079例，两者交集1,076例，重叠率98.4%。")

# Table 2.2: Integration
tbl_integration <- data.frame(
  Metric = c("mRNA-miRNA common patients", "Overlap rate", "Patients with mutation data", "Patients with complete clinical annotation"),
  Value = c("1,076", "98.4%", "983 (91.4%)", "1,076 (100%)"),
  stringsAsFactors = FALSE
)
doc <- add_table(doc, tbl_integration, "表2.2 多组学整合队列特征")

doc <- add_heading(doc, "2.3 数据预处理", 3)
doc <- add_text(doc, "mRNA低表达过滤：原始计数矩阵包含60,660个基因，过滤标准为至少在10%样本中counts>=10，保留25,981个基因（42.8%）。TPM归一化通过Counts Per Million（CPM）校正测序深度，再除以基因有效长度进行基因长度校正。基因标识符通过org.Hs.eg.db建立Ensembl Gene ID到Gene Symbol的映射。临床数据从94个原始变量中提取25个核心变量。分子亚型分类中修正了逻辑运算符优先级导致的分类不一致问题，修正后亚型分布为：Luminal A 950例（86.3%）、HER2-enriched 37例（3.4%）、Triple Negative 115例（10.4%）。分期分布：Stage I 183例（16.6%）、Stage II 651例（58.9%）、Stage III 251例（22.7%）、Stage IV 20例（1.8%）。中位总生存随访时间约2.7年，死亡事件87例（7.9%）。缺失值处理：数值型采用中位数填补，分类型采用众数填补。")

# Table 2.3: Clinical features
tbl_clinical <- data.frame(
  Feature = c("Molecular Subtype", "", "", "AJCC Stage", "", "", "", "ER Status", "", "Survival Status", "", "Median OS Follow-up"),
  Category = c("Luminal A", "HER2-enriched", "Triple Negative", "Stage I", "Stage II", "Stage III", "Stage IV", "Positive", "Negative", "Alive", "Deceased", "—"),
  N = c(950, 37, 115, 183, 651, 251, 20, 912, 188, 1018, 87, NA),
  Pct = c("86.3%", "3.4%", "10.4%", "16.6%", "58.9%", "22.7%", "1.8%", "82.9%", "17.1%", "92.1%", "7.9%", "2.7 years"),
  stringsAsFactors = FALSE
)
doc <- add_table(doc, tbl_clinical, "表2.3 临床特征分布汇总")
doc <- body_add_break(doc)

# ---- 3. 数据挖掘算法与结果 ----
cat("\n--- Section 3: Results ---\n")
doc <- add_heading(doc, "3 数据挖掘算法与结果", 2)

# 3.1 DEG
doc <- add_heading(doc, "3.1 差异表达分析", 3)
doc <- add_text(doc, "采用DESeq2进行肿瘤组织（n=1,105）与正常组织（n=113）的差异表达分析[5]。显著性阈值取FDR<0.05（Benjamini-Hochberg校正）且|log2FC|>1。")
doc <- add_text(doc, "共鉴定6,768个显著差异表达基因，其中肿瘤相对正常上调表达4,334个、下调表达2,434个。上调基因数量约为下调的1.78倍。图3.1以火山图展示差异表达基因的全局分布：横轴为log2 Fold Change（截断±6），纵轴为-log10(p-value)（截断50），NPG色板区分上调（红）、下调（蓝）和不显著基因（灰）。上调侧红色散点密度明显高于下调侧蓝色散点，与上调基因数量占优的定量结果一致。ggrepel自动标注了差异最显著的前10个上调基因和前10个下调基因的Gene Symbol，火山图两侧的极端离群点——log2FC绝对值极大且p值极小的基因——构成了候选生物标志物的优先筛选池。两条虚线分别标记|log2FC|=1和p=0.05的显著性阈值。")

add_fig(doc, file.path(FIG_DIR, "fig1_volcano.png"),
        "图3.1 差异表达火山图（Tumor vs Normal, 6,768 DEGs, NPG红/蓝色板, ggrepel基因标注）", 5.5, 4.2)

# DEG summary table
tbl_deg <- data.frame(
  Metric = c("Total genes analyzed", "Tumor samples", "Normal samples", "Total DEGs (|log2FC|>1, FDR<0.05)", "Up-regulated", "Down-regulated"),
  Value = c("25,981", "1,105", "113", "6,768", "4,334", "2,434"),
  stringsAsFactors = FALSE
)
doc <- add_table(doc, tbl_deg, "表3.1 差异表达分析结果汇总")

# 3.2 Enrichment
doc <- add_heading(doc, "3.2 功能富集分析", 3)
doc <- add_text(doc, "采用g:Profiler在线API（gprofiler2 R包）进行功能富集分析[8]，覆盖GO、KEGG、Reactome和WikiPathway六个数据库（FDR<0.05）。分别对上调基因（4,334个）、下调基因（2,434个）和全部差异基因（6,768个）三组进行富集。")
doc <- add_text(doc, "上调基因富集到643条显著条目，集中于细胞周期和DNA复制通路（Reactome Cell Cycle, 187个基因, p=1.5×10⁻²¹）。下调基因富集到1,414条，集中于细胞外基质组织、细胞黏附和脂质代谢通路。")

tbl_enrich <- data.frame(
  Gene_Set = c("Up-regulated", "Down-regulated", "All DEGs"),
  N_Genes = c("4,334", "2,434", "6,768"),
  Significant_Terms = c("643", "1,414", "1,280"),
  Top_Pathway = c("Cell Cycle (Reactome)", "ECM Organization (Reactome)", "Cell Cycle (Reactome)"),
  Top_p_value = c("1.5×10⁻²¹", "~10⁻¹⁵", "~10⁻¹⁸"),
  stringsAsFactors = FALSE
)
doc <- add_table(doc, tbl_enrich, "表3.2 功能富集分析结果汇总")

add_fig(doc, file.path(FIG_DIR, "fig12_enrichment_up.png"),
        "图3.2 上调基因GO/KEGG富集气泡图（Top 25, FDR<0.05）", 5.8, 3.5)
add_fig(doc, file.path(FIG_DIR, "fig13_enrichment_down.png"),
        "图3.3 下调基因GO/KEGG富集气泡图（Top 25, FDR<0.05）", 5.8, 3.5)

# 3.3 Classification
doc <- add_heading(doc, "3.3 分子亚型分类", 3)
doc <- add_text(doc, "比较三种监督学习算法对BRCA分子亚型的分类性能。随机森林（ntree=500），LASSO多分类回归（5折交叉验证），XGBoost（nrounds=100）。输入特征为Top 500可变基因的log2(TPM+1)，训练集与测试集按7:3分层划分。分类标签为临床注释中的3类分子亚型（Luminal A、HER2-enriched、Triple Negative）。")
doc <- add_text(doc, "LASSO以89.4%准确率（Kappa=0.528）居首，仅使用35个区分性基因，验证了亚型间转录组差异的稀疏性。随机森林以86.7%（Kappa=0.365）次之。XGBoost准确率仅为35.8%，主要因默认超参数未针对本数据调优。")

tbl_class <- data.frame(
  Model = c("LASSO (Multinomial)", "Random Forest (500 trees)", "XGBoost"),
  Accuracy = c("0.894", "0.867", "0.358"),
  Kappa = c("0.528", "0.365", "0.002"),
  F1_Macro = c("0.612", "0.706", "0.232"),
  Features = c("35", "500", "500"),
  stringsAsFactors = FALSE
)
doc <- add_table(doc, tbl_class, "表3.3 分子亚型分类模型性能对比")

add_fig(doc, file.path(FIG_DIR, "fig6_model_comparison.png"),
        "图3.4 分类模型性能对比（Accuracy/Kappa/F1）", 5, 4)

# 3.4 Clustering
doc <- add_heading(doc, "3.4 聚类分析", 3)
doc <- add_text(doc, "综合运用PCA、t-SNE、层次聚类、K-means及共识聚类五种方法。PCA基于Top 2,000可变基因进行中心化和标准化。t-SNE采用Rtsne包（perplexity=30）。共识聚类（ConsensusClusterPlus）通过80%重采样评估K=2-8的聚类稳定性[10]。")
doc <- add_text(doc, "PCA分析显示PC1（15.3%）和PC2（8.5%）解释了最大方差，样本沿PC1主要按ER状态分离。t-SNE降维进一步展示了Luminal A与Triple Negative之间的清晰分离。K-means轮廓系数在K=2时达到最大值0.18，支持ER+/ER-二分是BRCA转录组数据中最强的自然分组结构。")

add_fig(doc, file.path(FIG_DIR, "fig2_pca_subtype.png"),
        "图3.5 PCA主成分分析（NPG四色分子亚型，PC1 15.3% × PC2 8.5%）", 5.5, 4.2)

# 3.5 WGCNA
doc <- add_heading(doc, "3.5 WGCNA共表达网络分析", 3)
doc <- add_text(doc, "使用Top 5,000可变基因在肿瘤样本中的log2(TPM+1)矩阵进行WGCNA分析[6]。软阈值选择分析显示power=8时达到scale-free R²=0.887。经动态树切割和相似模块合并后，共识别9个共表达模块。")

tbl_wgcna <- data.frame(
  Module = c("grey (unassigned)", "turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink"),
  Genes = c("2,817", "584", "553", "336", "260", "224", "119", "63", "44"),
  Top_Hub_Gene = c("ENSG00000119866", "ENSG00000115163", "ENSG00000167208", "ENSG00000123500", "ENSG00000136997", "ENSG00000075624", "ENSG00000089157", "ENSG00000131143", "ENSG00000117525"),
  MM = c("0.836", "0.889", "0.959", "0.911", "0.874", "0.903", "0.881", "0.854", "0.847"),
  stringsAsFactors = FALSE
)
doc <- add_table(doc, tbl_wgcna, "表3.4 WGCNA共表达模块概况（9个模块）")

add_fig(doc, file.path(FIG_DIR, "fig7_wgcna_sft.png"),
        "图3.6 WGCNA软阈值选择（SFT R² × 平均连接度）", 5.8, 2.4)
add_fig(doc, file.path(FIG_DIR, "fig8_wgcna_modules.png"),
        "图3.7 WGCNA模块大小分布", 4, 3)

# 3.6 Survival
doc <- add_heading(doc, "3.6 生存分析", 3)
doc <- add_text(doc, "KM曲线分析显示病理分期呈现清晰有序的逐级分离趋势（log-rank p<0.001），而分子亚型之间未检测到显著的OS差异（log-rank p=0.493）。多变量Cox模型C-index为0.771（似然比p=8.62×10⁻¹²）。在控制分期和年龄后，分子亚型与总生存无显著独立关联，验证了病理分期在BRCA预后中相对于分子亚型的预测主导地位。")

add_fig(doc, file.path(FIG_DIR, "fig3_km_stage.png"),
        "图3.8 Kaplan-Meier生存曲线（按Stage I-IV分层，log-rank p<0.0001）", 5, 4.5)

tbl_cox <- data.frame(
  Variable = c("Age", "Lymph Nodes", "Stage II (vs I)", "Stage III (vs I)", "Stage IV (vs I)", "HER2-enriched", "Triple Negative"),
  HR = c("1.02", "1.15", "1.82", "2.96", "8.74", "1.21", "1.02"),
  CI95 = c("1.00-1.04", "0.80-1.65", "0.99-3.35", "1.43-6.12", "3.61-21.14", "0.62-2.36", "0.57-1.82"),
  p_value = c("0.065", "0.453", "0.054", "0.004", "<0.001", "0.578", "0.947"),
  Significance = c("", "", "", "**", "***", "", ""),
  stringsAsFactors = FALSE
)
doc <- add_table(doc, tbl_cox, "表3.5 多变量Cox回归结果（C-index=0.771）")

add_fig(doc, file.path(FIG_DIR, "fig9_cox_forest.png"),
        "图3.9 Cox回归森林图（HR及95%CI）", 5, 4)

# 3.7 Mutation
doc <- add_heading(doc, "3.7 体细胞突变分析", 3)
doc <- add_text(doc, "采用maftools R包分析MuTect2流程输出的突变数据。突变数据集涵盖990例样本、15,413个突变基因、67,546个变异位点，平均每例约68个突变。高频突变基因为PIK3CA（369例, 37.3%）、TP53（348例, 35.2%）及TTN（268例, 27.1%）。PIK3CA与TP53之间呈现互斥趋势。")

tbl_mut <- data.frame(
  Rank = 1:10,
  Gene = c("PIK3CA", "TP53", "TTN", "CDH1", "GATA3", "MUC16", "MAP3K1", "KMT2C", "RYR2", "FLG"),
  Samples = c(369, 348, 268, 164, 148, 131, 119, 107, 98, 93),
  Frequency = c("37.3%", "35.2%", "27.1%", "16.6%", "14.9%", "13.2%", "12.0%", "10.8%", "9.9%", "9.4%"),
  Function = c("PI3K/AKT activation", "Tumor suppressor loss", "Passenger (gene length)", "Cell adhesion", "Transcription factor", "Mucin", "MAPK signaling", "Histone methyltransferase", "Ca2+ channel", "Structural protein"),
  stringsAsFactors = FALSE
)
doc <- add_table(doc, tbl_mut, "表3.6 BRCA Top 10 高频突变基因")

add_fig(doc, file.path(FIG_DIR, "fig4_oncoplot.png"),
        "图3.10 突变瀑布图（Top 12基因，990例样本）", 6, 3)

# 3.8 miRNA-mRNA
doc <- add_heading(doc, "3.8 关联规则与调控网络", 3)
doc <- add_text(doc, "采用Apriori算法进行关联规则挖掘，将Top 50差异基因按中位数离散化为High/Low，与临床特征组合构建事务数据集。miRNA-mRNA调控网络分析以|r|>0.4阈值识别出2,065对潜在调控关系对。")

add_fig(doc, file.path(FIG_DIR, "fig10_mirna_mrna_corr.png"),
        "图3.11 miRNA-mRNA调控相关性热图（Top 20, |r|>0.4）", 5, 3.8)
add_fig(doc, file.path(FIG_DIR, "fig15_association_rules.png"),
        "图3.12 关联规则散点图（Support × Confidence, lift着色）", 5, 3.5)

doc <- body_add_break(doc)

# ---- 4. Visualization ----
cat("\n--- Section 4: Visualization ---\n")
doc <- add_heading(doc, "4 可视化分析", 2)
doc <- add_text(doc, "本研究运用R语言生态系统中的多种专业绘图工具（ggplot2、pheatmap、maftools、survminer等），对各项分析结果进行了系统性的可视化表征。所有图表遵循统一规范：PNG格式、300 DPI、NPG（Nature Publishing Group）色板、英文标注。")

doc <- add_text(doc, "差异表达可视化以火山图为核心，展示6,768个DEGs的全局分布，Top 25 DEGs热图展示肿瘤与正常组织的表达模式差异。分类模型评估采用分组柱状图和ROC曲线。降维可视化通过PCA和t-SNE展示样本在高维基因空间中的自然分布，ER状态是PC1方向的最强驱动因素。WGCNA软阈值选择图和模块-性状关联热图构成网络分析的核心可视化输出。生存分析采用KM曲线和Cox森林图。突变分析以Oncoplot为核心，直观展示Top 12突变基因在990例样本中的分布。功能富集分析通过气泡图和来源数据库比较图展示基因集的功能全景。")

add_fig(doc, file.path(FIG_DIR, "fig5_deg_heatmap.png"),
        "图4.1 Top 25 DEGs表达热图（Tumor vs Normal, Z-score标准化）", 5.5, 4.2)
add_fig(doc, file.path(FIG_DIR, "fig14_enrichment_comparison.png"),
        "图4.2 富集来源数据库对比（上调 vs 下调）", 5, 3)
add_fig(doc, file.path(FIG_DIR, "fig16_enrichment_top10.png"),
        "图4.3 Top 10最显著富集通路", 5.5, 2.8)
add_fig(doc, file.path(FIG_DIR, "fig11_summary.png"),
        "图4.4 综合分析六面板汇总图", 6, 4)

doc <- body_add_break(doc)

# ---- 5. Conclusion ----
cat("\n--- Section 5: Conclusion ---\n")
doc <- add_heading(doc, "5 结论与展望", 2)
doc <- add_text(doc, "本研究以TCGA-BRCA乳腺癌队列为数据基础，整合mRNA、miRNA和体细胞突变三个组学层次，构建了覆盖差异表达分析、功能富集、分子亚型分类、聚类分析、WGCNA共表达网络、生存分析、突变分析及关联规则挖掘等九个分析模块的完整数据挖掘流程。")

doc <- add_text(doc, "主要发现包括：（1）DESeq2鉴定6,768个显著差异基因，功能富集揭示上调基因集中于细胞周期通路，下调基因集中于ECM组织通路；（2）LASSO以89.4%准确率仅用35个基因实现分子亚型三分类；（3）PCA和K-means一致表明ER状态是BRCA转录组结构的最强驱动因素；（4）WGCNA识别9个共表达模块，blue模块枢纽基因MM达0.959；（5）多变量Cox证实病理分期为BRCA预后最强预测因子（C-index=0.771, Stage IV HR=8.74）；（6）PIK3CA和TP53确认为BRCA最高频驱动突变；（7）关联规则挖掘识别了基因表达与临床特征之间的系统性关联模式。")

doc <- add_text(doc, "本研究存在以下局限性：DNA甲基化数据因数据量过大（11.7 GB）未纳入分析；生存分析事件率较低（7.9%）限制了多变量预后模型的统计功效；miRNA-mRNA调控关系基于统计相关而非因果验证。未来研究方向包括：采用降采样策略纳入DNA甲基化数据实现四组学整合；利用独立验证队列进行外部验证；对WGCNA blue模块枢纽基因进行功能实验验证；探索深度学习方法在多组学数据整合中的应用。")

doc <- body_add_break(doc)

# ---- References ----
doc <- add_heading(doc, "参考文献", 2)
refs <- c(
  "[1] Sung H, et al. Global cancer statistics 2022. CA Cancer J Clin, 2024.",
  "[2] Hanahan D, Weinberg RA. Hallmarks of cancer. Cell, 2011.",
  "[3] TCGA Network. Comprehensive molecular portraits of human breast tumours. Nature, 2012.",
  "[4] Perou CM, et al. Molecular portraits of human breast tumours. Nature, 2000.",
  "[5] Love MI, et al. DESeq2. Genome Biology, 2014.",
  "[6] Langfelder P, Horvath S. WGCNA. BMC Bioinformatics, 2008.",
  "[7] Simon N, et al. Cox LASSO. J Stat Softw, 2011.",
  "[8] Kolberg L, et al. g:Profiler 2023 update. NAR, 2023.",
  "[9] Agrawal R, Srikant R. Apriori. VLDB, 1994.",
  "[10] Wilkerson MD, Hayes DN. ConsensusClusterPlus. Bioinformatics, 2010.",
  "[11] Colaprico A, et al. TCGAbiolinks. NAR, 2016.",
  "[12] Mayakonda A, et al. Maftools. Genome Research, 2018."
)
for (ref in refs) {
  doc <- body_add_par(doc, ref, style = "Normal")
}

# ---- Save ----
cat("\nSaving document...\n")
print(doc, target = OUTPUT)
cat(sprintf("Saved: %s (%.1f KB)\n", OUTPUT, file.size(OUTPUT) / 1024))
cat("\n========== Done ==========\n")
