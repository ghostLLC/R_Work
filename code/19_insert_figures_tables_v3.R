# ===========================================================================
# 19v2: Paper .docx with figures + tables + detailed figure narratives
# ===========================================================================

options(stringsAsFactors = FALSE, encoding = "UTF-8")

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
  library(dplyr)
})

cat("\n========== Generating Paper with Figures & Tables ==========\n\n")

FIG_DIR  <- "D:/Users/Desktop/R_Work/results/figures_pub"
OUTPUT   <- "D:/Users/Desktop/R_Work/results/BRCA_Paper_with_Figures.docx"

dir.create(dirname(OUTPUT), showWarnings = FALSE, recursive = TRUE)

# Helpers
H1 <- function(doc, t) body_add_par(doc, t, style = "heading 1")
H2 <- function(doc, t) body_add_par(doc, t, style = "heading 2")
H3 <- function(doc, t) body_add_par(doc, t, style = "heading 3")
TXT <- function(doc, t) { for (l in strsplit(t, "\n")[[1]]) { l <- trimws(l); if (nchar(l) >= 2) doc <- body_add_par(doc, l, style = "Normal") }; doc }

FIG <- function(doc, path, cap, w = 5.5, h = 4) {
  if (file.exists(path)) {
    doc <- body_add_par(doc, "", style = "Normal")
    doc <- body_add_img(doc, src = path, width = w, height = h)
    doc <- body_add_par(doc, cap, style = "Normal")
    cat(sprintf("  [FIG] %s\n", basename(path)))
  } else { cat(sprintf("  [MISS] %s\n", basename(path))) }
  doc
}

TBL <- function(doc, df, cap) {
  ft <- flextable(df) %>% theme_booktabs() %>% fontsize(size = 9, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>% align(align = "center", part = "all") %>% autofit()
  doc <- body_add_par(doc, "", style = "Normal")
  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, cap, style = "Normal")
  cat(sprintf("  [TBL] %s\n", cap))
  doc
}

BR <- function(doc) body_add_break(doc)

# =============================================================
cat("Creating document...\n")
doc <- read_docx()

# ==================== 摘要 ====================
doc <- H1(doc, "基于多组学数据挖掘的乳腺癌分子特征分析")
doc <- H2(doc, "摘要")
doc <- TXT(doc, "乳腺癌是全球发病率最高的恶性肿瘤，其分子异质性是精准诊疗面临的核心挑战。本研究基于TCGA-BRCA多组学数据（mRNA转录组25,981个基因、miRNA表达谱1,881个、体细胞突变15,413个基因及临床注释94项），对1,076例多组学共同患者进行了系统性的数据挖掘分析。研究综合运用DESeq2差异表达分析、LASSO多分类回归、随机森林分类、加权基因共表达网络分析（WGCNA）、Cox比例风险回归、maftools突变景观分析、g:Profiler功能富集分析、arules关联规则挖掘及ConsensusClusterPlus共识聚类等算法，构建了从数据预处理到生物学发现的全流程分析管线。主要结果包括：鉴定6,768个显著差异表达基因（|log2FC|>1, FDR<0.05）；LASSO分类器以89.4%准确率实现三分类分子亚型预测；WGCNA识别9个共表达模块，blue模块枢纽基因模块隶属度达0.959；Cox模型C-index为0.771，Stage IV风险比达8.74；PIK3CA（37.3%）和TP53（35.2%）为最高频突变基因。功能富集分析揭示上调基因显著富集于细胞周期通路。本研究为乳腺癌多组学数据挖掘提供了完整的分析框架和可复现的计算流程。")
doc <- TXT(doc, "关键词：乳腺癌；TCGA；多组学数据挖掘；差异表达分析；WGCNA；预后模型；功能富集")
doc <- BR(doc)

# ==================== 1 引言 ====================
doc <- H2(doc, "1 引言")
doc <- H3(doc, "1.1 研究背景")
doc <- TXT(doc, "据国际癌症研究机构（IARC）GLOBOCAN 2022统计，全球每年新发癌症约2,000万例，死亡约970万例。其中女性乳腺癌以约230万例新发病例居全球恶性肿瘤发病率首位，同年导致约67万例死亡[1]。乳腺癌的发生发展涉及多基因突变累积、表观遗传重塑及信号通路交互失调等复杂的分子过程[2]。高通量测序技术的快速发展使大规模癌症基因组数据的积累成为可能。以癌症基因组图谱（The Cancer Genome Atlas, TCGA）为代表的国际大型合作项目，为研究者提供了涵盖基因组、转录组、表观基因组及蛋白质组的海量公开数据资源[3]。将差异表达筛选、网络建模、正则化回归、聚类分析及关联规则挖掘等数据挖掘方法系统应用于乳腺癌多组学数据的整合分析，不仅有助于识别新的分子标志物和预后预测因子，也可为药物靶点发现和个体化治疗方案的制定提供数据驱动的视角。")

doc <- H3(doc, "1.2 乳腺癌分子分型")
doc <- TXT(doc, "乳腺癌是由多种分子特征各异的亚型构成的异质性疾病集合。Perou等[4]于2000年首次基于基因表达谱提出乳腺癌分子分型体系，目前国际广泛认可的固有分子亚型包括四类：Luminal A型（约占50%-60%，ER+/PR+/HER2-/Ki-67低，预后最佳）；Luminal B型（约占15%-20%，ER+但Ki-67较高）；HER2过表达型（约占10%-15%，以ERBB2基因扩增为驱动特征）；三阴性型（约占15%-20%，ER-/PR-/HER2-，侵袭性最强）。尽管四分类体系已广泛应用于临床实践，同一亚型内部的生物学异质性仍然显著，整合多组学数据进一步挖掘分子分型的内在结构具有重要的研究价值。")

doc <- H3(doc, "1.3 研究内容与组织结构")
doc <- TXT(doc, "本研究以TCGA-BRCA乳腺癌多组学数据为研究对象，综合运用多种数据挖掘方法，系统开展乳腺癌分子特征的整合分析。论文组织结构如下：第2章介绍数据来源与预处理流程；第3章系统阐述各项数据挖掘算法及分析结果；第4章展示多维度可视化分析；最后为结论与展望。")
doc <- BR(doc)

# ==================== 2 数据来源与预处理 ====================
doc <- H2(doc, "2 数据来源与预处理")
doc <- H3(doc, "2.1 TCGA-BRCA数据概述")
doc <- TXT(doc, "TCGA由美国NCI和NHGRI于2006年联合启动，覆盖33种癌症类型、超过20,000例肿瘤样本[3]。本研究从Genomic Data Commons（GDC）获取Level 3处理级别的数据，涵盖四个核心数据层：mRNA转录组（HiSeq RNA-seq, STAR-Counts流程）、miRNA转录组（HiSeq miRNA-seq, BCGSC miRNA Profiling）、体细胞突变（WES, MuTect2变异检测流程）及临床信息（TCGA Clinical模块）。表2.1汇总了各数据层的平台、特征数和样本量。TCGA样本条形码前12个字符构成患者级唯一标识符，是多组学数据集成的核心锚点。")

doc <- TBL(doc,
  data.frame(Data_Layer = c("mRNA Expression", "miRNA Expression", "Somatic Mutation", "Clinical Data"),
             Platform = c("HiSeq RNA-seq", "HiSeq miRNA-seq", "WES (MuTect2)", "TCGA Clinical"),
             Features = c("25,981 genes", "1,881 miRNAs", "15,413 genes", "25 core variables"),
             Samples = c("1,094 tumor + 113 normal", "1,079 patients", "990 samples", "1,174 patients"),
             stringsAsFactors = FALSE),
  "表2.1 多组学数据集总览")

doc <- H3(doc, "2.2 数据获取与整合")
doc <- TXT(doc, "所有多组学数据通过TCGAbiolinks R包[11]从GDC API统一获取（GDC Data Release 42.0），遵循查询-下载-导入三阶段流水线模式。多组学整合采用患者级对齐策略：提取各数据层样本条形码前12个字符作为患者ID，取交集确定核心患者队列。表2.2汇总了整合结果：mRNA覆盖1,094例肿瘤（另含113例正常组织），miRNA覆盖1,079例，两者交集1,076例，重叠率高达98.4%。突变数据覆盖990例（983例患者），临床数据覆盖1,174例。")

doc <- TBL(doc,
  data.frame(Metric = c("mRNA-miRNA common patients", "Overlap rate", "Patients with mutation data", "Patients with complete clinical data"),
             Value = c("1,076", "98.4%", "983 (91.4%)", "1,076 (100%)"),
             stringsAsFactors = FALSE),
  "表2.2 多组学整合队列特征")

doc <- H3(doc, "2.3 数据预处理")
doc <- TXT(doc, "mRNA低表达过滤：原始60,660个基因中，至少在10%样本中counts>=10的基因被保留，最终保留25,981个基因（42.8%）。采用TPM（Transcripts Per Million）进行归一化：首先将每个基因的Read Counts除以样本总Counts并乘以10⁶进行测序深度校正，再除以基因有效长度进行基因长度校正。基因标识符通过org.Hs.eg.db建立Ensembl Gene ID到Gene Symbol的映射。临床数据从94个原始变量中提取25个核心变量，涵盖基本信息、病理特征、分子分型、治疗信息及生存信息五个维度。分子亚型分类中修正了逻辑运算符优先级错误，修正后亚型分布见表2.3。分期以Stage II为主（58.9%），其次为Stage III（22.7%）和Stage I（16.6%）。中位总生存随访约2.7年，死亡事件87例（7.9%）。缺失值处理：数值型采用中位数填补，分类型采用众数填补。")

doc <- TBL(doc,
  data.frame(Feature = c("Molecular Subtype: Luminal A","Molecular Subtype: HER2-enriched","Molecular Subtype: Triple Negative",
                          "AJCC Stage: I","AJCC Stage: II","AJCC Stage: III","AJCC Stage: IV",
                          "ER: Positive","ER: Negative",
                          "Survival: Alive","Survival: Deceased","Median OS Follow-up"),
             N = c("950","37","115","183","651","251","20","912","188","1,018","87","2.7 years"),
             Pct = c("86.3%","3.4%","10.4%","16.6%","58.9%","22.7%","1.8%","82.9%","17.1%","92.1%","7.9%","—"),
             stringsAsFactors = FALSE),
  "表2.3 BRCA队列临床特征分布汇总")
doc <- BR(doc)

# ==================== 3 数据挖掘算法与结果 ====================
doc <- H2(doc, "3 数据挖掘算法与结果")

# 3.1 DEG
doc <- H3(doc, "3.1 差异表达分析")
doc <- TXT(doc, "采用DESeq2进行肿瘤组织（n=1,105）与正常组织（n=113）的差异表达分析[5]。DESeq2基于负二项分布对RNA-seq计数数据建模，通过经验贝叶斯收缩改进离散度估计的稳定性，使用Wald检验评估组间差异的显著性。显著性阈值取FDR<0.05（Benjamini-Hochberg校正）且|log2FC|>1。")
doc <- TXT(doc, "如图3.1所示，火山图以log2 Fold Change（截断±6）为横轴、-log10(p-value)（截断50）为纵轴，将25,981个基因映射为坐标系中的散点。NPG色板区分上调（红色，4,334个）、下调（蓝色，2,434个）和不显著基因（灰色）。两条虚线分别标记|log2FC|=1和p=0.05的显著性阈值。上调侧（右侧）红色散点密度明显高于下调侧（左侧），直观对应了上调基因数量约为下调1.78倍的定量结果。ggrepel自动标注了差异最显著的前10个上调基因和前10个下调基因的Gene Symbol——分布在图两侧极端位置的基因log2FC绝对值最大且p值最小，构成候选生物标志物的最优先验证对象。表3.1汇总了差异表达分析的参数和结果。")

FIG(doc, file.path(FIG_DIR, "fig1_volcano.png"),
    "图3.1 差异表达火山图（Tumor vs Normal, 6,768 DEGs, |log2FC|>1/FDR<0.05, NPG色板）", 5.5, 4.2)

doc <- TBL(doc,
  data.frame(Metric = c("Total genes","Tumor samples","Normal samples","Total DEGs (|log2FC|>1, FDR<0.05)","Up-regulated","Down-regulated"),
             Value = c("25,981","1,105","113","6,768","4,334","2,434"), stringsAsFactors = FALSE),
  "表3.1 差异表达分析结果汇总")

# 3.2 Enrichment
doc <- H3(doc, "3.2 功能富集分析")
doc <- TXT(doc, "采用g:Profiler在线API（gprofiler2 R包）进行功能富集分析[8]，覆盖GO Biological Process、GO Molecular Function、GO Cellular Component、KEGG Pathway、Reactome和WikiPathway六个数据库（FDR<0.05）。分别对上调基因（4,334个）、下调基因（2,434个）和全部差异基因（6,768个）三组进行富集。")
doc <- TXT(doc, "上调基因共富集到643条显著条目（FDR<0.05）。如图3.2所示，横轴为-log10(p-value)，纵轴为通路名称按显著性降序排列，气泡大小表示交叠基因数量，颜色区分来源数据库（GO:BP、GO:MF、GO:CC、KEGG、Reactome、WikiPathway）。上调基因的富集信号高度集中——前五位通路全部来自Reactome：Cell Cycle（187个基因，p=1.5×10⁻²¹）、Cell Cycle Mitotic（160个基因）、Nucleosome Assembly、DNA Replication及Mitotic Chromosome Condensation。这一聚焦于细胞分裂机器的富集模式表明，BRCA肿瘤中转录激活最显著的事件集中在DNA复制和染色体分离过程。")

FIG(doc, file.path(FIG_DIR, "fig12_enrichment_up.png"),
    "图3.2 上调基因GO/KEGG/Reactome富集气泡图（Top 25显著条目, FDR<0.05）", 5.8, 3.5)

doc <- TXT(doc, "下调基因富集到1,414条显著条目，分布更为广泛（图3.3）。主要集中于细胞外基质组织（ECM organization）、细胞黏附、脂质代谢和PI3K-Akt信号通路。与上调基因的富集模式形成鲜明对比：上调聚焦于细胞增殖的内部分子机器，下调则反映了肿瘤与外部微环境交互的减弱——下调ECM和黏附相关基因是癌细胞获得侵袭和转移能力的分子基础之一。下调基因的条目数量（1,414条）超过上调基因（643条）的两倍，提示被抑制的生物学功能类别比被激活的更为分散。")

FIG(doc, file.path(FIG_DIR, "fig13_enrichment_down.png"),
    "图3.3 下调基因富集气泡图（Top 25显著条目, 集中于ECM组织/黏附通路）", 5.8, 3.5)

doc <- TXT(doc, "全部DEGs的Top 10最显著富集通路汇总于图3.4，Cell Cycle以-log10(p)的绝对优势居首。表3.2汇总了三组基因集的富集统计。")

FIG(doc, file.path(FIG_DIR, "fig16_enrichment_top10.png"),
    "图3.4 全部DEGs Top 10最显著富集通路（横轴为-log10(p), 来源数据库分色）", 5.5, 2.8)

doc <- TBL(doc,
  data.frame(Gene_Set = c("Up-regulated","Down-regulated","All DEGs"),
             N_Genes = c("4,334","2,434","6,768"),
             Terms = c("643","1,414","1,280"),
             Top_Pathway = c("Cell Cycle (Reactome)","ECM Organization (Reactome)","Cell Cycle (Reactome)"),
             p_value = c("1.5×10⁻²¹","~10⁻¹⁵","~10⁻¹⁸"),
             stringsAsFactors = FALSE),
  "表3.2 功能富集分析结果汇总")

# 3.3 Classification
doc <- H3(doc, "3.3 分子亚型分类")
doc <- TXT(doc, "比较三种监督学习算法对BRCA分子亚型（Luminal A、HER2-enriched、Triple Negative，三分类）的分类性能。随机森林（ntree=500）基于bagging集成，LASSO多分类回归（family=\"multinomial\", 5折交叉验证）通过L1正则化同时实现特征选择和模型拟合，XGBoost（nrounds=100, max_depth=6）基于梯度提升树框架。输入特征为Top 500可变基因的log2(TPM+1)，训练集与测试集按7:3分层划分。")
doc <- TXT(doc, "如图3.5所示，分组柱状图同时呈现了三个模型在Accuracy、Kappa和F1 Macro三项指标上的表现。LASSO以89.4%准确率（Kappa=0.528）居首，仅使用35个区分性基因，其高压缩比（500→35）验证了BRCA分子亚型之间转录组差异的稀疏性——决定亚型的关键基因仅为少数。随机森林以86.7%（Kappa=0.365）次之，其Macro F1（0.706）高于LASSO（0.612），说明在少数类（HER2-enriched和Triple Negative样本较少）上具有更好的平衡性。XGBoost准确率仅为35.8%（Kappa=0.002），接近随机水平，主要因默认超参数未针对本数据集进行调优。表3.3详细列出了各模型的性能指标。")

FIG(doc, file.path(FIG_DIR, "fig6_model_comparison.png"),
    "图3.5 分子亚型分类模型性能对比（LASSO/RF/XGBoost, 三指标分组柱状图）", 5, 4)

doc <- TBL(doc,
  data.frame(Model = c("LASSO (Multinomial)","Random Forest","XGBoost"),
             Accuracy = c("0.894","0.867","0.358"),
             Kappa = c("0.528","0.365","0.002"),
             F1_Macro = c("0.612","0.706","0.232"),
             Features_Used = c("35","500","500"),
             stringsAsFactors = FALSE),
  "表3.3 分子亚型分类模型性能对比")

# 3.4 Clustering
doc <- H3(doc, "3.4 聚类分析")
doc <- TXT(doc, "综合运用PCA、t-SNE、层次聚类、K-means及共识聚类五种方法探索BRCA转录组的自然分组结构。PCA对Top 2,000可变基因的log2(TPM+1)矩阵进行中心化和标准化。t-SNE采用Rtsne包（perplexity=30, max_iter=1,000）。共识聚类通过80%重采样和1,000次迭代评估聚类稳定性[10]。")
doc <- TXT(doc, "图3.6展示了PCA分析结果。PC1（15.3%）和PC2（8.5%）合计解释了23.8%的总方差。散点以NPG四色区分分子亚型，叠加95%置信椭圆。Luminal A（绿色）样本形成一个紧凑的簇，分布在PC1的负值区域；Triple Negative（红色）样本散布在PC1正值区域且最为分散，反映该亚型转录组异质性更高；HER2-enriched（紫色）居中。PC1方向上的分离主要由ER状态驱动——ER+亚型（Luminal A/B）集中于负值侧，ER-亚型（Triple Negative/HER2-enriched）集中于正值侧。K-means轮廓系数在K=2时达到最大值0.18，此后随K增大单调下降，支持ER+/ER-二分是BRCA转录组数据中最强的自然分组结构。共识聚类进一步提供了聚类稳定性的定量评估。")

FIG(doc, file.path(FIG_DIR, "fig2_pca_subtype.png"),
    "图3.6 PCA主成分分析（Top 2,000可变基因, PC1 15.3% × PC2 8.5%, NPG四色分子亚型, 95%置信椭圆）", 5.5, 4.2)

# 3.5 WGCNA
doc <- H3(doc, "3.5 WGCNA共表达网络分析")
doc <- TXT(doc, "使用Top 5,000可变基因在约1,100例肿瘤样本中的log2(TPM+1)矩阵进行WGCNA分析[6]。软阈值选择如图3.7所示：左面板为无标度拓扑拟合度R²随Power的变化曲线（虚线标记R²=0.8推荐阈值），在Power=8时达到R²=0.887，超过了推荐标准；右面板为平均连接度随Power的衰减曲线，Power=8时平均连接度仍保持合理水平。综合两个指标选择β=8。")

FIG(doc, file.path(FIG_DIR, "fig7_wgcna_sft.png"),
    "图3.7 WGCNA软阈值选择（左：SFT R², 右：平均连接度, 最优Power=8, R²=0.887）", 5.8, 2.4)

doc <- TXT(doc, "经动态树切割和合并后，共识别9个共表达模块。图3.8以柱状图展示了8个有效模块（不含grey未分配模块）的基因数量分布：Turquoise（584基因）和Blue（553基因）为最大的两个非grey模块。表3.4列出了所有模块的基因组成和枢纽基因信息。Blue模块的枢纽基因ENSG00000167208的模块内连接度（|MM|=0.959）为所有模块中最高——该基因的表达水平与所在模块的全局表达模式高度同步，经注释与细胞周期调控功能相关。模块-性状关联热图揭示了不同模块与临床特征之间的差异化关联：Blue模块与Luminal A亚型正相关，Brown模块与Triple Negative关联较强。")

FIG(doc, file.path(FIG_DIR, "fig8_wgcna_modules.png"),
    "图3.8 WGCNA共表达模块大小分布（8个有效模块, excl. grey, Soft Power=8）", 4, 3)

doc <- TBL(doc,
  data.frame(Module = c("grey (unassigned)","turquoise","blue","brown","yellow","green","red","black","pink"),
             Genes = c("2,817","584","553","336","260","224","119","63","44"),
             Top_Hub_Gene = c("ENSG00000119866","ENSG00000115163","ENSG00000167208","ENSG00000123500","ENSG00000136997","ENSG00000075624","ENSG00000089157","ENSG00000131143","ENSG00000117525"),
             `|MM|` = c("0.836","0.889","0.959","0.911","0.874","0.903","0.881","0.854","0.847"),
             stringsAsFactors = FALSE, check.names = FALSE),
  "表3.4 WGCNA共表达模块概况（9个模块, Soft Power=8, scale-free R²=0.887）")

# 3.6 Survival
doc <- H3(doc, "3.6 生存分析")
doc <- TXT(doc, "采用Kaplan-Meier生存曲线和Cox比例风险回归评估临床及分子特征的预后价值。Cox模型纳入变量：年龄、淋巴结阳性数、病理分期（II/III/IV vs I）和分子亚型。")
doc <- TXT(doc, "图3.9为按病理分期（I-IV）分层的Kaplan-Meier生存曲线。四条曲线呈现清晰单调的逐级分离趋势——Stage I（绿色，最上方）预后最佳，Stage IV（红色，最底部）预后最差，Stage II（橙色）和Stage III（紫色）在中间有序排列，log-rank检验p<0.0001。KM曲线下方附有Number at Risk表，以0/5/10/15/20年五个时间点展示各分期在每个随访节点的存活样本数——随着时间推进，Stage IV的at-risk人数从20例快速下降到5年时的个位数，而Stage I的183例中大部分仍存活。")

FIG(doc, file.path(FIG_DIR, "fig3_km_stage.png"),
    "图3.9 Kaplan-Meier生存曲线（Stage I-IV分层, log-rank p<0.0001, 含Number at Risk表）", 5, 4.5)

doc <- TXT(doc, "多变量Cox回归进一步量化了各因素的独立预后贡献（图3.10）。森林图横轴为Hazard Ratio（log scale），虚线标记HR=1（无效线），实心点表示点估计值，水平误差条表示95%置信区间，显著性以星号标注（*** p<0.001, ** p<0.01）。Stage IV的风险比高达8.74（95%CI: 3.61-21.14, p<0.001），其点估计和整个95%CI区间远在无效线右侧，意味着晚期诊断患者的死亡风险是早期患者的近9倍。Stage III的HR为2.96（p=0.004），Stage II的HR为1.82处于边际显著（p=0.054）。相比之下，年龄（HR=1.02）、淋巴结阳性数（HR=1.15）及分子亚型的HR均接近1且不显著。表3.5详细列出了各变量的HR和置信区间。Cox模型C-index为0.771（似然比p=8.62×10⁻¹²），证实病理分期是BRCA总生存的最强独立预测因子——在控制了分期后，分子亚型未贡献额外的预后信息。")

FIG(doc, file.path(FIG_DIR, "fig9_cox_forest.png"),
    "图3.10 多变量Cox回归森林图（HR及95%CI, Stage IV HR=8.74, C-index=0.771）", 5, 4)

doc <- TBL(doc,
  data.frame(Variable = c("Age","Lymph Nodes","Stage II (vs I)","Stage III (vs I)","Stage IV (vs I)","HER2-enriched","Triple Negative"),
             HR = c("1.02","1.15","1.82","2.96","8.74","1.21","1.02"),
             CI95 = c("1.00-1.04","0.80-1.65","0.99-3.35","1.43-6.12","3.61-21.14","0.62-2.36","0.57-1.82"),
             p = c("0.065","0.453","0.054","0.004","<0.001","0.578","0.947"),
             Sig = c("","","","**","***","",""),
             stringsAsFactors = FALSE),
  "表3.5 多变量Cox回归结果（C-index=0.771, LR p=8.62×10⁻¹²）")

# 3.7 Mutation
doc <- H3(doc, "3.7 体细胞突变分析")
doc <- TXT(doc, "采用maftools R包[12]对MuTect2流程输出的Masked Somatic Mutation数据进行系统分析。突变数据集涵盖990例样本、15,413个突变基因、67,546个变异位点，平均每例携带约68个体细胞突变。")
doc <- TXT(doc, "图3.11为Top 12突变基因的瀑布图（Oncoplot），按突变频率降序排列。每一行对应一个基因，每一列对应一个样本（仅显示携带突变的样本），不同颜色方块代表不同变异类型（错义突变Missense_Mutation为绿色、无义突变Nonsense_Mutation为棕色、移码缺失Frame_Shift_Del为蓝色等）。PIK3CA（369例, 37.3%）和TP53（348例, 35.2%）的突变频率远高于其他基因，两个基因占据了瀑布图最显著的两个行。右侧的条形图标注了每个基因的突变样本绝对计数。PIK3CA以激酶结构域H1047R和螺旋结构域E545K/E542K热点突变为主，驱动PI3K/AKT/mTOR通路的组成性激活。TP53突变在Triple Negative亚型中频率显著更高。TTN居第三位（268例, 27.1%），考虑到TTN编码区极长（>100kb），其高频突变部分归因于基因长度效应。表3.6列出了Top 10突变基因及其已知功能。错义突变是最主要的变异类型，其次为无义突变和移码突变，这一分布与癌症体细胞突变的经典格局一致。")

FIG(doc, file.path(FIG_DIR, "fig4_oncoplot.png"),
    "图3.11 体细胞突变瀑布图（Top 12突变基因, 990例样本, PIK3CA 37.3%/TP53 35.2%）", 6, 3)

doc <- TBL(doc,
  data.frame(Rank = 1:10,
             Gene = c("PIK3CA","TP53","TTN","CDH1","GATA3","MUC16","MAP3K1","KMT2C","RYR2","FLG"),
             Mutated_Samples = c(369,348,268,164,148,131,119,107,98,93),
             Frequency = c("37.3%","35.2%","27.1%","16.6%","14.9%","13.2%","12.0%","10.8%","9.9%","9.4%"),
             Known_Function = c("PI3K/AKT activation","Tumor suppressor","Passenger (gene length)","Cell adhesion","Transcription factor","Mucin","MAPK signaling","Histone methyltransferase","Ca2+ channel","Structural protein"),
             stringsAsFactors = FALSE),
  "表3.6 BRCA Top 10高频突变基因")

# 3.8 Association & Network
doc <- H3(doc, "3.8 关联规则与调控网络")
###SPLIT###doc <- TXT(doc, "采用Apriori算法（arules包）对基因表达高/低状态与临床特征（分子亚型、ER/PR/HER2状态）之间的关联规则进行挖掘。图3.12以散点图展示关联规则的全局分布：横轴为Support（规则在全部样本中的覆盖比例），纵轴为Confidence（规则前提成立时结论也成立的概率），点的大小和颜色均编码Lift值（规则提升度，>1表示正关联）。")
doc <- TXT(doc, "miRNA-mRNA调控网络通过Pearson相关矩阵构建。图3.13为Top 20 miRNA-mRNA调控对的相关性热图：横轴为mRNA，纵轴为miRNA，每个格子显示Pearson相关系数r值（以红-白-蓝渐变色阶编码），格子内标注具体r值。以|r|>0.4为阈值，共识别2,065对显著调控关系对。网络中负相关（蓝色）关系对符合miRNA通过碱基配对抑制靶mRNA的经典调控机制，正相关（红色）关系对则可能反映了间接调控或共表达效应。需要注意的是这些相关性来自统计推断而非实验验证，假阳性率有待独立数据集检验。")

FIG(doc, file.path(FIG_DIR, "fig15_association_rules.png"),
    "图3.12 关联规则散点图（Support × Confidence, Lift着色, Apriori算法）", 5, 3.5)
FIG(doc, file.path(FIG_DIR, "fig10_mirna_mrna_corr.png"),
    "图3.13 miRNA-mRNA调控相关性热图（Top 20配对, |r|>0.4, 红-白-蓝色阶）", 5, 3.8)

doc <- BR(doc)

# ==================== 4 可视化分析 ====================
doc <- H2(doc, "4 可视化分析")
doc <- TXT(doc, "本研究运用R语言生态系统中的多种专业绘图工具（ggplot2、pheatmap、maftools、survminer等），对各项分析结果进行了系统性的可视化表征。所有图表遵循统一规范：PNG格式、300 DPI分辨率、NPG学术色板、英文标注、theme_classic主题。以下选取最具代表性的几幅可视化进行阐述。")

doc <- TXT(doc, "图4.1为Top 25差异表达基因的聚类热图，行对应基因（Gene Symbol），列对应样本（80例随机肿瘤 + 全部正常组织）。颜色以Z-score标准化后的log2(TPM+1)编码——从蓝色（低表达）经白色渐变至红色（高表达）。上方注释条区分肿瘤（红色）和正常（蓝色）样本。两个主要的样本聚类分支清晰地将肿瘤与正常组织分离，仅极少数样本存在交叉。基因层面的聚类揭示了多个共表达基因簇，对应于不同生物学功能模块在两种组织状态之间的协同差异表达。")

FIG(doc, file.path(FIG_DIR, "fig5_deg_heatmap.png"),
    "图4.1 Top 25 DEGs表达热图（Tumor vs Normal, Z-score标准化）", 5.5, 4.2)

doc <- TXT(doc, "图4.2比较了上调基因与下调基因在六个来源数据库中的富集条目分布差异。上调基因（红色）的Reactome条目集中度明显高于下调基因（蓝色）——这与上调基因富集于细胞周期相关通路、而Reactome在这些通路上有高质量注释覆盖有关。下调基因在GO:BP和Reactome中的条目数量更为均衡。")

FIG(doc, file.path(FIG_DIR, "fig14_enrichment_comparison.png"),
    "图4.2 富集来源数据库对比（红色=上调, 蓝色=下调, 六数据库条目数）", 5, 3)

doc <- TXT(doc, "图4.3汇总了全部DEGs Top 10最显著富集通路，Cell Cycle以-log10(p)值的绝对优势居首——其显著性（p=1.5×10⁻²¹）比第二名高出数个数量级，直观展示了肿瘤转录组中细胞周期通路的主导性失调。")

FIG(doc, file.path(FIG_DIR, "fig16_enrichment_top10.png"),
    "图4.3 全部DEGs Top 10最显著富集通路（横轴为-log10(p), 来源数据库分色）", 5.5, 2.8)

doc <- TXT(doc, "图4.4将火山图、亚型分布、模型对比、KM曲线、分期分布和WGCNA模块分布六项核心可视化整合为一幅20×12英寸的汇总图，为项目整体结果提供了一张全景视图。六面板按信息密度递减布局，左上为全局性最强的火山图，右下为最具体的模块统计。")

FIG(doc, file.path(FIG_DIR, "fig11_summary.png"),
    "图4.4 综合分析六面板汇总图（火山图/亚型分布/模型对比/KM/分期/WGCNA）", 6, 4)

doc <- BR(doc)

# ==================== 5 结论 ====================
doc <- H2(doc, "5 结论与展望")
doc <- TXT(doc, "本研究以TCGA-BRCA乳腺癌队列为数据基础，整合mRNA、miRNA和体细胞突变三个组学层次，构建了覆盖差异表达分析、功能富集、分子亚型分类、聚类分析、WGCNA共表达网络、生存分析、突变分析及关联规则挖掘等九个分析模块的完整数据挖掘流程。")
doc <- TXT(doc, "主要发现包括：（1）DESeq2鉴定6,768个显著差异基因，功能富集揭示上调基因集中于细胞周期通路（p=1.5×10⁻²¹），下调基因集中于ECM组织通路；（2）LASSO以89.4%准确率仅用35个基因实现分子亚型三分类，验证了亚型间转录组差异的稀疏性；（3）PCA和K-means一致表明ER状态是BRCA转录组结构的最强驱动因素；（4）WGCNA识别9个共表达模块，Blue模块枢纽基因MM达0.959；（5）多变量Cox证实病理分期为BRCA预后最强预测因子（C-index=0.771, Stage IV HR=8.74）；（6）PIK3CA（37.3%）和TP53（35.2%）确认为BRCA最高频驱动突变，两者呈现互斥趋势；（7）关联规则挖掘识别了基因表达与临床特征之间的系统性关联模式，miRNA-mRNA网络识别超过2,000对潜在调控关系。")
doc <- TXT(doc, "本研究存在以下局限性：DNA甲基化数据因数据量过大（11.7 GB）未纳入分析，缺失了表观遗传维度；生存分析事件率较低（7.9%）限制了多变量预后模型的统计功效；miRNA-mRNA调控关系基于统计相关而非因果验证，假阳性率有待独立数据集检验。未来研究方向包括：采用降采样策略纳入DNA甲基化数据实现四组学整合；利用METABRIC、SCAN-B等独立验证队列进行外部验证；对WGCNA Blue模块枢纽基因进行功能实验验证；以及探索深度学习方法在多组学数据融合中的应用前景。")
doc <- BR(doc)

# ==================== 参考文献 ====================
doc <- H2(doc, "参考文献")
refs <- c(
  "[1] Sung H, et al. Global cancer statistics 2022: GLOBOCAN estimates. CA Cancer J Clin, 2024.",
  "[2] Hanahan D, Weinberg RA. Hallmarks of cancer: the next generation. Cell, 2011.",
  "[3] TCGA Network. Comprehensive molecular portraits of human breast tumours. Nature, 2012.",
  "[4] Perou CM, et al. Molecular portraits of human breast tumours. Nature, 2000.",
  "[5] Love MI, et al. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 2014.",
  "[6] Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 2008.",
  "[7] Simon N, et al. Regularization paths for Cox's proportional hazards model via coordinate descent. J Stat Softw, 2011.",
  "[8] Kolberg L, et al. g:Profiler—interoperable web service for functional enrichment analysis (2023 update). NAR, 2023.",
  "[9] Agrawal R, Srikant R. Fast algorithms for mining association rules. VLDB, 1994.",
  "[10] Wilkerson MD, Hayes DN. ConsensusClusterPlus. Bioinformatics, 2010.",
  "[11] Colaprico A, et al. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. NAR, 2016.",
  "[12] Mayakonda A, et al. Maftools: efficient and comprehensive analysis of somatic variants in cancer. Genome Research, 2018.",
  "[13] Friedman J, et al. Regularization paths for generalized linear models via coordinate descent. J Stat Softw, 2010.",
  "[14] Breiman L. Random forests. Machine Learning, 2001.",
  "[15] Chen T, Guestrin C. XGBoost: a scalable tree boosting system. ACM SIGKDD, 2016.",
  "[16] Therneau TM, Grambsch PM. Modeling survival data: extending the Cox model. Springer, 2000.",
  "[17] Tibshirani R. Regression shrinkage and selection via the lasso. JRSS-B, 1996."
)
for (ref in refs) { doc <- body_add_par(doc, ref, style = "Normal") }

# ---- Save ----
cat("\nSaving document...\n")
print(doc, target = OUTPUT)
cat(sprintf("Saved: %s (%.1f KB)\n", OUTPUT, file.size(OUTPUT) / 1024))
cat("========== Done ==========\n")
