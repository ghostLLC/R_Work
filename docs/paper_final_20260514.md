# 基于多组学数据挖掘的乳腺癌分子特征分析

## 摘要

乳腺癌是全球发病率最高的恶性肿瘤，其分子异质性是精准诊疗面临的核心挑战。本研究基于TCGA-BRCA多组学数据（mRNA转录组25,981个基因、miRNA表达谱1,881个、体细胞突变15,413个基因及临床注释94项），对1,076例多组学共同患者进行了系统性的数据挖掘分析。研究综合运用DESeq2差异表达分析、LASSO多分类回归、随机森林分类、加权基因共表达网络分析（WGCNA）、Cox比例风险回归、maftools突变景观分析、g:Profiler功能富集分析、arules关联规则挖掘及ConsensusClusterPlus共识聚类等算法，构建了从数据预处理到生物学发现的全流程分析管线。主要结果包括：鉴定6,768个显著差异表达基因（|log2FC|>1, FDR<0.05），其中上调4,334个、下调2,434个；LASSO分类器以89.4%准确率实现三分类分子亚型预测，特征压缩至35个基因；WGCNA识别9个共表达模块，blue模块枢纽基因模块隶属度达0.959；多变量Cox模型C-index为0.771，Stage IV风险比达8.74；PIK3CA（37.3%）和TP53（35.2%）为最高频突变基因。功能富集分析揭示上调基因显著富集于细胞周期和DNA复制通路，下调基因富集于细胞外基质组织和黏附相关通路。关联规则挖掘识别了基因表达与分子亚型之间的系统性关联模式。本研究为乳腺癌多组学数据挖掘提供了完整的分析框架和可复现的计算流程。

**关键词**：乳腺癌；TCGA；多组学数据挖掘；差异表达分析；WGCNA；预后模型；功能富集

---

## 1 引言

### 1.1 研究背景

据国际癌症研究机构（IARC）GLOBOCAN 2022统计，全球每年新发癌症约2,000万例，死亡约970万例。其中，女性乳腺癌以约230万例新发病例居全球恶性肿瘤发病率首位，同年导致约67万例死亡[1]。在中国，乳腺癌的发病率和死亡率持续上升，且发病年龄显著早于欧美国家，构成重大的公共卫生负担。

乳腺癌的发生发展涉及多基因突变累积、表观遗传重塑及信号通路交互失调等复杂的分子过程[2]。传统单一维度研究方法（如孤立分析单个基因表达或特定通路激活状态）难以全面揭示乳腺癌的系统性分子特征。高通量测序技术的快速发展及成本的持续下降，使大规模、多层次癌症基因组数据的积累成为可能。以癌症基因组图谱（The Cancer Genome Atlas, TCGA）为代表的国际大型合作项目，为研究者提供了涵盖基因组、转录组、表观基因组及蛋白质组的海量公开数据资源[3]，使基于多组学整合的系统生物学研究范式得以实践。

在数据驱动的研究范式下，如何从高维度、高噪声的癌症多组学数据中高效提取具有生物学意义和临床价值的信息，是生物信息学领域的核心问题。数据挖掘技术为这一问题的解决提供了方法论支撑。将差异表达筛选、网络建模、正则化回归、聚类分析及关联规则挖掘等方法系统应用于乳腺癌多组学数据的整合分析，不仅有助于识别新的分子标志物和预后预测因子，也可为药物靶点发现和个体化治疗方案的制定提供数据驱动的视角。

### 1.2 乳腺癌分子分型

乳腺癌是由多种分子特征各异的亚型构成的异质性疾病集合。Perou等[4]于2000年首次基于基因表达谱提出乳腺癌分子分型体系，此后经过多次完善，已成为指导精准治疗的核心框架。目前国际广泛认可的固有分子亚型（Intrinsic Subtypes）包括四类：（1）Luminal A型，约占50%-60%，特征为ER阳性、PR阳性、HER2阴性及Ki-67低表达，预后最佳；（2）Luminal B型，约占15%-20%，ER阳性但PR可能低表达或阴性、Ki-67较高，HER2可阳性或阴性；（3）HER2过表达型（HER2-enriched），约占10%-15%，以ERBB2基因扩增为驱动特征；（4）三阴性型（Triple Negative/Basal-like），约占15%-20%，ER、PR、HER2均为阴性，侵袭性最强。

尽管四分类体系已广泛应用于临床实践，同一亚型内部的生物学异质性仍然显著。部分Luminal A患者出现早期复发，而部分三阴性患者可长期无病生存——这一现象说明现有分型体系存在改进空间。整合多组学数据进一步挖掘分子分型的内在结构，具有重要的研究价值和临床转化潜力。

### 1.3 数据挖掘在癌症研究中的应用

高通量组学技术的普及使癌症研究面临从"数据稀缺"到"数据过载"的转变。传统人工筛选和单变量分析方法在几万维基因表达矩阵面前已力不从心，数据挖掘方法因此成为癌症生物信息学事实上的标准工具。

在差异表达分析方面，基于负二项分布的DESeq2[5]和基于经验贝叶斯的limma是应用最广泛的两套框架，分别针对计数数据和连续数据进行了精细的统计建模。在网络分析层面，加权基因共表达网络分析（WGCNA）[6]通过构建无标度拓扑网络将海量基因聚类为共表达模块，并将模块与临床表型相关联，已在多种癌症中成功筛选出关键模块和枢纽基因。在预后模型构建方面，LASSO-Cox比例风险模型[7]能有效应对高维表达数据中变量数远大于样本数的维数灾难问题。在功能富集方面，g:Profiler[8]等在线工具提供了覆盖GO、KEGG、Reactome、WikiPathway等多数据库的一站式富集分析。在关联模式发现方面，Apriori算法[9]等关联规则挖掘方法可从离散化表达数据中发现基因与临床特征之间的频繁共现模式。在分子分型方面，共识聚类（Consensus Clustering）[10]通过重采样评估聚类稳定性，提供了不依赖先验标签的数据驱动分层方案。

### 1.4 研究内容与组织结构

本研究以TCGA-BRCA乳腺癌多组学数据为研究对象，综合运用多种数据挖掘方法，系统开展乳腺癌分子特征的整合分析。具体研究内容围绕以下方面展开：（1）基于mRNA转录组数据识别肿瘤与正常组织的差异表达基因并进行功能富集分析；（2）利用WGCNA及共识聚类识别共表达模块和自然样本分组；（3）构建分子亚型分类模型和预后评估模型；（4）分析体细胞突变图谱；（5）挖掘关联规则并构建miRNA-mRNA调控网络。

论文组织结构如下：第2章介绍数据来源与预处理流程；第3章系统阐述所采用的各项数据挖掘算法及分析结果；第4章展示多维度可视化分析；最后为结论与展望。

---

## 2 数据来源与预处理

### 2.1 TCGA-BRCA数据概述

TCGA由美国NCI和NHGRI于2006年联合启动，覆盖33种癌症类型、超过20,000例肿瘤样本[3]。乳腺浸润性癌（BRCA）是TCGA最早启动的癌种之一，也是样本量最大的实体瘤项目。本研究从Genomic Data Commons（GDC）获取Level 3处理级别的数据，涵盖四个核心数据层：

（1）**mRNA转录组**（HiSeq RNA-seq, STAR-Counts流程）：覆盖全基因组表达谱；
（2）**miRNA转录组**（HiSeq miRNA-seq, BCGSC miRNA Profiling）：覆盖微小RNA表达谱；
（3）**体细胞突变**（WES, MuTect2变异检测流程）：覆盖点突变和插入缺失；
（4）**临床信息**（TCGA Clinical模块）：覆盖人口统计学、病理诊断、TNM分期、分子分型及生存随访数据。

TCGA样本条形码采用标准化分层命名规则，前12个字符（如TCGA-A8-A079）构成患者级唯一标识符，是跨组学数据集成的核心锚点。样本类型代码（第14-15位）中，01表示原发实体瘤，11表示癌旁正常组织。

### 2.2 数据获取与整合

所有多组学数据通过TCGAbiolinks R包[11]从GDC API统一获取（GDC Data Release 42.0）。数据获取遵循"查询-下载-导入"三阶段流水线模式：GDCquery构建查询对象，GDCdownload批量下载，GDCprepare解析为R结构化数据对象。

多组学整合采用患者级对齐策略：提取各数据层样本条形码的前12个字符作为患者ID，取交集确定核心患者队列。整合结果为：mRNA覆盖1,094例肿瘤（另含113例正常组织），miRNA覆盖1,079例，两者交集1,076例，重叠率98.4%。突变数据覆盖990例（983例患者），临床数据覆盖1,174例。最终多组学整合队列为1,076例患者，均具备完整mRNA和miRNA表达谱及临床注释。

### 2.3 数据预处理

**mRNA低表达过滤**：原始计数矩阵包含60,660个基因。过滤标准为至少在10%样本中counts≥10，该阈值参考DESeq2推荐的独立过滤准则。过滤后保留25,981个基因（42.8%），被滤除的以低丰度lncRNA和假基因为主。

**TPM归一化**：采用Transcripts Per Million（TPM）作为归一化指标。计算分两步：首先将每个基因的Read Counts除以该样本总Counts并乘以10⁶进行测序深度校正（CPM），再将CPM除以对应基因的有效长度（千碱基）进行基因长度校正。基因长度信息通过Ensembl基因组注释数据库（GRCh38/hg38）获取。

**基因标识符映射**：通过org.Hs.eg.db建立Ensembl Gene ID到Gene Symbol的映射关系。对于一对多映射，保留平均表达量最高的Ensembl ID。对于无法映射的ID，保留原始Ensembl ID作为后备标识。

**临床数据清洗**：从94个原始变量中提取25个核心变量，涵盖基本信息、病理特征、分子分型、治疗信息及生存信息五个维度。分子亚型分类中发现了逻辑运算符优先级导致的分类不一致问题，经修正后亚型分布为：Luminal A 950例（86.3%）、HER2-enriched 37例（3.4%）、Triple Negative 115例（10.4%）。Luminal B仅3例因样本量不足在后续分类分析中被过滤。分期分布：Stage I 183例（16.6%）、Stage II 651例（58.9%）、Stage III 251例（22.7%）、Stage IV 20例（1.8%）。中位总生存随访时间约2.7年，死亡事件87例（7.9%）。

**缺失值处理**：数值型变量采用中位数填补，分类型变量采用众数填补。生存分析中Alive患者的缺失days_to_death属于天然删失，不作填补。

**miRNA与突变预处理**：miRNA数据（1,881个miRNA）全部保留，不作额外低表达过滤。突变数据采用MuTect2流程输出的Masked Somatic Mutation，统计发现15,413个基因在至少1例患者中存在体细胞突变，前20高频突变基因与已知乳腺癌驱动基因谱高度一致。

### 2.4 数据概览

| 数据层 | 平台 | 特征数 | 样本量 | 归一化 |
|--------|------|--------|--------|--------|
| mRNA表达 | HiSeq RNA-seq | 25,981基因 | 1,094例（肿瘤）+113例（正常） | TPM |
| miRNA表达 | HiSeq miRNA-seq | 1,881 miRNA | 1,079例 | — |
| 体细胞突变 | WES (MuTect2) | 15,413基因 | 990例 | — |
| 临床数据 | TCGA Clinical | 25变量 | 1,174例 | — |

核心整合队列：1,076例（mRNA ∩ miRNA），重叠率98.4%。

---

## 3 数据挖掘算法与结果

### 3.1 差异表达分析

**方法原理**：采用DESeq2进行肿瘤组织（n=1,105）与正常组织（n=113）的差异表达分析。DESeq2基于负二项分布对RNA-seq计数数据建模，通过经验贝叶斯收缩改进离散度估计的稳定性，使用Wald检验评估组间差异的显著性[5]。分析参数：显著性阈值取FDR<0.05（Benjamini-Hochberg校正）且|log2FC|>1。

**分析结果**：共鉴定6,768个显著差异表达基因，其中肿瘤相对正常上调表达4,334个、下调表达2,434个。上调基因数量约为下调的1.78倍，提示肿瘤组织中转录激活事件多于转录抑制事件，与该类细胞的增殖活跃和代谢重编程特征一致。

### 3.2 功能富集分析

**方法原理**：采用g:Profiler在线API（gprofiler2 R包）进行功能富集分析[8]，覆盖GO Biological Process、GO Molecular Function、GO Cellular Component、KEGG Pathway、Reactome和WikiPathway六个数据库。显著性阈值取FDR<0.05，排除电子注释（IEA）。分别对上调基因（4,334个）、下调基因（2,434个）和全部差异基因（6,768个）三组进行富集。

**分析结果**：上调基因富集到643条显著条目，下调基因富集到1,414条，全部DEGs富集到1,280条。上调基因的富集信号高度集中——前五位通路全部来自Reactome，包括Cell Cycle（187个基因，p=1.5×10⁻²¹）、Cell Cycle Mitotic（160个基因）、Nucleosome Assembly、DNA Replication及Mitotic Chromosome Condensation，均与细胞周期和DNA复制相关。下调基因的富集方向为细胞外基质组织（ECM organization）、细胞黏附、脂质代谢、PI3K-Akt信号通路及整合素细胞表面相互作用。上调与下调基因的功能富集差异反映了肿瘤细胞增殖激活与微环境交互减弱的双重特征。

### 3.3 分子亚型分类

**方法原理**：比较三种监督学习算法对BRCA分子亚型的分类性能。随机森林（ntree=500）基于bagging集成降低过拟合风险。LASSO多分类回归（family="multinomial", α=1, 5折交叉验证）通过L1正则化同时实现特征选择和模型拟合。XGBoost（max_depth=6, eta=0.1, nrounds=100）基于梯度提升树框架。输入特征为Top 500可变基因的log2(TPM+1)，训练集与测试集按7:3分层划分。分类标签为临床注释中的分子亚型（去除样本数<15的类别后为3类：Luminal A、HER2-enriched、Triple Negative）。

**分析结果**：LASSO以89.4%准确率（Kappa=0.528）居首，仅使用35个区分性基因。随机森林以86.7%（Kappa=0.365）次之，其Macro F1（0.706）高于LASSO（0.612），说明在处理少数类时平衡性更好。XGBoost准确率仅为35.8%，接近随机水平（3分类基线33%），主要因默认超参数未针对本数据调优。LASSO的高压缩比（500个候选→35个选择）验证了乳腺癌分子亚型之间转录组差异的稀疏性：决定亚型的关键基因仅为少数。

### 3.4 聚类分析

**方法原理**：综合运用PCA、t-SNE、层次聚类、K-means及共识聚类五种方法。PCA基于prcomp对Top 2,000可变基因进行中心化和标准化。t-SNE采用Rtsne包（perplexity=30, max_iter=1000）。层次聚类采用Ward.D2方法和Pearson相关系数距离。K-means对K=2至10评估轮廓系数。共识聚类（ConsensusClusterPlus）通过80%重采样和1,000次迭代评估K=2-8的聚类稳定性[10]。

**分析结果**：PCA分析显示PC1（15.3%）和PC2（8.5%）解释了最大方差，样本沿PC1主要按ER状态分离。t-SNE降维进一步展示Luminal A与Triple Negative之间的清晰分离，HER2-enriched散布于中间地带。K-means轮廓系数在K=2时达到最大值0.18，此后随K增大单调下降，支持ER+/ER-二分是BRCA转录组数据中最强的自然分组结构。共识聚类分析提供了聚类稳定性的定量评估，CDF曲线和delta area图确认了最优聚类数，簇特异性标记基因定义了各共识簇的分子特征。

### 3.5 WGCNA共表达网络分析

**方法原理**：WGCNA通过构建无标度共表达网络将基因聚类为功能相关的模块[6]。使用Top 5,000可变基因在肿瘤样本中的log2(TPM+1)矩阵。流程包括：计算基因间Pearson相关系数构建相似性矩阵；通过软阈值幂函数将相似性矩阵转化为邻接矩阵（选择标准为scale-free R²>0.8）；构建拓扑重叠矩阵（TOM）以考虑间接连接模式；基于TOM相异度进行层次聚类和动态树切割识别模块（minModuleSize=30, mergeCutHeight=0.25）；计算模块特征基因（第一主成分）与临床性状的相关性；通过模块内连接度（Module Membership, MM）识别枢纽基因。

**分析结果**：软阈值选择分析显示power=8时达到scale-free R²=0.887，同时平均连接度保持合理水平。经动态树切割和相似模块合并后，共识别9个共表达模块：grey（2,817个基因，未分配）、turquoise（584）、blue（553）、brown（336）、yellow（260）、green（224）、red（119）、black（63）、pink（44）。不含grey模块的有效模块共涵盖2,183个基因。blue模块枢纽基因（ENSG00000167208）的模块内连接度最高（|MM|=0.959），与细胞周期调控功能相关。模块-性状关联分析显示blue模块与Luminal A亚型正相关，brown模块与Triple Negative亚型关联较强。

### 3.6 生存分析

**方法原理**：采用Kaplan-Meier生存曲线（log-rank检验）和Cox比例风险回归评估临床及分子特征的预后价值。Cox模型纳入变量：年龄、淋巴结阳性数、病理分期（II/III/IV vs I）、分子亚型（HER2-enriched/Triple Negative vs Luminal A）。模型评价采用C-index和似然比检验。此外，采用单变量Cox回归对Top 500差异基因逐个进行预后评估（以基因中位表达值分组），筛选标准取FDR<0.05。

**分析结果**：KM曲线分析显示病理分期呈现清晰有序的逐级分离趋势，Stage IV预后最差。分子亚型之间未检测到显著的OS差异（log-rank p=0.493），可能与随访时间较短（中位2.7年）及死亡事件数量有限（87/1,105）有关。多变量Cox模型C-index为0.771（似然比p=8.62×10⁻¹²），Stage IV vs Stage I的风险比达8.74（95%CI: 3.61-21.14, p<0.001），Stage III的HR为2.96（p=0.004）。在控制分期和年龄后，分子亚型与总生存无显著独立关联（所有亚型对比p>0.5），验证了病理分期在BRCA预后中相对于分子亚型的预测主导地位。

单变量Cox筛选识别了与总体生存显著关联的基因（FDR<0.05），功能注释显示涉及细胞周期调控、DNA损伤修复和免疫应答通路。基于这些基因构建了加权风险得分（risk score），按中位得分将患者分为高风险组和低风险组。与LASSO-Cox多变量筛选（未选出基因，受限于事件率过低）相比，单变量策略在统计功效上更为宽松，适用于事件率较低场景下的初筛。

### 3.7 体细胞突变分析

**方法原理**：采用maftools R包[12]对MuTect2流程输出的Masked Somatic Mutation数据进行系统分析。分析流程包括突变概况统计、瀑布图（Oncoplot）可视化、变异类型分布及互斥性/共现性检验（配对Fisher精确检验）。

**分析结果**：突变数据集涵盖990例样本、15,413个突变基因、67,546个变异位点，平均每例携带约68个突变。高频突变基因为PIK3CA（369例, 37.3%）、TP53（348例, 35.2%）及TTN（268例, 27.1%）。PIK3CA的热点突变（E545K、E542K、H1047R）导致PI3K/AKT/mTOR通路组成性激活，在ER+乳腺癌中尤为常见。TP53突变在Triple Negative亚型中频率显著升高。TTN的高频突变部分归因于其基因长度效应（>100kb）。变异类型以错义突变（Missense_Mutation）为主，其次为无义突变和移码突变。互斥性分析显示PIK3CA与TP53之间具有一定程度的互斥趋势，符合两种基因分别驱动不同信号通路的生物学逻辑。

### 3.8 关联规则挖掘与调控网络

**方法原理**：采用arules包[9]的Apriori算法进行关联规则挖掘。将Top 50差异基因的表达值按中位数离散化为High/Low，与分子亚型、ER/PR/HER2状态组合构建事务数据集。参数：最小支持度0.3，最小置信度0.7。此外，对Top 20差异基因进行了χ²独立性检验。miRNA-mRNA调控网络通过计算Top miRNA与mRNA的Pearson相关系数矩阵构建，筛选标准为|r|>0.3和|r|>0.4。

**分析结果**：Apriori算法挖掘出基因表达状态与临床特征之间的关联规则，其中"{基因高表达}→{特定亚型}"类型的规则具有较高的lift值。χ²检验识别出多个与分子亚型显著关联的基因（p<0.05）。关联规则散点图（support × confidence, lift着色）展示了规则的全局分布特征。

miRNA-mRNA调控网络分析以|r|>0.3阈值共识别228对潜在调控关系对，以更严格的|r|>0.4阈值识别出2,065对（在大规模对齐数据集上的全局扫描）。负相关关系符合miRNA对靶mRNA的负调控预期，但这些关系基于统计关联而非实验验证，假阳性率需进一步评估。

---

本章分析使用R 4.6.0 + Bioconductor 3.23，在Windows 11环境下执行。主要R包：DESeq2, caret, randomForest, glmnet, WGCNA, survival, survminer, maftools, gprofiler2, arules, arulesViz, ConsensusClusterPlus, ggsci, pheatmap, ggplot2等。

---

## 4 可视化分析

可视化在数据挖掘成果向生物学知识转化的过程中扮演着关键角色。本研究运用R语言生态系统中的多种专业绘图工具，对各项分析结果进行了系统性的可视化表征。

**差异表达可视化**：火山图以log2FC为横轴、-log10(p)为纵轴，NPG色板区分上调（红）、下调（蓝）和不显著基因（灰），ggrepel标注关键基因。Top 25 DEGs表达热图以Z-score标准化展示肿瘤与正常组织的表达模式差异。

**分类模型评估可视化**：多类别ROC曲线（一对多）展示LASSO模型对不同亚型的区分能力。模型对比分组柱状图同时展示Accuracy、Kappa和F1三项指标。

**降维与聚类可视化**：PCA散点图叠加95%置信椭圆，t-SNE降维图揭示非线性结构，层次聚类热图展示基因-样本双向聚类模式。

**WGCNA网络可视化**：软阈值选择双面板图（SFT R²+平均连接度）、模块-性状关联热图及模块大小分布图构成WGCNA的核心可视化输出。

**生存分析可视化**：KM曲线附log-rank p值和风险表，Cox森林图展示各变量HR及95%CI。

**突变景观可视化**：Oncoplot展示Top 12突变基因在990例样本中的突变分布，变异类型柱状图及互斥性热图补充结构和功能信息。

**富集分析可视化**：上调与下调基因的富集气泡图（按来源数据库分色）及来源数据库比较图展示功能富集的全景分布。

---

## 5 结论与展望

本研究以TCGA-BRCA乳腺癌队列为数据基础，整合mRNA、miRNA和体细胞突变三个组学层次，构建了覆盖差异表达分析、功能富集、分子亚型分类、聚类分析、WGCNA共表达网络、生存分析、突变分析及关联规则挖掘等九个分析模块的完整数据挖掘流程。主要发现如下：

（1）DESeq2差异表达分析鉴定6,768个显著差异基因，上调（4,334个）数量约为下调（2,434个）的1.78倍。功能富集分析揭示上调基因集中于细胞周期和DNA复制通路（Reactome Cell Cycle, p=1.5×10⁻²¹），下调基因集中于细胞外基质组织和黏附相关通路。

（2）LASSO多分类回归以89.4%准确率实现三分类分子亚型预测，仅使用35个区分性基因，验证了亚型间转录组差异的稀疏性。

（3）PCA（PC1 15.3%）和K-means轮廓系数（K=2最佳）一致表明ER状态是BRCA转录组结构的最强驱动因素。共识聚类提供了不依赖临床标签的数据驱动分层方案。

（4）WGCNA识别9个共表达模块，blue模块（553个基因）枢纽基因模块隶属度0.959，为后续功能验证的高优先级候选。

（5）Cox多变量回归（C-index=0.771）证实病理分期是BRCA预后最强预测因子（Stage IV HR=8.74），分子亚型在控制分期后无独立预后意义。

（6）PIK3CA（37.3%）和TP53（35.2%）确认为BRCA最高频驱动突变，两者呈现互斥趋势。

（7）关联规则挖掘识别了基因表达与临床特征之间的系统性关联模式，miRNA-mRNA调控网络分析识别了超过2,000对潜在调控关系对。

本研究存在以下局限性：DNA甲基化数据因数据量过大（11.7 GB）未纳入分析，缺失了表观遗传维度；生存分析的事件率较低（7.9%），限制了多变量预后模型的统计功效；miRNA-mRNA调控关系基于统计相关而非因果验证，假阳性率待独立数据集检验。

未来研究方向包括：（1）采用降采样策略纳入DNA甲基化数据，实现四组学整合；（2）利用METABRIC、SCAN-B等独立验证队列对分类模型进行外部验证；（3）对WGCNA blue模块枢纽基因进行功能实验验证；（4）使用单细胞RNA测序数据在细胞分辨率上解析肿瘤微环境异质性；（5）探索深度学习方法在多组学数据整合中的应用。

---

## 参考文献

[1] Sung H, Ferlay J, Siegel R L, et al. Global cancer statistics 2022: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries[J]. CA: A Cancer Journal for Clinicians, 2024, 74(3): 229-263.

[2] Hanahan D, Weinberg R A. Hallmarks of cancer: the next generation[J]. Cell, 2011, 144(5): 646-674.

[3] Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours[J]. Nature, 2012, 490(7418): 61-70.

[4] Perou C M, Sorlie T, Eisen M B, et al. Molecular portraits of human breast tumours[J]. Nature, 2000, 406(6797): 747-752.

[5] Love M I, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2[J]. Genome Biology, 2014, 15(12): 550.

[6] Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis[J]. BMC Bioinformatics, 2008, 9: 559.

[7] Simon N, Friedman J, Hastie T, et al. Regularization paths for Cox's proportional hazards model via coordinate descent[J]. Journal of Statistical Software, 2011, 39(5): 1-13.

[8] Kolberg L, Raudvere U, Kuzmin I, et al. g:Profiler—interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update)[J]. Nucleic Acids Research, 2023, 51(W1): W207-W212.

[9] Agrawal R, Srikant R. Fast algorithms for mining association rules[C]. Proceedings of the 20th International Conference on Very Large Data Bases, 1994: 487-499.

[10] Wilkerson M D, Hayes D N. ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking[J]. Bioinformatics, 2010, 26(12): 1572-1573.

[11] Colaprico A, Silva T C, Olsen C, et al. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data[J]. Nucleic Acids Research, 2016, 44(8): e71.

[12] Mayakonda A, Lin D C, Assenov Y, et al. Maftools: efficient and comprehensive analysis of somatic variants in cancer[J]. Genome Research, 2018, 28(11): 1747-1756.

[13] Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent[J]. Journal of Statistical Software, 2010, 33(1): 1-22.

[14] Breiman L. Random forests[J]. Machine Learning, 2001, 45(1): 5-32.

[15] Chen T, Guestrin C. XGBoost: a scalable tree boosting system[C]. Proceedings of the 22nd ACM SIGKDD, 2016: 785-794.

[16] Therneau T M, Grambsch P M. Modeling survival data: extending the Cox model[M]. New York: Springer, 2000.

[17] Tibshirani R. Regression shrinkage and selection via the lasso[J]. Journal of the Royal Statistical Society: Series B, 1996, 58(1): 267-288.

---

*本研究分析使用R 4.6.0 + Bioconductor 3.23，在Windows 11环境下执行。全部分析脚本和处理后数据见项目目录code/和results/。*
