# 第4章 可视化分析

数据挖掘做完了，东西跑出来了——但一堆数字和矩阵放在那里是没人看得懂的。可视化的任务就是把高维的、抽象的统计结果翻译成肉眼能读的信息。这一章逐个梳理了前面分析中用到的可视化策略：差异表达的火山图和热图、分类模型的ROC和混淆矩阵、聚类的PCA和t-SNE降维图、WGCNA的网络图和模块-性状热图、生存分析的KM曲线和森林图、突变分析的瀑布图——以及最后把它们拼成一幅综合汇总图。以下按这个顺序来。

## 4.1 差异表达可视化

### 4.1.1 火山图

火山图是差异表达分析里用得最多的全局视图。横轴log2FC、纵轴-log10(p)，每个基因一个点。两条阈值线把图切成几块——横轴±1（表达量变化2倍以上）、纵轴-log10(0.05)——落在右上和左上的就是显著差异基因。

BRCA肿瘤 vs 正常的火山图用ggplot2绘制，6,768个差异基因中4,334个标红（上调侧）、2,434个标蓝（下调侧）。一眼看过去，右侧的红色密度明显高于左侧——上调数量确实是下调的1.78倍。上调侧偏重不是偶然的：肿瘤细胞的转录组整体处于激活状态，癌基因、细胞周期基因、代谢重编程相关通路都在往上走。

有一件事值得提一句：火山图两侧那些log2FC特别大、p值特别小的极端离群点，往往就是关键驱动基因或候选生物标志物。后续功能验证可以优先从这些两极的点里挑。

### 4.1.2 表达热图

Top 50 DEG的聚类热图（pheatmap，Z-score标准化行）把肿瘤和正常分成了两个主要的样本分支——几乎没有混在一起的。这说明筛出来的差异基因确实能稳稳地把两组分开。颜色编码很直观：红色高表达、蓝色低表达，第一眼看过去就能判断哪边是肿瘤、哪边是对照。

## 4.2 分类模型性能可视化

### 4.2.1 多类ROC曲线

用pROC包为LASSO模型画了一对多（one-vs-all）的ROC曲线。每个亚型一条曲线，AUC值反映了模型对该亚型的区分能力。三条曲线都明显高于对角线（随机基线），LASSO对三类的区分能力在ROC空间里是可视化的——不是勉强过关，而是确实学到了东西。

### 4.2.2 混淆矩阵与模型对比

混淆矩阵以比例热图呈现——对角线是正确分类的比例，非对角线是误分类的模式。LASSO的混淆矩阵对角线占主导，误分类主要集中在Luminal A和HER2-enriched之间的边界——这两个亚型在表达谱上确实更接近。

模型对比的柱状图则很直白地把LASSO（89.4%）、RF（86.7%）和XGBoost（35.8%）排在一张图里。三个柱子一摆，哪个好用一目了然。XGBoost那根几乎比另外两根矮了一半——这个视觉冲击比任何文字描述都有说服力：参数不调，再好的模型也是废的。

## 4.3 聚类与降维可视化

### 4.3.1 PCA散点图

PCA散点图（PC1 15.3% × PC2 8.5%），四个分子亚型用NPG四种颜色区分，外加95%置信椭圆。Luminal A缩在最紧的一团里，Triple Negative散得最开——后者转录组异质性明显更强。HER2-enriched居中。PC1方向上ER+和ER-沿着整条轴拉开了距离。

### 4.3.2 t-SNE降维

t-SNE（perplexity=30, Rtsne包）在两维空间里把Luminal A和Triple Negative推得更远。HER2-enriched散布在中间地带——不偏不倚，恰好落在两种极端之间。t-SNE在捕捉局部结构上比PCA给了更多信息：在PCA里互相靠近但实际来自不同亚型的一些样本，在t-SNE里就被拉开了。

### 4.3.3 层次聚类热图

Top 500可变基因的双向热图用Ward.D2方法聚类。两大样本分支分别对应ER+和ER-——和PCA、t-SNE的结论一致。基因层面能看出细胞周期簇、激素信号簇和免疫应答簇各自成团。

## 4.4 WGCNA网络可视化

### 4.4.1 软阈值选择

双面板图：左半边是无标度拓扑拟合度随power的变化曲线（R²=0.887在power=8时过0.8线），右半边是平均连接度的衰减曲线。两个面板放在一起看，8是平衡点——R²够高，连接度也没掉到太低。

### 4.4.2 模块-性状关联热图

9个模块×临床性状的Pearson相关矩阵热图——颜色从蓝（负相关）到白（无关）到红（正相关）。blue模块和Luminal A之间的红色方块、brown模块和Triple Negative之间的正相关——这些都是后续功能验证要追的线索。

### 4.4.3 共表达网络

对blue模块导出了edge list（邻接矩阵>0.3阈值），枢纽基因在网络的中心位置一目了然。

## 4.5 生存分析可视化

### 4.5.1 KM生存曲线

按分期分组的KM曲线用survminer绘制——Stage I在最上面，Stage IV压在最底部，中间两条（II和III）有序排列。log-rank p值标在图里，风险表附在下方。按分子亚型分组的KM（p=0.493）四条线互相缠绕，几乎分不开——和分期的清晰分层形成了很直白的对比。

### 4.5.2 Cox森林图

每个变量的HR以水平误差条展示，虚线画在HR=1的位置。Stage IV那条误差条的HR点远远甩到了右边（约8.74），CI很宽但下限也远超1——哪怕是CI的下限，也意味着至少3.6倍的风险增加。年龄和分子亚型的HR点几乎压在无效线上。

## 4.6 突变景观可视化

### 4.6.1 瀑布图

maftools的oncoplot——Top 20基因在990例样本中的突变分布，用不同颜色区分错义突变、无义突变、移码突变等类型。PIK3CA那一行红色方块密密麻麻（369例），TP53紧随其后（348例）。瀑布图在突变分析中的地位就相当于火山图在差异分析中的地位——给一个"全局第一眼"。

### 4.6.2 变异类型和互斥性

变异类型分布柱状图——错义突变占压倒性多数，这是所有癌症体细胞突变的标准模式。互斥性热图展示了PIK3CA和TP53之间的互斥趋势——两个基因倾向于不在同一个肿瘤中共突变，支持它们驱动不同信号通路的生物学逻辑。

## 4.7 综合可视化汇总

### 4.7.1 六面板汇总图

用gridExtra把六张核心图拼成一幅20×12英寸的大图——火山图、亚型分布、分期分布、KM曲线、模型对比和WGCNA模块大小。布局遵循信息密度从高到低、从左到右、从上到下的自然阅读顺序。这张图单独拿出来就能作为整个项目的"一页纸总结"。

### 4.7.2 关键基因跨亚型热图

30个关键基因（含Top DEG和经典BRCA标志物如ESR1、PGR、ERBB2、MKI67、TP53、BRCA1、BRCA2等）按分子亚型和分期排列的热图，红蓝表达色阶——Luminal A的ER/PR高表达区和非Luminal型的低表达区形成清晰的分界线。

### 4.7.3 HTML交互报告

`brca_analysis_report.html`——把所有分析结果和数据卡片打包进一个浏览器能打开的网页。在不需要R环境的情况下，合作者或导师可以直接浏览完整的分析流程和结果。数据卡片（1,094 mRNA患者、6,768 DEGs、LASSO 89.4%、Cox C-index 0.771）以醒目的数字块呈现在报告顶部。

---

# 结论

本研究用TCGA-BRCA队列做了一次比较完整的多组学数据挖掘尝试。数据层面整合了三大块：mRNA转录组（25,981个基因，1,094例肿瘤）、miRNA表达谱（1,881个miRNAs，1,079例样本）和体细胞突变（990例样本，15,413个基因）。方法层面覆盖了差异表达、机器学习分类、无监督聚类、WGCNA共表达网络、生存分析和突变特征分析——基本把当前癌症生物信息学中常用的几类数据挖掘工具都过了一遍。

先说差异表达的结果。DESeq2在肿瘤和正常对照之间筛出了6,768个显著差异基因（|log2FC|>1, FDR<0.05），其中4,334个在肿瘤中上调、2,434个下调。上调占多数的全局趋势并不意外——癌细胞的转录组本身就处于一种广泛激活的状态。真正有筛选价值的是那些差异幅度特别大的"极端离群基因"，它们落在火山图的左右两端，构成后续功能验证和标志物筛选的最优先候选池。

分子亚型分类上，LASSO的表现明显优于另外两个模型。89.4%的准确率靠的是L1正则化从500个候选基因里自动挑出了35个真正有区分度的特征——在高维小样本场景下，这种自动变量筛选加同时建模的策略确实比随机森林（86.7%）更有效。XGBoost只拿到了35.8%，基本可以断定是默认参数的问题。这个对比本身也说明了一个实践教训：在癌症组学数据集上使用梯度提升方法，不做参数调优就不要指望结果能用。

无监督聚类的发现比较有意思。PCA就把事情讲得挺清楚了——PC1（15.3%）和PC2（8.5%）合起来，主要按ER状态把样本分成两大块，Luminal A缩在最紧的一个团里，Triple Negative则散得很开，说明后者的转录组异质性确实更明显。t-SNE进一步看到Luminal B和HER2-enriched像是横跨在两种极端之间的一道"过渡带"。K-means的最优K就是2——这和PCA的结论互相印证：ER状态就是BRCA转录组结构中最强的那根主轴，其他因素都在它下面。

WGCNA识别出了9个共表达模块，soft-threshold取8（scale-free R²=0.887）。最有意思的是blue模块：553个基因，枢纽基因的模块隶属度高达0.959——从这个模块里挑出来的前20个枢纽基因，后续做功能注释和实验验证的优先级应该排到最高。其他几个模块也和不同的临床特征呈现出了有规律的关联模式，不过这些关联的效应量偏小，独立验证的价值可能有限。

生存分析方面，Kaplan-Meier曲线给出的信息很直白：分期就是BRCA预后的最大boss。I期到IV期的曲线分离有序且显著，而分子亚型之间的差异反而没达到统计学显著性（log-rank p=0.493）。多变量Cox回归确认了这一点——C-index 0.771，Stage IV的风险比飙到约8.74。不过话说回来，这个队列87例死亡事件分布在1,105个患者里——事件率不到8%，中位随访也才2.7年。事件太少、随访太短，LASSO-Cox自然选不出基因来（结果是零），分子亚型的KM曲线也测不到差异——这些都是统计功效不足的问题，不是生物学上真没有差异。

突变分析的结果和已知的BRCA驱动突变谱高度吻合。PIK3CA（369例）和TP53（348例）稳坐前两名——一个是PI3K通路的激活，一个是基因组监护的崩坏。TTN排第三很大程度上是基因长度带来的背景噪音——这一点在用maftools做突变频次排序时一定要心中有数。前20个高频基因的互斥性和共现性分析额外提供了一些信号通路层级关系的线索。

多组学整合这一层，miRNA和mRNA的关联分析找到了2,065对相关系数超过0.4的配对关系——数量不算少，但Pearson相关毕竟不是因果推断，这里面假阳性的比例有多高，得靠独立数据集来验证。1,076例患者同时拥有mRNA和miRNA数据，这个重叠度（98.4%）本身说明TCGA的数据结构对多组学分析是比较友好的。

总结一下这个项目的方法论价值：它跑的其实是一条比较常规的分析管线——从数据下载到预处理到算法评估再到可视化——但每一步的代码、参数和输出都是公开和可复现的。13个R脚本串起了整个流程，中间踩过的坑（subtype分类的运算符优先级bug、XGBoost默认参数的翻车、甲基化数据大到没法下载）反而比那些顺利跑通的部分更有参考意义。

至于局限性，有几件事需要老老实实承认：甲基化数据因为11.7 GB的体积实在太大，下载和加载都超出了当前环境的承受范围，最终没有纳入——缺失了表观遗传这个维度，这是最遗憾的一点。预后模型方面，事件率太低直接导致LASSO-Cox没选出基因，所以本文没有一个可以拿去外部验证的多基因风险评分。XGBoost的调参工作也没来得及系统做。miRNA-mRNA的调控关系目前还停留在相关层面。

后面可以推的方向包括：(1) 用降采样或启动子探针子集的方式，把甲基化数据加进来，实现真正意义上的四组学整合；(2) 把分类模型拿到METABRIC或SCAN-B这样的外部数据集上去跑一跑，看看89.4%的准确率到底能不能复现；(3) 从blue模块的枢纽基因里挑几个做功能验证——CRISPR敲除或者过表达之后测一下表型变化，这是最直接的下一步；(4) 如果有条件拿到BRCA的单细胞数据，就可以把bulk水平的分子亚型拆解到细胞类型的分辨率上去；(5) 深度学习方法在处理多组学数据融合时可能比LASSO和随机森林有更多的发挥空间，但这个方向需要更大的样本量和更扎实的计算资源。

## 参考文献

[1] Hanahan D, Weinberg R A. Hallmarks of cancer: the next generation[J]. Cell, 2011, 144(5): 646-674.

[2] Sung H, Ferlay J, Siegel R L, et al. Global cancer statistics 2022: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries[J]. CA: A Cancer Journal for Clinicians, 2024, 74(3): 229-263.

[3] Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours[J]. Nature, 2012, 490(7418): 61-70.

[4] Perou C M, Sorlie T, Eisen M B, et al. Molecular portraits of human breast tumours[J]. Nature, 2000, 406(6797): 747-752.

[5] Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis[J]. BMC Bioinformatics, 2008, 9: 559.

[6] Love M I, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2[J]. Genome Biology, 2014, 15(12): 550.

[7] Tibshirani R. Regression shrinkage and selection via the lasso[J]. Journal of the Royal Statistical Society: Series B, 1996, 58(1): 267-288.

[8] Mayakonda A, Lin D C, Assenov Y, et al. Maftools: efficient and comprehensive analysis of somatic variants in cancer[J]. Genome Research, 2018, 28(11): 1747-1756.

[9] Colaprico A, Silva T C, Olsen C, et al. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data[J]. Nucleic Acids Research, 2016, 44(8): e71.

[10] Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent[J]. Journal of Statistical Software, 2010, 33(1): 1-22.

[11] Breiman L. Random forests[J]. Machine Learning, 2001, 45(1): 5-32.

[12] Chen T, Guestrin C. XGBoost: a scalable tree boosting system[C]. Proceedings of the 22nd ACM SIGKDD, 2016: 785-794.

[13] Therneau T M, Grambsch P M. Modeling survival data: extending the Cox model[M]. New York: Springer, 2000.

[14] Simon N, Friedman J, Hastie T, et al. Regularization paths for Cox's proportional hazards model via coordinate descent[J]. Journal of Statistical Software, 2011, 39(5): 1-13.
