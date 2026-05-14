# 结论

本研究用TCGA-BRCA队列做了一次比较完整的多组学数据挖掘尝试。数据层面整合了三大块：mRNA转录组（25,981个基因，1,094例肿瘤）、miRNA表达谱（1,881个miRNAs，1,079例样本）和体细胞突变（990例样本，15,413个基因）。方法层面覆盖了差异表达、机器学习分类、无监督聚类、WGCNA共表达网络、生存分析和突变特征分析——基本把当前癌症生物信息学中常用的几类数据挖掘工具都过了一遍。

先说差异表达的结果。DESeq2在肿瘤和正常对照之间筛出了6,768个显著差异基因（|log2FC|>1, FDR<0.05），其中4,334个在肿瘤中上调、2,434个下调。上调占多数的全局趋势并不意外——癌细胞的转录组本身就处于一种广泛激活的状态。真正有筛选价值的是那些差异幅度特别大的"极端离群基因"，它们落在火山图的左右两端，构成后续功能验证和标志物筛选的最优先候选池。

分子亚型分类上，LASSO的表现明显优于另外两个模型。89.4%的准确率（Kappa=0.528）靠的是L1正则化从500个候选基因里自动挑出了35个真正有区分度的特征——在高维小样本场景下，这种"自动变量筛选+同时建模"的策略确实比随机森林（86.7%）更有效。XGBoost只拿到了35.8%，基本可以断定是默认参数的问题，不是模型本身不行。这个对比本身也说明了一个实践教训：在癌症组学数据集上使用梯度提升方法，不做参数调优就不要指望结果能用。

无监督聚类的发现比较有意思。PCA就把事情讲得挺清楚了——PC1（15.3%）和PC2（8.5%）合起来，主要按ER状态把样本分成两大块，Luminal A缩在最紧的一个团里，Triple Negative则散得很开，说明后者的转录组异质性确实更明显。t-SNE进一步看到Luminal B和HER2-enriched像是横跨在两种极端之间的一道"过渡带"。K-means的最优K就是2——这和PCA的结论互相印证，ER状态就是BRCA转录组结构中最强的那根主轴，其他因素都在它下面。

WGCNA识别出了9个共表达模块，soft-threshold取8（scale-free R²=0.887）。最有意思的是blue模块：553个基因，枢纽基因的模块隶属度高达0.959——从这个模块里挑出来的前20个枢纽基因，后续做功能注释和实验验证的优先级应该排到最高。其他几个模块（turquoise 584个基因、brown 336个基因、yellow 260个基因）也和不同的临床特征呈现出了有规律的关联模式，不过这些关联的效应量偏小，独立验证的价值可能有限。

生存分析方面，Kaplan-Meier曲线给出的信息很直白：分期就是BRCA预后的最大boss。I期到IV期的曲线分离有序且显著，而分子亚型之间的差异反而没达到统计学显著性（log-rank p=0.493）。多变量Cox回归确认了这一点——C-index 0.771，Stage III和Stage IV的风险比分别飙到约5.7和约30.4，相比之下分子亚型的作用几乎被稀释干净了。不过话说回来，这个队列87例死亡事件分布在1,105个患者里——事件率不到8%，中位随访也才2.7年。事件太少、随访太短，LASSO-Cox自然选不出基因来（结果是零），分子亚型的KM曲线也测不到差异——这些都是统计功效不足的问题，不是生物学上真没有差异。

突变分析的结果和已知的BRCA驱动突变谱高度吻合。PIK3CA（369例）和TP53（348例）稳坐前两名——一个是PI3K通路的激活，一个是基因组监护的崩坏，这两件事在BRCA的分子发病机制里就是两条主线。TTN排第三（268例）很大程度上是基因长度带来的背景噪音——这一点在用maftools做突变频次排序时一定要心中有数。前20个高频基因的互斥性和共现性分析额外提供了一些信号通路层级关系的线索，不过这些计算证据需要更精细的实验来坐实。

多组学整合这一层，miRNA和mRNA的关联分析找到了2,065对相关系数超过0.4的配对关系——数量不算少，但Pearson相关毕竟不是因果推断，这里面假阳性的比例有多高，得靠独立数据集来验证。1,076例患者同时拥有mRNA和miRNA数据，这个重叠度（98.4%）本身说明TCGA的数据结构对多组学分析还是比较友好的。

总结一下这个项目的方法论价值：它跑的其实是一条比较常规的分析管线——从数据下载到预处理到算法评估再到可视化——但每一步的代码、参数和输出都是公开和可复现的。13个R脚本串起了整个流程，中间踩过的坑（subtype分类的运算符优先级bug、XGBoost默认参数的翻车、甲基化数据大到没法下载）反而比那些顺利跑通的部分更有参考意义。

至于局限性，有几件事需要老老实实承认：甲基化数据因为11.7 GB的体积实在太大，下载和加载都超出了当前环境的承受范围，最终没有纳入——缺失了表观遗传这个维度，这是最遗憾的一点。预后模型方面，事件率太低直接导致LASSO-Cox没选出基因，所以本文没有一个可以拿去外部验证的多基因风险评分。XGBoost的调参工作也没来得及系统做。miRNA-mRNA的调控关系目前还停留在相关层面。

后面可以推的方向包括：(1) 用降采样或启动子探针子集的方式，把甲基化数据加进来，实现真正意义上的四组学整合；(2) 把分类模型拿到METABRIC或SCAN-B这样的外部数据集上去跑一跑，看看89.4%的准确率到底能不能复现；(3) 从blue模块的枢纽基因里挑几个做功能验证——CRISPR敲除或者过表达之后测一下表型变化，这是最直接的下一步；(4) 如果有条件拿到BRCA的单细胞数据，就可以把bulk水平的分子亚型拆解到细胞类型的分辨率上去；(5) 深度学习方法（自编码器、图神经网络）在处理多组学数据融合时可能比LASSO和随机森林有更多的发挥空间，但这个方向需要更大的样本量和更扎实的计算资源。

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
