# CUT-RUN-snakemake
CUT &amp; RUN analysis 

### CUT & RUN原理
![cut_run](https://github.com/SitaoZ/CUT-RUN-snakemake/assets/29169319/052e8504-6a29-4e65-9b79-9978efed0016)

CUT & RUN (Cleavage Under Targets and Release Using Nuclease)是一种体内的方法，使用抗体和蛋白的特异性结合获取特定的蛋白和DNA复合体进行下游分析。
根据其原理，CUT & RUN可以用于转录因子绑定序列的识别，RBP蛋白绑定序列的识别，聚集小体等蛋白的识别。
CUT&RUN实验周期短，高通量化(CUT&RUN-Seq)可以很好的替代ChIP-Seq技术。技术的原始文献为 [cut & run](https://pubmed.ncbi.nlm.nih.gov/28079019/)。
###  优点
- 需要更少的上样量。ChIP需要500K cells, 而CUT&RUN最少可低至5K cells。
- 测序深度要求不高。CUT&RUN 每个样本测3-8M read/sample,而ChIP-Seq则需要大于等于30M read/sample。
- 背景噪音显著减少，因为CUN&RUN是目标片段测序。
- 成本减少。因为抗体的使用减少，建库测序的深度减少。
### 缺点
- 由于钙依赖性MNase反应时间不合适，可能导致DNA过度消化。类似的限制存在于现在的ChIP-seq操作中，酶或超声DNA打断必须优化。
- 并不是所有的蛋白质都适用于该方案。根据感兴趣的蛋白质，可能需要摸索实验条件。
- 可能是染色质复合物太大而无法扩散出去，或者蛋白质-蛋白质相互作用保留了被切割的复合物。在这种情况下，需要消化后才能提取总DNA。
### 分析步骤
CUT & RUN 适合分析转录因子结合位点

### 注意事项

### Reference
[Henikoff lab paper](https://elifesciences.org/articles/21856)
[CUT & RUN参考](https://nf-co.re/cutandrun/3.2/docs/output#6-peak-calling)
