---
title: Data Overview
author: Shaojie Ma
date: July 11, 2022
---


#### Generate dendrogram
mf1.dendro.gen.v3.R \
It also requires subtype markers, as calculated in this script: mar.calc.R; zfun.wilcox.R

#### Do circular plot
mf1.plot.v8.R \
It also requires subtype average expression, as calculated in this script: avg.calc.R


#### Plot data quality
- nGenes, nUMIs, percent.mito, etc: \
sf1.quality.v3.R


- Sample contribution: \
sf1.sample.contribution.R


#### Calculate cluster robustness
- Source functions: \
robust.fun.R \

- Randomly remove one sample and integrate the rest three samples \
inte.3donor.exn.R \
inte.3donor.inn.R \
inte.3donor.nnc.R \

- Then, Calculate AUC scores evaluating cluster separability followed by visualization \
calc.robust.auc.R \
plot.robust.auc.R \

- Also checked the # markers per subtype \
plot.robust.mar.R


#### Alternative integration methods
- Harmony \
harmony.ctp.all.R

- MNN \
mnn.ctp.all.R

- Seurat \
seurat.ctp.all.R



#### Robustness to different gene annotation models
- Extract count data \
extract_sp.R \


- Cluster in each species and visualize \
cluster.insp.R \
cluster.insp.vis.R \


- Visualize species-specific subtypes \
subcluster.insp.R \
subcluster.insp.cbn.vis.R \


- Visualize species-specific expression \
plot.dem.exn.R \
plot.dem.inn.R \
plot.dem.nnc.R \
plot.th.umap.R \
plot.foxp2.R \












