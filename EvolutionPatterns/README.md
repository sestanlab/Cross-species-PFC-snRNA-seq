---
title: Overview of cellular and molecular changes across the primate dlPFC 
author: Shaojie Ma
date: July 11, 2022
---


#### Evaluate species-specificity in this study
plot.auc.arhgap18.R \
plot.auc.astro.R \
plot.auc.immune.R \
plot.auc.syt10.R \


#### Plot the expression of representative markers of species-specific subtypes
plot.marexpr.vln.arhgap18.R \
plot.marexpr.vln.rAstro.R \
plot.marexpr.vln.syt10.R \


#### Evaluate species-specificity in motor cortex
plot.syt10.motor.dlPFC.R \
plot.arhgap18.motor.dlPFC.R \


#### Cell type proportion comparisons (scCODA)
- Prepare data for scCODA analysis \
prep_sccoda_v1.R \

- scCODA analysis \
calc_sccoda_v2_defref.py \
plot.sccoda.hres.fdr02.v2.R \

- Plot pie summary for scCODA results \
plot.sccoda.summary.R \

- Plot subtype proportions \
plot.prop.R \

- Plot main figure network plot \
source functions: \
network.fun.v2.R \

- visualization: \
plot.MF2B.R \


#### Transcriptomic heterogeneity
- Entropy analysis \
calc.heter.fun.v3.R \
vis.heter.fun.R \
calc.bin.entropy.v4.R \
vis.dif.entropy.v4.R \

- Bootstrap PCA analysis
calc.bootstrap.fun.R \
calc.bootstrap.pca.R \
calc.boot.meandist.v1.R \
vis.boot.pca.v1.R \


#### Transcriptomic divergence
- Transcriptomic divergence \
Different types: balanced HVGs from each class; top HVGs across all cells; subset cell \
get_hvg_byctp.R \
get_hvg_noctp.R \
calc.subset.avg.R \
plot.cor.bal.ctp.v2.R \
plot.cor.nobal.v2.R \
plot.cor.subset.v2.R \


- Plot # DEGs \
plot.ndex.R \







