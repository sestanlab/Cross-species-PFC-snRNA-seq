## Plot the number of DEX genes for each clusters 
source("../scripts/pfc.fun.R")
library(ggpubr)


## Downsampling the dataset and get the avgs (by species)
pfc <- readRDS(file = paste0(dataDir, "Final_PFC_HPRC.All.07102021.rds"))


## save <= 100 cells per cluster
pfc@meta.data$spcls <- paste0(as.character(pfc@meta.data$species), "|", as.character(pfc@meta.data$fig1cluster))


cls_order <- readRDS(file = paste0(dataDir, "PFC.cluster.order.rds"))
cls_order <- setdiff(cls_order, c("L2-3 CUX2 ARHGAP18", "InN LAMP5 SYT10", "rAstro AQP4 OSMR","nOligo MOG CDH7", "Micro P2RY12 CCL3", "Micro P2RY12 GLDN", "B EBF1 IGKC", "cOPC PDGFRA TOP2A")) ##"ABC SLC47A1 SLC4A4", 
cls_order <- cls_order[c(1:13, 15:21, 14, 22:length(cls_order))]


#homo_cls <- readRDS(file = paste0(dataDir, "homolog.cluster.rds"))
all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
all_cls <- rep(all_sps, each = length(cls_order)) %>%
				paste0(., "|", rep(cls_order, times = 4))


set.seed(20210804)
cells <- lapply(all_cls, function(x) {
	subc <- colnames(pfc)[pfc@meta.data$spcls == x]
	if (length(subc) >= 100){
		subc <- sample(subc, 100)
	}
	return(subc)
	}) %>%
		unlist()


seu <- pfc[, cells]
Idents(seu) <- "spcls"
avgs <- log(AverageExpression(seu)$RNA + 1)
saveRDS(avgs, file = paste0(inputdir, "Subset.avg.rds"))























