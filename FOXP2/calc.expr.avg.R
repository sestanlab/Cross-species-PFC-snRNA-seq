library(Seurat)
library(dplyr)
library(Matrix)
library(Signac)
library(ggplot2)


obj <- readRDS(file = paste0("~/myshare/AdultHumanMultiome/CombinedSeurat/Four.smps.adultPFC.addUMAP.v05092022.rds"))


ctdata <- obj$ARCRNA@counts
Idents(obj) <- "mres"
logavgs <- log(AverageExpression(object = obj, assay = "ARCRNA", slot = "data")$ARCRNA + 1) %>% 
			as.matrix()


identity <- as.character(obj@meta.data[, "mres"])
all_cls <- levels(as.factor(obj@meta.data[, "mres"]))
ratios <- lapply(all_cls, function(x) {
			message(paste0("Calculate expression ratio for cluster: ", x))	
			sub_cdt <- ctdata[, identity == x, drop = FALSE]
			er <- Matrix::rowMeans(sub_cdt != 0)
			er
			}) %>% setNames(., all_cls) %>%
				as.data.frame(., check.names = FALSE) %>% as.matrix()
save(logavgs, ratios, file = paste0("./load_files/", "Expr_avg_ratio_multiome_4sps.rds"))


