## Get the Average expression & expression ratio (composite & sep by species & sep by samplename)
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")


##------------------------------------------------------------------------
hprc <- readRDS(file = paste0(dataDir, "Final_PFC_HPRC.All.rm1680.07102021.rds"))
ii <- args[1] ##split.by c("species", "samplename", "all") 


group.by <- c("fig1cluster", "fig7cluster", "mres")
ctdata <- hprc$RNA@counts


res_list <- list()
for (jj in group.by){
	if (ii != "all"){
		hprc@meta.data$newcls <- paste0(as.character(hprc@meta.data[, ii]), "|", as.character(hprc@meta.data[, jj]))
		tt <- "newcls"
	} else {
		tt <- jj
	}


	Idents(hprc) <- tt
	logavgs <- log(AverageExpression(object = hprc, assay = "RNA", slot = "data")$RNA + 1) %>% as.matrix()


	identity <- as.character(hprc@meta.data[, tt])
	all_cls <- levels(as.factor(hprc@meta.data[, tt]))
	exp_ratio <- lapply(all_cls, function(x) {
			message(paste0("Calculate expression ratio for cluster: ", x))	
			sub_cdt <- ctdata[, identity == x, drop = FALSE]
			er <- Matrix::rowMeans(sub_cdt != 0)
			er
			}) %>% setNames(., all_cls) %>%
				as.data.frame(., check.names = FALSE) %>% as.matrix()
	res_list[[jj]] <- list(avg = logavgs, ratio = exp_ratio)
}
saveRDS(res_list, file = paste0(dataDir, "Expr_avg_ratio_by_", ii, ".rds"))







