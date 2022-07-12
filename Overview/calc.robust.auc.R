args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")
source("./robust.fun.R")
library(ROCR)



sp_list <- list(Human = c("HSB106", "HSB189", "HSB340", "HSB628"), 
				Chimpanzee = c("PTB165", "PTB166", "PTB1841", "PTB2169"), 
				Rhesus = c("RMB161", "RMB196", "RMB295", "RMB307"), 
				Marmoset = c("CJB1435", "CJB1540", "CJB1577", "CJB1680"))


sp <- args[1]
ctp <- args[2]
smps <- sp_list[[sp]]



aucFile <- paste0(inputdir, "RobustAUC_", sp, "_", ctp, ".rds")
if (!file.exists(aucFile)){
	meta <- readRDS(file = paste0(dataDir, "PFC.filtered.meta.05082022.rds"))


	auc_list <- list()
	for (smp in smps){
		print(paste0("Calculating robustness: ", sp, "::", ctp, " remove ", smp))
		seu <- readRDS(file = paste0(inputdir, "ThreeDonorInte/ThreeSmpInte_", sp, "_", ctp, "_rm_", smp, ".rds"))

		pcs <- seu$pca@cell.embeddings
		idents <- setNames(as.character(meta[colnames(seu), "fig1cluster"]), colnames(seu))
		auc_list[[smp]] <- subcell_robustness(dr = pcs, identity = idents, dims = 40, k = 15)
	}


	df <- lapply(auc_list, function(x) x[names(auc_list[[1]])]) %>%
	        setNames(., smps) %>%
	        as.data.frame(., check.names = FALSE) %>%
	        as.matrix()
	saveRDS(df, file = aucFile)
}









