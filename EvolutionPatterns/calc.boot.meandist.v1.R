args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")


## Load bootstrap Dimension reduction results

cls_ord <- c("L2_3_IT", "L3_5_IT_1","L3_5_IT_2","L3_5_IT_3", "L5_PT", "L5_6_NP", "L6_CT", "L6_IT_1", "L6_IT_2", "L6b",
        "LAMP5_LHX6", "LAMP5_RELN", "VIP", "ADARB2", "SST", "TH", "PVALB","PVALB_CHC",
        "Astro", "OPC","Oligo", "Micro","immune", "blood","Endo","PC","SMC","VLMC") ##"SST_NPY", 
all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")


cls_use <- cls_ord[as.numeric(args[1])]



spres_list <- readRDS(file = paste0(inputdir, "Bootstrap_PCA_v1/Bootstrap_DR_", cls_use, "_res.rds"))
boot_list <- lapply(spres_list, function(xx) {
		lapply(xx, function(yy) yy$boot)
		})


## Calculate mean euclidean distance (for bootstraps)
boot_dis <- lapply(1:length(boot_list), function(idx) {
		print(paste0("Working on replicate ", idx, " / ", length(boot_list)))
		sub_list <- boot_list[[idx]]
		sub_dis <- lapply(names(sub_list), function(sp) {
				adis <- sapply(sub_list[[sp]], function(xx) { 
					mdis <- xx$pca@cell.embeddings[, 1:10] %>%
							dist() %>% 
							mean()
					return(mdis)
					})
				df <- data.frame(idx = idx, 
						meanDis = adis, 
						species = sp, 
						cluster = cls_use,
						type = "bootstrap",
						stringsAsFactors = FALSE)
		return(df)
		}) %>%
		do.call(rbind, .)

		return(sub_dis)
		}) %>%
		do.call(rbind, .)




## Calculate mean euclidean distance (for raw data)
raw_list <- lapply(spres_list, function(xx) {
		lapply(xx, function(yy) yy$raw)
		})
raw_dis <- lapply(1:length(boot_list), function(idx) {
		sub_list <- raw_list[[idx]]
		sub_dis <- lapply(names(sub_list), function(sp) {
				mdis <- sub_list[[sp]]$pca@cell.embeddings[, 1:10] %>%
							dist() %>% 
							mean()
				df <- data.frame(idx = idx, 
						meanDis = mdis, 
						species = sp, 
						cluster = cls_use,
						type = "observed",
						stringsAsFactors = FALSE)
		return(df)
		}) %>%
		do.call(rbind, .)


		return(sub_dis)
		}) %>%
		do.call(rbind, .)


final_df <- rbind(boot_dis, raw_dis)

saveRDS(final_df, file = paste0(inputdir, "Bootstrap_PCA_v1/Bootstrap_PCA-MeanDist_", cls_use, "_res.rds"))










