## Find the HVG
source("../scripts/pfc.fun.R")


if (FALSE){
	pfc <- readRDS(file = paste0(dataDir, "PFC_filtered_seu.11052020.rds"))

	hvg.list <- SplitObject(pfc, split.by = "samplename") %>%
						lapply(., function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 4500, verbose = FALSE)) %>%
						lapply(., function(x) VariableFeatures(x))
	saveRDS(hvg.list, file = paste0(inputdir, "HVG.noCTP.rds"))
}



## Select the top HVG for each major cell type and get the union
sp_list <- list(Human = c("HSB106", "HSB189", "HSB340", "HSB628"),
				Chimpanzee = c("PTB165", "PTB166", "PTB1841", "PTB2169"), 
				Rhesus = c("RMB161", "RMB196", "RMB295", "RMB307"),
				Marmoset = c("CJB1435", "CJB1540", "CJB1577"))

flist <- list()
for (nfeatures in c(2000, 2500, 3000, 3500, 4000, 4500)){
	hvg.list <- readRDS(file = paste0(inputdir, "HVG.noCTP.rds"))
	final.hvg <- lapply(names(sp_list), function(sp) 
			SelectHVG(hvg.list = hvg.list[sp_list[[sp]]], nfeatures = nfeatures)
			) %>% 
			unlist() %>% unique() %>% unname()
	flist[[as.character(nfeatures)]] <- final.hvg
}
saveRDS(final.hvg, file = paste0(inputdir, "HVG.noCTP.combine.rds"))








