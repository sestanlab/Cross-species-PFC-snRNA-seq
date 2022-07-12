## Find the HVG
source("../scripts/pfc.fun.R")


pfc <- readRDS(file = paste0(dataDir, "PFC_filtered_seu.11052020.rds"))


for (ctp in c("ExN", "InN", "NNC")){
	subseu <- subset(pfc, group == ctp)

	hvg.list <- SplitObject(subseu, split.by = "samplename") %>%
					lapply(., function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)) %>%
					lapply(., function(x) VariableFeatures(x))
	saveRDS(hvg.list, file = paste0(inputdir, "HVG.byCTP.", ctp, ".rds"))
}




## Select the top HVG for each major cell type and get the union
sp_list <- list(Human = c("HSB106", "HSB189", "HSB340", "HSB628"),
				Chimpanzee = c("PTB165", "PTB166", "PTB1841", "PTB2169"), 
				Rhesus = c("RMB161", "RMB196", "RMB295", "RMB307"),
				Marmoset = c("CJB1435", "CJB1540", "CJB1577")) ## not include CJB1680

flist <- list()
for (nfeatures in c(1000, 1250, 1500)){
	final.hvg <- lapply(c("ExN", "InN", "NNC"), function(ctp) {
		hvg.list <- readRDS(file = paste0(inputdir, "HVG.byCTP.", ctp, ".rds"))


		new.hvg <- lapply(names(sp_list), function(sp) 
			SelectHVG(hvg.list = hvg.list[sp_list[[sp]]], nfeatures = nfeatures)
			) %>% 
				unlist() %>% unique() %>% unname()
		return(new.hvg)
		}) %>% 
			unlist() %>% unique() %>% unname()

	flist[[as.character(nfeatures)]] <- final.hvg
}
saveRDS(flist, file = paste0(inputdir, "HVG.byCTP.combine.rds"))








