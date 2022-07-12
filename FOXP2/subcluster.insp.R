args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")
source("../realign/inte.fun.R")


library(future)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 18*1000*1024^2)



gp <- args[1]
nfeatures <- 1000
sel_gps <- switch(gp,
				Astro = "Astro", 
				Immune = c("immune", "Micro"))
seu <- readRDS(file = paste0("../Multiome_preprocess/load_files/", "InteAll_FINAL.rds"))
seu <- subset(seu, mres %in% sel_gps)
if (gp == "Immune"){
	seu@meta.data$raw_samplename <- seu@meta.data$samplename
	seu@meta.data$samplename <- gsub("RT00385N", "RT00390N", seu@meta.data$samplename) ## RT00385N has limited microglia, merge it with another sample for integration
}



match_gp <- switch(gp,
				Astro = "Astro", 
				Immune = "Immune")
hvg <- readRDS(file = paste0("../recluster_v2/load_files/Cls_bysp_detail_", match_gp, "_Human.rds"))$integrated %>%
		rownames() %>% intersect(., rownames(seu))



cbn <- SeuInte(object = seu, hvg = hvg, file_name = paste0("Multiome-subcluster_insp_", gp), input_dir = inputdir, do_cluster = FALSE, merge.order = NULL)

