args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")
source("../realign/inte.fun.R")



library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 45*1000*1024^2)


sp <- args[1]

seu <- readRDS(file = paste0(dataDir, "Final_PFC_HPRC.All.03022022.rds"))
meta <- readRDS(file = paste0(dataDir, "PFC.filtered.meta.05082022.rds"))
meta <- meta[meta$species == sp & meta$group %in% c("Glia", "NNC"), ]

seu <- seu[, rownames(meta)]
seu@meta.data$fig1cluster <- as.character(meta[colnames(seu), "fig1cluster"])
seu@meta.data$mres <- as.character(meta[colnames(seu), "mres"])


## Prepare HVGs
seu_list <- SplitObject(seu, split.by = "samplename")
hvg <- lapply(seu_list[setdiff(names(seu_list), "CJB1680")], function(x) FindVariableFeatures(x, nfeatures = 1500)) %>%
				SelectIntegrationFeatures(., nfeatures = 2000) ## CJB1680 has limited cells and was not included for HVG calculation
rm(seu_list)


cbn <- SeuInte(object = seu, hvg = hvg, file_name = paste0("Cluster_NNC_insp_", sp), input_dir = inputdir, do_cluster = FALSE, merge.order = NULL)


DimFig(cbn, group.by = c("fig1cluster", "mres", "samplename"), file_name = paste0("Cluster_NNC_insp_", sp))





