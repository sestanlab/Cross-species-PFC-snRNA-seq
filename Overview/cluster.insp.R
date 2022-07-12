args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")
source("./inte.fun.R")


library(future)
plan("multiprocess", workers = 2)
options(future.globals.maxSize = 90*1000*1024^2)


sp <- args[1]

seu <- readRDS(file = paste0(inputdir, "RawSeurat_", sp, ".rds"))


## Prepare HVGs
seu_list <- SplitObject(seu, split.by = "samplename")
hvg <- lapply(seu_list[setdiff(names(seu_list), "CJB1680")], function(x) FindVariableFeatures(x, nfeatures = 1500)) %>%
				SelectIntegrationFeatures(., nfeatures = 3000)
## CJB1680 has small sample size and was not included in HVG selection
rm(seu_list)


cbn <- SeuInte(object = seu, hvg = hvg, file_name = paste0("Cluster_insp_", sp), input_dir = inputdir, do_cluster = FALSE, merge.order = NULL)


