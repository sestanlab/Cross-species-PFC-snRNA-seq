args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")
source("./preprocess.fun.R")


all_samples <- c("2RT00374N", "RT00382N", "RT00383N", "RT00385N", "RT00390N")
gp <- args[1]


seu <- readRDS(file = paste0(inputdir, "InteAll_withHSB106_", gp, ".filtered.rds"))
seu_list <- SplitObject(seu, split.by = "samplename")



load(file = paste0(inputdir, "Human_HVGs.Rdata"))
## hvg_all, hvg_exn, hvg_inn, hvg_nnc
hvg <- switch(gp, 
			ExN = hvg_exn,
			InN = hvg_inn,
			NNC = hvg_nnc)


## Perfrom integration
meta <- lapply(all_samples, function(smp) {
	print(paste0("Working on sample: ", smp))
	subseu <- QuickTransfer(ref_obj = seu_list[["HSB106"]], query_obj = seu_list[[smp]], inte.dims = 1:30, group.by = "mres", hvg = hvg)
	submeta <- subseu@meta.data[, c("mres", "predicted.id", "samplename")]
	submeta
	}) %>%
		do.call(rbind, .)
saveRDS(meta, file = paste0(inputdir, "Label_transferred_", gp, "_meta.rds"))



