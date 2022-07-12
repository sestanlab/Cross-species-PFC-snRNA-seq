args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")
source("./calc.heter.fun.v3.R")
library(tidyr)



## Prepare cell type name
cls_list <- list(L23 = "L2_3_IT",
    L35 = c("L3_5_IT_1","L3_5_IT_2","L3_5_IT_3"),
    L56 = c("L5_PT", "L5_6_NP", "L6_CT", "L6_IT_1", "L6_IT_2", "L6b"),
    CGE = c("LAMP5_LHX6", "LAMP5_RELN", "VIP", "ADARB2"), 
    MGE = c("SST", "TH", "PVALB","PVALB_CHC"),
    Astro = c("Astro"),
    Oligo = c("OPC","Oligo"),
    Immune = c("Micro","immune"),
    Vas = c("blood","Endo","PC","SMC","VLMC")) ## SST_NPY is not included due to low number of cells
hres2lres <- lapply(names(cls_list), function(gp) setNames(rep(gp, length(cls_list[[gp]])), cls_list[[gp]])) %>%
            do.call(c, .)


cls_use <- names(hres2lres)[as.numeric(args[1])]
gp <- hres2lres[cls_use]



## Load seurat objects
all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
hvg_use <- readRDS("../recluster_v4/load_files/Subclass_hvgs_union_species.rds")[[cls_use]]
seu_list <- lapply(all_sps, function(sp) {
    seu <- readRDS(paste0("../recluster_v4/load_files/", "Cls_bysp_detail_", gp, "_", sp, ".rds"))
    ## Subset to the given cell type
    seu <- subset(seu, mres == cls_use)

    ## Regress out nUMI, percent.mt for the integrated data
    DefaultAssay(seu) <- "integrated"
    seu <- ScaleData(seu, features = intersect(hvg_use, rownames(seu$integrated@data)), vars.to.regress = c("percent.mt", "nUMI"))
    return(seu)
    }) %>%
        setNames(., all_sps)



## Determine the number of cells in each subset
msize <- sapply(seu_list, ncol) %>% min()
if (msize > 1000){
	ncells <- 1000
	bin_par <- c(20, 15)
} else if (msize <= 1000 & msize > 500){
	ncells <- 500
	bin_par <- c(15, 10)
} else if (msize <= 500 & msize > 200){
	ncells <- 200
	bin_par <- c(10, 7)
} else {
	ncells <- 100
	bin_par <- c(7, 5)
}



set.seed(0)
## Prepare cell list, 100 replicates
full_list <- lapply(seu_list, colnames)
cell_list <- replicate(100, lapply(full_list, function(x) sample(x, ncells)), simplify = FALSE)



res_list <- list()
for (ii in 1:length(cell_list)){
	print(paste0("Working on iteration: ", ii, " / ", length(cell_list)))
	cur_list <- cell_list[[ii]]
	seu_list2 <- lapply(names(cur_list), function(sp) {
		seu <- seu_list[[sp]][, cur_list[[sp]]]
	    return(seu)
	    }) %>%
	        setNames(., all_sps)

	## Get average list
	avg_res <- CalcUMAPbinAvg_allspecies(object_list = seu_list2, dims = 1:10, nbins = bin_par[1], features = hvg_use, verbose = FALSE, npermutations = 30L, ncores = 8L)
	res_list[[ii]] <- CalcEntropyDiffs_allspecies(avg_res = avg_res, ncores = 8L, nbins = bin_par[2])
}
saveRDS(res_list, file = paste0(inputdir, "UMAP_Bin_Res_v4/Entropy_list_regressed_", cls_use, ".v3.rds"))






