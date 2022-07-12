## Calculate PCA for observed data & bootstrap replicates (with replacement)
## Output
## 1. raw data
## 2. subseted data (balanced cells across species)
## 3. HVGs, bootstrap replicates of HVGs 
## 4. PCA


args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")
source("./calc.bootstrap.fun.R")
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

    ## Do the scaling using hvg genes
    DefaultAssay(seu) <- "integrated"
    seu <- ScaleData(seu, features = intersect(hvg_use, rownames(seu$integrated@data)), vars.to.regress = c("percent.mt", "nUMI"))
    return(seu)
    }) %>%
        setNames(., all_sps)



## Subset (balance cell numbers across species)
## Determine the number of cells to use
msize <- sapply(seu_list, ncol) %>% min()
if (msize > 1000){
	ncells <- 1000
} else if (msize <= 1000 & msize > 500){
	ncells <- 500
} else if (msize <= 500 & msize > 200){
	ncells <- 200
} else {
	ncells <- 100
}



set.seed(0)
## Prepare cell list, 100 replicates
full_list <- lapply(seu_list, colnames)
cell_list <- replicate(100, lapply(full_list, function(x) sample(x, ncells)), simplify = FALSE)


## Bootstrap replicates of HVG genes
set.seed(42)
nsp_size <- ceiling(length(hvg_use) * 0.8)
hvg_boot <- replicate(100, sample(hvg_use, nsp_size, replace = FALSE), simplify = FALSE)


res_list <- list()
for (ii in 1:length(cell_list)){
	print(paste0("Working on iteration: ", ii, " / ", length(cell_list)))
	cur_list <- cell_list[[ii]]
	seu_list2 <- lapply(names(cur_list), function(sp) {
		seu <- seu_list[[sp]][, cur_list[[sp]]]
		DefaultAssay(seu) <- "integrated"
	    return(seu)
	    }) %>%
	        setNames(., all_sps)
	

	res_list[[ii]] <- lapply(seu_list2, function(seu) {
			boot_res <- BootDR(object = seu, hvg_list = hvg_boot, nCores = 12, do_UMAP = FALSE)
			raw_res <- SingleDR(object = seu, hvg = hvg_use, do_UMAP = FALSE)
			return(list(boot = boot_res, raw = raw_res))
			})
}

saveRDS(res_list, file = paste0(inputdir, "Bootstrap_PCA_v1/Bootstrap_DR_", cls_use, "_res.rds"))

