## Harmony integration of all cells
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")



##---------------------------------------------------------------------------------------------
## Function to perform the harmony integration
library(harmony)
InteAllSp.harmony <- function(object, split.by, hvg, file_name, input_dir = inputdir, inte.dims = 30, theta = 2, lambda = 1, sigma = 0.1) {

    inte.slim.file <- paste0(input_dir, file_name, ".harmony.rds")
    if (!file.exists(inte.slim.file)){
        ## Do the Harmony integration
        object <- ScaleData(object, split.by = split.by, do.center = FALSE, features = hvg)%>%
                    RunPCA(., features = hvg, verbose = FALSE) %>%
                    RunHarmony(., group.by.vars = split.by, lambda = lambda, theta = theta, dims.use = 1:inte.dims, sigma = sigma)
        object <- RunUMAP(object, dims = 1:inte.dims, reduction = "harmony", umap.method = "umap-learn", metric = "correlation")
        object <- FindNeighbors(object, dims = 1:inte.dims, reduction = "harmony", k.param = 25) %>%
                    FindClusters(., resolution = 1.2, n.iter = 20)
        saveRDS(object, file = inte.slim.file)
    }else {
        object <- readRDS(file = inte.slim.file)
    }
    return(object)
}



##---------------------------------------------------------------------------------------------
## Use the above function to perform harmony integration
ctp <- args[1] ## ExN, InN, NNC
pfc <- readRDS(file = paste0(dataDir, "Final_PFC_HPRC.", ctp, ".rm1680.07102021.rds"))



## Get the highly variable genes
hvg.list <- readRDS(file = paste0(dataDir, "HVG.list.", ctp, ".rds"))
sp_infor <- list(Human = c("HSB628", "HSB189", "HSB340", "HSB106"),
				Chimpanzee = c("PTB165", "PTB166", "PTB2169", "PTB1841"),
				Rhesus = c("RMB161", "RMB196", "RMB295", "RMB307"),
				Marmoset = c("CJB1435", "CJB1540", "CJB1577"))
hvg <- lapply(sp_infor, function(x) SelectHVG(hvg.list[x], nfeatures = 3000)) %>%
			unlist() %>%
			unique()


seu <- InteAllSp.harmony(object = pfc, split.by = "samplename", hvg = hvg, file_name = paste0("Recls_all_", ctp), input_dir = inputdir, inte.dims = 40, theta = 2, lambda = 0.8, sigma = 0.1)














