## mnn integration of all cells
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")



##---------------------------------------------------------------------------------------------
## Function to perform the MNN integration
RunMNNCustom <- function(object, split.by, assay = "RNA", hvg, k, ndims = 40, merge.order = NULL) {
    ## Split data by batches
    seu_list <- SplitObject(object, split.by = split.by)


    ## Do MNN
    data <- lapply(seu_list, function(xx) xx[[assay]]@data[hvg,])
    mnn_res <- do.call(batchelor::fastMNN, c(data, list(k = k, cos.norm = TRUE, d = ndims, merge.order = merge.order, subset.row = NULL)))
    sub_mnn <- reducedDim(mnn_res)

    ## Modify row and col names
    rownames(sub_mnn) <- lapply(data, colnames) %>% unlist(., use.names = FALSE)
    colnames(sub_mnn) <- paste0("MNN_", as.character(1:ndims))
    return(sub_mnn)
}



InteAllSp.mnn <- function(object, split.by, hvg, cls.dims = 35, do.cluster = TRUE, k =16, merge.order = NULL, input_dir = inputdir, file_name) {
    inte.slim.file <- paste0(input_dir, file_name, ".mnn.rds")
    if (!file.exists(inte.slim.file)){
        ## Perform MNN dimension reduction
        sub_mnn <- RunMNNCustom(object = object, split.by = split.by, assay = "RNA", hvg = hvg, k = k, ndims = 50, merge.order = merge.order)
        object[["mnn"]] <- CreateDimReducObject(embeddings = sub_mnn, key = "MNN_", assay = "RNA")

    
        ## RunUMAP and cluster
        object <- RunUMAP(object, dims = 1:cls.dims, reduction = "mnn", umap.method = "umap-learn", metric = "correlation")
        object <- FindNeighbors(object, dims = 1:cls.dims, reduction = "mnn", k.param = 25) %>%
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


seu <- InteAllSp.mnn(object = pfc, split.by = "samplename", hvg = hvg, file_name = paste0("Recls_all_", ctp), input_dir = inputdir, cls.dims = 40, k =16, do.cluster = TRUE, merge.order = sp_infor)














