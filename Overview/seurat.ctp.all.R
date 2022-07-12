## This script will do the Seurat integration on  all HPRC cells. 

source("./pfc.fun.R")
library(future)
plan("multiprocess", workers = 6)
options(future.globals.maxSize = 150*1000*1024^2)


##------------------------------------------------------------------------
## Integrate all samples 
inte.file <- paste0(dataDir, "Integrate.hprc.full.rds")
anchor.file <- paste0(dataDir, "Integrate.hprc.anchor.rds")

hprc <- readRDS(file = paste0(dataDir, "PFC_filtered_seu.11052020.rds"))
hprc <- hprc[, hprc@meta.data$samplename != "CJB1680"] 
## CJB1680 has a relatively larger batch effect and has small amount of cells. This sample will be later project to all cells based on the transcriptomic similarity with other marmoset cells.


if (!file.exists(inte.file)){
	if (!file.exists(anchor.file)){
		hprc_list <- SplitObject(hprc, split.by = "samplename") %>%
					lapply(., function(x) FindVariableFeatures(x, nfeatures = 2000))

    	data.anchors <- FindIntegrationAnchors(object.list = hprc_list, dims = 1:30, assay = NULL, anchor.features = 2500)
    	rm(hprc_list)
    	saveRDS(data.anchors, file = anchor.file)
	} else {
		data.anchors <- readRDS(file = anchor.file)
	}


    hprc.inte <- IntegrateData(anchorset = data.anchors, dims = 1:30)
    DefaultAssay(hprc.inte) <- "integrated"
    hprc.inte <- ScaleData(hprc.inte, verbose = FALSE) %>%
                                    RunPCA(., npcs = 40, verbose = FALSE)
    saveRDS(hprc.inte, file = inte.file)
}


##------------------------------------------------------------------------
## Make the data slim (remove the "integration assay)
inte.slim.file <- paste0(dataDir, "Integrate.hprc.seurat.slim.rds")
if (!file.exists(inte.slim.file)){
    hprc.inte <- readRDS(file = inte.file)
    
    hprc[["pca"]] <- CreateDimReducObject(embeddings = hprc.inte$pca@cell.embeddings[colnames(hprc), ], loadings = hprc.inte$pca@feature.loadings, stdev = hprc.inte$pca@stdev, key = "PC_", assay = "RNA")

    hprc <- RunUMAP(hprc, dims = 1:40, umap.method = "umap-learn", metric = "correlation")
    saveRDS(hprc, file = inte.slim.file)

    DimFig(hprc, group.by = c("fig1cluster", "mres", "species"), file_name = "HPRC", plot.scale = 2)
}           




















