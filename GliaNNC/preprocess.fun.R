library(stringr) 
library(dplyr)
library(tibble)
library(Seurat)
library(parallel)
library(ggplot2)


Integratelist.seurat <- function(obj.list, hvg, file_name, input_dir = inputdir, inte.dims = 1:30, cluster.dims = 1:30, reference = NULL, do.cluster = FALSE) {
    if (length(obj.list) == 2){
        newseu <- merge(x = obj.list[[1]], y = obj.list[[2]])
    } else {
        newseu <- merge(x = obj.list[[1]], y = obj.list[2:length(obj.list)])
    }


    inte.slim.file <- paste0(input_dir, file_name, ".slim.rds")
    if (!file.exists(inte.slim.file)){
        dg.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = inte.dims, assay = NULL, anchor.features = hvg, reference = reference)
        seuinte <- IntegrateData(anchorset = dg.anchors, dims = inte.dims)
        DefaultAssay(seuinte) <- "integrated"
        seuinte <- ScaleData(seuinte, verbose = FALSE) %>%
                                RunPCA(., npcs = 50, verbose = FALSE)
        
        newseu[["pca"]] <- CreateDimReducObject(embeddings = seuinte$pca@cell.embeddings[colnames(newseu), ], loadings = seuinte$pca@feature.loadings, stdev = seuinte$pca@stdev, key = "PC_", assay = "RNA")
        rm(seuinte)
        newseu <- RunUMAP(newseu, dims = cluster.dims, umap.method = "umap-learn", metric = "correlation")
        if (do.cluster){
            newseu <- FindNeighbors(newseu, dims = cluster.dims,  k.param = 25) %>%
                    FindClusters(., resolution = 1.2, n.iter = 20)
        } else {
            newseu@meta.data$seurat_clusters <- "empty"
        }
        saveRDS(newseu, file = inte.slim.file)
    }else {
        newseu <- readRDS(file = inte.slim.file)
    }
    return(newseu)
}



