FindHVGUnion <- function(seu_list, nfeatures = 400) {
	hvg <- lapply(seu_list, function(xx){
		sublist <- SplitObject(xx, split.by = "samplename")
		smps <- setdiff(names(sublist), "CJB1680")
		subhvg <- lapply(smps, function(smp) {
			smpobj <- sublist[[smp]]
			## Remove lowly-expressed genes
			DefaultAssay(smpobj) <- "RNA"
			smpobj <- smpobj[Matrix::rowSums(smpobj$RNA@data != 0) >= 5, ]

			## Find HVGs
			smpobj <- FindVariableFeatures(smpobj, nfeatures = (nfeatures + 100))
			return(smpobj)
			}) %>%
					SelectIntegrationFeatures(., nfeatures = nfeatures)
		return(subhvg)
		}) %>%
			unlist() %>% 
			unique()
	return(hvg)
}




library(foreach)
library(parallel)

SingleDR <- function(object, hvg, do_UMAP = FALSE) {
	DefaultAssay(object) <- "integrated"
	object <- RunPCA(object, npcs = 10, verbose = FALSE, features = intersect(hvg, rownames(object$integrated@scale.data)))
	dr <- list(pca = object$pca)
	if (do_UMAP) {
		object <- RunUMAP(object, dims = 1:10, n.neighbors = 25, min.dist = 0.25, n.components = 1L, verbose = TRUE)
		dr$umap <- object$umap
	}
	return(dr)
}




BootDR <- function(object, hvg_list, nCores, do_UMAP = FALSE) {
	message(Sys.time())
    cl = makeCluster(nCores, type = "FORK", outfile="")
    doParallel::registerDoParallel(cl);

    res_list <- foreach(idx = 1:length(hvg_list), .inorder = TRUE, .export = c("SingleDR"), .packages = c("Seurat"), .verbose = TRUE) %dopar% {
    		print(paste0("Working on replicate: ", idx, " / ", length(hvg_list)))
    		dr <- SingleDR(object = object, hvg = hvg_list[[idx]], do_UMAP = do_UMAP)
			return(dr)
        }
    stopCluster(cl)
    message("Finish Bootstrap dimension reduction analysis")
    message(Sys.time())
	return(res_list)
}








