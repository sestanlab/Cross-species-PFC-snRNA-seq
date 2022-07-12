## Perform data integration
FastIntePCA <- function(object, split.by = "samplename", nfeatures = 1000){
	hprc_list <- SplitObject(object, split.by = "samplename")
	data.anchors <- FindIntegrationAnchors(object.list = hprc_list, dims = 1:30, assay = NULL, anchor.features = nfeatures)
    rm(hprc_list)
 
	hprc.inte <- IntegrateData(anchorset = data.anchors, dims = 1:30)
    DefaultAssay(hprc.inte) <- "integrated"
	hprc.inte <- ScaleData(hprc.inte, verbose = FALSE) %>%
                        RunPCA(., npcs = 40, verbose = FALSE)


    object[["pca"]] <- CreateDimReducObject(embeddings = hprc.inte$pca@cell.embeddings[colnames(object), ], loadings = hprc.inte$pca@feature.loadings, stdev = hprc.inte$pca@stdev, key = "PC_", assay = "RNA")
    object <- RunUMAP(object, dims = 1:30)
    return(object)
}




##-------------------------------------------------------
## Calculate AUC scores
neighbor_voting <- function(dr, identity, dims.use = 30, group.by, k = 20){
    #Find the k-nearst neighbor of each inquiry cell in the references datase via the input_dr dimensions
    data.use <- ifelse_check(length(dims.use) == 1, dr[, 1:dims.use], dr[, dims.use])

    get_most <- function(x){
        return(names(sort(table(x), decreasing = TRUE))[1])
    }


    if (sum(is.na(identity)) > 0){
    	##Assume the NA values are the inquiry cell
	    ref_cells <- names(identity)[!is.na(identity)]
	    inquiry_cells <- setdiff(names(identity), ref_cells)


	    ## Get the KNN cells for all the inquiry cells
	    knn_cells <- FNN::get.knnx(data = data.use[ref_cells, ,drop = FALSE], query = data.use[inquiry_cells, ,drop = FALSE], k = k) %>%
	                    .$nn.index
	    ref_anno <- identity[ref_cells] %>% setNames(., ref_cells)

	    ## Do the annotation transfer
	    new_ident <- apply(knn_cells, 1, function(x) get_most(ref_anno[x])) %>% 
	                            setNames(., inquiry_cells) %>%
	                            c(., ref_anno) %>%
	                            .[names(identity)]
	} else {
		knndata <- FNN::get.knn(data = data.use, k = k)
	    knn_idx <- knndata$nn.index;
	    rm(knndata)

	    new_ident <- setNames(apply(knn_idx, 1, function(x) get_most(identity[x])), names(identity))
	}
    return(new_ident) 
}



subcell_robustness <- function(dr, identity, dims = 40, k = 15) {
	## Check the consistency of cell order
	if (sum(names(identity) %in% rownames(dr)) != length(identity)){
		stop("Cell numbers are inconsistent between the reduced dimensions and identity")
	}
	dr <- dr[names(identity), ,drop = FALSE]


	pred_ident <- neighbor_voting(dr = dr, identity = identity, dims.use = dims, k = k)
    all_cls <- levels(as.factor(as.character(identity)))


    simi <- sapply(all_cls, function(x) {
        pred <- prediction(ifelse(pred_ident == x, 1, 0), 
                                    ifelse(identity == x, 1, 0))
        perf <- performance(pred, measure = "auc")@y.values
        return(perf[[1]])
        })
    return(simi)
}







