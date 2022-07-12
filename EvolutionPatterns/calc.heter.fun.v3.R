CalcUMAPbinAvg <- function(object, dims = 1:10, nbins, features) {
	DefaultAssay(object) <- "integrated"
	object <- RunPCA(object, npcs = 30, verbose = FALSE, features = rownames(object$integrated@scale.data)) %>%
				RunUMAP(., dims = dims, n.neighbors = 25, min.dist = 0.25, n.components = 1L, verbose = FALSE)
	qs <- quantile(object$umap@cell.embeddings[, 1], c(0.01, 0.99))
	pdata <- data.frame(UMAP = object$umap@cell.embeddings[, 1], 
					cell = colnames(object),
					stringsAsFactors = FALSE) %>%
				filter(UMAP <= qs[2] & UMAP >= qs[1]) %>%
				mutate(UMAPbin = as.character(as.numeric(cut(UMAP, nbins))))

	## Remove bins with less than 5 cells to avoid variability
	kp_bins <- table(pdata$UMAPbin) %>% .[. >= 5] %>% names()
	pdata <- filter(pdata, UMAPbin %in% kp_bins)


	subobj <- object[, pdata$cell]
	DefaultAssay(subobj) <- "RNA"
	subobj$UMAPbin <- pdata$UMAPbin
	Idents(subobj) <- "UMAPbin"
	avgs <- log(AverageExpression(subobj, features = features, assays = "RNA", verbose = FALSE)$RNA + 1) %>%
			as.matrix() %>%
			.[features, ,drop = FALSE]
	return(avgs)
}



CalcEntropyDiffs_allspecies <- function(avg_res, ncores = 8L, nbins = 15) {
	nperms <- length(avg_res$perm)
	enp_res <- mclapply(1:nperms, function(idx) {
		enp <- c(avg_res$raw, avg_res$perm[[idx]]) %>%
				CalcEntropyDiffs(avg_list = ., nbins = nbins)
		return(enp)
		}, mc.cores = ncores)
	return(enp_res)
}




CalcUMAPbinAvg_allspecies <- function(object_list, dims, nbins, features, verbose = FALSE, npermutations = 8L, ncores = 8L) {
	raw_list <- lapply(names(object_list), function(sp) {
		if (verbose) {
			print(paste0("Averaging raw data in species: ", sp))
		}
		
		raw <- CalcUMAPbinAvg(object = object_list[[sp]], dims = dims, nbins = nbins, features = features)
		return(raw)
		}) %>%
			setNames(., paste0("raw_", names(object_list)))


	perm_reps <- mclapply(1:npermutations, function(idx) {
		perm_list <- lapply(names(object_list), function(sp) {
			obj <- object_list[[sp]]

			if (verbose) {
				print(paste0("Averaging permutated data in species: ", sp))
			}
			set.seed(idx)
			obj$integrated@scale.data <- apply(obj$integrated@scale.data, 1, sample) %>% t() 
			per <- CalcUMAPbinAvg(object = obj, dims = dims, nbins = nbins, features = features)
			return(per)
			}) %>%
			setNames(., paste0("perm_", names(object_list)))
		return(perm_list)
		}, mc.cores = ncores)
	return(list(raw = raw_list, perm = perm_reps))
}






CalcEntropy <- function(explist, nbins = 10) {
	xx <- unlist(explist)
	if (min(xx) == max(xx)){
		enp_list <- sapply(explist, function(x) 0)
	} else {
		bks <- seq(min(xx), max(xx), length.out = nbins + 1)
		enp_list <- sapply(explist, function(x1) {
			y1 <- as.numeric(cut(x1, bks, include.lowest = TRUE))
			fq1 <- table(y1)/length(y1)
			enp1 <- -sum(sapply(fq1, function(x) x*log(x)))
			return(enp1)
			})
	}
	
	return(enp_list)
}




CalcEntropyDiffs <- function(avg_list, nbins = 20) {
	all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
	features <- rownames(avg_list[[1]])
	enp_res <- lapply(features, function(gg) {
		enp <- lapply(avg_list, function(xx) xx[gg, ]) %>%
					CalcEntropy(explist = ., nbins = nbins)
		res <- sapply(all_sps, function(sp){
			enp[paste0("raw_", sp)] - enp[paste0("perm_", sp)]
			}) %>%
				setNames(., all_sps)
		return(res)
		}) %>%
			setNames(., features) %>%
			as.data.frame(., check.names = FALSE) %>%
			t()
	return(enp_res)
}





