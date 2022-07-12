args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")
source("./inte.fun.R")


sp <- args[1]

sprawFile <- paste0(inputdir, "RawSeurat_", sp, ".rds")
if (!file.exists(sprawFile)){
	meta <- readRDS(file = paste0(dataDir, "PFC.filtered.meta.05082022.rds"))
	meta <- meta[meta$species == sp, ]
	reps <- levels(as.factor(meta$repname))
	genes <- readRDS(file = paste0(inputdir, "filtered_counts/Filtered.", reps[1], ".rds")) %>%
				rownames()
	counts <- lapply(reps, function(x) {
		print(paste0("Working on replicate: ", x))
		ctm <- readRDS(file = paste0(inputdir, "filtered_counts/Filtered.", x, ".rds"))
		return(ctm[genes, ])
		}) %>%
			do.call(cbind, .)

	if (sum(rownames(meta) %in% colnames(counts)) != nrow(meta)){
		stop(paste0("Cell numbers are inconsistent for species: ", sp))
	}


	seu <- seu_prepare(counts = counts, min.cells = 0, normalization.method = "LogNormalize", hvg.method = NULL)
	sel_cols <- setdiff(colnames(meta), colnames(seu@meta.data))
	seu@meta.data <- cbind(seu@meta.data, meta[colnames(seu), sel_cols])
	seu@meta.data$fig1cluster <- as.character(seu@meta.data$fig1cluster)
	saveRDS(seu, file = sprawFile)
}



