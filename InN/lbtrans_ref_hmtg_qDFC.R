## This script will generate the data for sankey plot visualization
source("../scripts/pfc.fun.R")
source("../scripts/inte.fun.R")


library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 28*1000*1024^2)


## Do the label transfer between allen and our data
hmtg <- readRDS(file = paste0(dataDir, "Human_MTG_org_05072021.rds")) %>%
			subset(class == "GABAergic")
pfc <- readRDS(file = paste0(dataDir, "Final_PFC_HPRC.InN.rm1680.07102021.rds"))
hdfc <- subset(pfc, species == "Human")


## Add some information
hmtg@meta.data$samplename <- "HSB000"
hmtg@meta.data$dataset <- "HMTG"

hdfc@meta.data$dataset <- "HDFC"


source("../scripts/inte.fun.R")

inte.list <- Inte.Prepare(data.list = list(MTG = hmtg, PFC = hdfc), inte.by = "samplename", group.by = c(MTG = "cluster", PFC = "fig1cluster"), split.by = "dataset", hvg.method = c("exphvg", "share")[1], nfeatures = 2000)


## Transfer the mouse annotation to each of the human samples
QuickTransfer <- function(ref_obj, query_obj, hvg, inte.dims = 1:30, k.filter = 200, k.score = 30, k.weight = 50, group.by = "cluster") {
	raw_ident <- ref_obj[[group.by]]
    seu.anchors <- FindTransferAnchors(reference = ref_obj, query = query_obj, dims = inte.dims, reduction = "cca", k.filter = k.filter, k.score = k.score, max.features = 300, features = hvg)
    predictions <- TransferData(anchorset = seu.anchors, weight.reduction = "cca", k.weight = k.weight, refdata = setNames(raw_ident[, 1], rownames(raw_ident)), dims = inte.dims)
	query_obj <- AddMetaData(query_obj, metadata = predictions)
	return(query_obj)
}



hsb_list <- inte.list$data[c("HSB106", "HSB189", "HSB340", "HSB628")]
meta <- lapply(hsb_list, function(x) {
	seu <- QuickTransfer(ref_obj = inte.list$data[["HSB000"]], query_obj = x, inte.dims = 1:30, k.filter = 200, k.score = 30, k.weight = 50, group.by = "cluster", hvg = inte.list$hvg)
	submeta <- seu@meta.data[, c("fig1cluster", "predicted.id", "samplename")]
	submeta
	}) %>%
		do.call(rbind, .)
saveRDS(meta, file = paste0(inputdir, "HDFC_transed_by_HMTG.meta.rds"))







