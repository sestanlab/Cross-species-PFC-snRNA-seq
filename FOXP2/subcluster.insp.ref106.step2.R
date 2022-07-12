source("../scripts/pfc.fun.R")
source("./inte.fun.R")



seu <- readRDS(file = "./load_files/Multiome-subcluster_insp_Immune.rds")

## Find clusters (remove outliers)
seu <- FindNeighbors(seu, dims = 1:40) %>%
                        FindClusters(., resolution = 1)
DimFig(seu, group.by = "seurat_clusters", file_name = paste0("Multiome-subcluster_insp_Immune"), plot.scale = 0.8)

##res1 <- FindMarkers(seu, ident.1 = "6", only.pos = TRUE, max.cells.per.ident = 500) # 6 is Oligo
##res2 <- FindMarkers(seu, ident.1 = "7", only.pos = TRUE, max.cells.per.ident = 500) # 7 is Vas
rm_cells <- colnames(seu)[as.character(seu$seurat_clusters) %in% c("6", "7") | seu$samplename == "HSB106"]
subseu <- seu[, setdiff(colnames(seu), rm_cells)]
rm(seu)




library(future)
plan("multiprocess", workers = 6)
options(future.globals.maxSize = 30*1000*1024^2)


bat1 <- readRDS(file = paste0("../recluster_v2/load_files/Cls_bysp_detail_", "Immune", "_Human.rds"))
DefaultAssay(bat1) <- 'RNA'

hvg <- bat1$integrated %>%
		rownames() %>% intersect(., rownames(subseu))
mars <- readRDS(file = paste0("../SpecclsMars/load_files/", "Immune.speccls.putative.markers.rds"))[c("Human CCL3", "Human B", "Human GLDN-v2")] %>%
			unlist() %>% unique()
hvg <- union(hvg, mars) %>%
		intersect(., rownames(subseu))


bat1 <- subset(bat1, samplename %in% c('HSB106', 'HSB189', 'HSB340'))
seu <- merge(x = subseu, y = bat1)
DefaultAssay(seu) <- 'RNA'


cbn <- SeuInte(object = seu, hvg = hvg, file_name = paste0("Multiome-subcluster_realign_Immune"), input_dir = inputdir, do_cluster = FALSE, merge.order = NULL)








