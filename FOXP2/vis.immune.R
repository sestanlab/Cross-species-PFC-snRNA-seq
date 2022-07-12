source("../scripts/pfc.fun.R")

seu <- readRDS(file = "./load_files/Multiome-subcluster_realign_Immune.rds")
seu <- FindNeighbors(seu, dims = 1:40) %>%
                        FindClusters(., resolution = 1.6)
seu$cluster <- get_ident(input_ident = seu, ident_col = "seurat_clusters", file_path = paste0(inputdir, "Immune.anno.txt"), sample_name = "Immune", label_names = NULL)



subseu <- subset(seu, samplename %in% c("2RT00374N", "RT00382N", "RT00383N", "RT00390N"))
subseu$samplename <- setNames(extract_field(colnames(subseu), 1, "_"), NULL)
subseu <- RunUMAP(subseu, dims = 1:40, umap.method = "umap-learn", metric = "correlation")
subseu <- FindNeighbors(subseu, dims = 1:40) %>%
                        FindClusters(., resolution = 1)
saveRDS(subseu, file = "./load_files/Multiome-subcluster_realign_Immune_filtered.rds")




## Load micro markers
mars <- readRDS(file = paste0("../SpecclsMars/load_files/", "Immune.speccls.putative.markers.rds"))
bat1 <- subset(seu, samplename %in% c("HSB106", "HSB189", "HSB340", "HSB628"))
Idents(bat1) <- "fig1cluster"
res <- FindMarkers(bat1, ident.1 = "T SKAP1 CD247", max.cells.per.ident = 500, only.pos = TRUE)
mars$Tcells <- rownames(res)[res$p_val_adj <= 0.01 & res$avg_logFC >= 0.5 & res$pct.1 >= 0.3 & res$pct.2 <= 0.3]



auc <- GetModuleScore(assay.data = subseu$RNA@data, features = mars, method = "aucell", input_dir = inputdir, file_name = paste0("Multiome-subcluster_realign_Immune_filtered"), output_dir = outputdir, rethreshold_list = NULL, cellbin.size = 8000)$auc


pdata <- subseu@meta.data[, c("samplename", "seurat_clusters"), drop = FALSE] %>%
			cbind(., auc[colnames(subseu), c("Human CCL3", "Human GLDN-v2", "Tcells", "Human B"), drop = FALSE]) %>%
			setNames(., c("donors", "seurat_clusters", "Micro_CCL3", "Micro_GLDN", "Tcells", "Bcells")) %>%
			cbind(., subseu$umap@cell.embeddings) %>%
			as.data.frame(., check.names = FALSE)
library(ggrastr)

## Random order
set.seed(42)
pdata <- pdata[sample(1:nrow(pdata)), ]


p1 <- ggplot(pdata, aes_string(x = "UMAP_1", y = "UMAP_2", color = "donors")) +
				rasterise(geom_point(size = 0.4, shape = 16), dpi = 300, scale = 1) +
	            coord_equal(ratio = 1) + 
	            theme_classic() + 
	            theme(legend.position = "bottom",
	                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
	                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
	                    plot.title = element_blank())
exp_list <- lapply(c("Micro_CCL3", "Micro_GLDN", "Tcells", "Bcells"), function(gg) {
	max_exp <- max(pdata[, gg]) * 1
	print(max_exp)
	pdata[, gg] <- MinMax(pdata[, gg], min = 0, max = max_exp)

	p <- ggplot(pdata, aes_string(x = "UMAP_1", y = "UMAP_2", color = gg)) +
				rasterise(geom_point(size = 0.4, shape = 16), dpi = 300, scale = 1) +
	            coord_equal(ratio = 1) + 
	            theme_classic() + 
	            scale_color_gradientn(colors = c("#BEBEBE", "#C4AAAA", "#E35050", "#FF0000"), limits = c(0, max_exp)) + 
	            theme(legend.position = "bottom",
	                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
	                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
	                    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
	return(p)
	})


pdf(paste0(outputdir, "SF.Multiome.Immune.pdf"), width = 4 * 5, height = 5.4)
plot <- patchwork::wrap_plots(plotlist = c(list(p1), exp_list), nrow = 1, ncol = 5) #& theme(legend.position = "bottom")
print(plot)
dev.off()






## Plot the cell type proportion across donors
meta <- readRDS(file = "../data/PFC.filtered.meta.05082022.rds") %>%
			filter(mres %in% c("Micro", "Immune") & species == "Human")
size_bat1 <- meta %>%
			mutate(fig1cluster = as.character(fig1cluster)) %>%
			group_by(samplename, fig1cluster) %>%
			summarize(ncells = n()) %>%
			ungroup() %>%
			group_by(samplename) %>%
			mutate(ratio = ncells * 100/sum(ncells), total_cells = sum(ncells)) %>%
			ungroup() %>%
			filter(fig1cluster %in% c("Micro P2RY12 CCL3", "T SKAP1 CD247", "Micro P2RY12 GLDN")) %>%
			mutate(cluster = fig1cluster) %>%
			select(samplename, ratio, total_cells, cluster) %>%
			mutate(batch = "RNA")
subseu@meta.data$cluster <- as.character(subseu@meta.data$seurat_clusters)
subseu@meta.data$cluster[subseu@meta.data$seurat_clusters %in% c("6", "8")] <- "T SKAP1 CD247"
subseu@meta.data$cluster[subseu@meta.data$seurat_clusters %in% c("7")] <- "Micro P2RY12 CCL3"
subseu@meta.data$cluster[subseu@meta.data$seurat_clusters %in% c("5")] <- "Micro P2RY12 GLDN"

size_bat2 <- subseu@meta.data %>%
			group_by(samplename, cluster, .drop = FALSE) %>%
			summarize(ncells = n()) %>%
			ungroup() %>%
			group_by(samplename) %>%
			mutate(ratio = ncells * 100/sum(ncells), total_cells = sum(ncells)) %>%
			ungroup() %>%
			filter(cluster %in% c("Micro P2RY12 CCL3", "T SKAP1 CD247", "Micro P2RY12 GLDN")) %>%
			select(samplename, ratio, total_cells, cluster) %>%
			mutate(batch = "Multiome")
sp_ord <- c(size_bat1$samplename, size_bat2$samplename) %>% unique()
size_data <- rbind(size_bat1, size_bat2) %>%
			tidyr::gather(., "type", "value", c("ratio", "total_cells"))
size_data$newtype <- ifelse(size_data$type == "ratio", paste0(size_data$type, "|", size_data$cluster), "total_cells")
size_data <- size_data %>%
			group_by(samplename, newtype) %>%
			summarize(value = mean(value), batch = unique(batch)) %>%
			mutate(samplename = factor(as.character(samplename), levels = sp_ord))

p <- ggplot(size_data, aes_string(x = "samplename", y = "value", fill = "batch")) +
			geom_bar(stat = "identity", size = 0.2) +
			theme_classic() +
			RotatedAxis() +
			facet_wrap(vars(newtype), nrow = 4, ncol = 1, scales = "free_y") +
			theme(strip.background = element_blank(),
				axis.line = element_line(size = 0.2),
				axis.ticks = element_line(size = 0.2),
				axis.text = element_text(size = rel(0.8)))
pdf(paste0(outputdir, "SF.Immune.ratio.twobatches.pdf"), width = 4, height = 4)
print(p)
dev.off()









