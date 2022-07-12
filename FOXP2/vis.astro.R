source("../scripts/pfc.fun.R")

seu <- readRDS(file = "./load_files/Multiome-subcluster_insp_Astro.rds")


## Find clusters (remove outliers)
seu <- FindNeighbors(seu, dims = 1:40) %>%
                        FindClusters(., resolution = 1)
DimFig(seu, group.by = "seurat_clusters", file_name = paste0("Multiome-subcluster_insp_Astro"), plot.scale = 0.8)


rm_cells <- colnames(seu)[seu$samplename == "HSB628" | as.character(seu$seurat_clusters) %in% c("8")]
subseu <- seu[, setdiff(colnames(seu), rm_cells)]
subseu <- RunUMAP(subseu, dims = 1:40, umap.method = "umap-learn", metric = "correlation")
DimFig(subseu, group.by = "seurat_clusters", file_name = paste0("Multiome-subcluster_insp_Astro_filtered"), plot.scale = 0.8)



## Load rAstro markers
mars <- readRDS(file = paste0("../SpecclsMars/load_files/", "Astro.speccls.putative.markers.rds"))
auc <- GetModuleScore(assay.data = subseu$RNA@data, features = mars, method = "aucell", input_dir = inputdir, file_name = paste0("Multiome-subcluster_insp_Astro_filtered"), output_dir = outputdir, rethreshold_list = NULL, cellbin.size = 8000)$auc



subseu <- FindNeighbors(subseu, dims = 1:40) %>%
                        FindClusters(., resolution = 1)
saveRDS(subseu, file = "./load_files/Multiome-subcluster_insp_Astro_filtered.rds")




pdata <- subseu@meta.data[, c("samplename", "seurat_clusters"), drop = FALSE] %>%
			cbind(., auc[colnames(subseu), "Human rAstro", drop = FALSE]) %>%
			setNames(., c("donors", "seurat_clusters", "auc")) %>%
			cbind(., subseu$umap@cell.embeddings) %>%
			cbind(., t(subseu$RNA@data[c("CHI3L1", "OSMR"), ])) %>%
			as.data.frame(., check.names = FALSE)
library(ggrastr)

p1 <- ggplot(pdata, aes_string(x = "UMAP_1", y = "UMAP_2", color = "donors")) +
				rasterise(geom_point(size = 0.5, shape = 16), dpi = 300, scale = 1) +
	            coord_equal(ratio = 1) + 
	            theme_classic() + 
	            theme(legend.position = "bottom",
	                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
	                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
	                    plot.title = element_blank())
exp_list <- lapply(c("auc", "CHI3L1", "OSMR"), function(gg) {
	max_exp <- max(pdata[, gg]) * 0.95
	pdata[, gg] <- MinMax(pdata[, gg], min = 0, max = max_exp)

	p <- ggplot(pdata, aes_string(x = "UMAP_1", y = "UMAP_2", color = gg)) +
				rasterise(geom_point(size = 0.5, shape = 16), dpi = 300, scale = 1) +
	            coord_equal(ratio = 1) + 
	            theme_classic() + 
	            scale_color_gradientn(colors = c("#BEBEBE", "#C4AAAA", "#E35050", "#FF0000"), limits = c(0, max_exp)) + 
	            theme(legend.position = "bottom",
	                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
	                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
	                    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
	return(p)
	})


pdf(paste0(outputdir, "SF.Multiome.rAstro.pdf"), width = 4 * 4, height = 5.4)
plot <- patchwork::wrap_plots(plotlist = c(list(p1), exp_list), nrow = 1, ncol = 4) #& theme(legend.position = "bottom")
print(plot)
dev.off()



## Plot the cell type proportion across donors
meta <- readRDS(file = "../data/PFC.filtered.meta.05082022.rds") %>%
			filter(mres == "Astro" & species == "Human")
size_bat1 <- meta %>%
			mutate(fig1cluster = as.character(fig1cluster)) %>%
			group_by(samplename, fig1cluster) %>%
			summarize(ncells = n()) %>%
			ungroup() %>%
			group_by(samplename) %>%
			mutate(ratio = ncells * 100/sum(ncells), total_cells = sum(ncells)) %>%
			ungroup() %>%
			filter(fig1cluster == "Astro AQP4 OSMR") %>%
			select(samplename, ratio, total_cells) %>%
			mutate(batch = "RNA")
size_bat2 <- subseu@meta.data %>%
			group_by(samplename, seurat_clusters, .drop = FALSE) %>%
			summarize(ncells = n()) %>%
			ungroup() %>%
			group_by(samplename) %>%
			mutate(ratio = ncells * 100/sum(ncells), total_cells = sum(ncells)) %>%
			ungroup() %>%
			filter(seurat_clusters == "6") %>%
			select(samplename, ratio, total_cells) %>%
			mutate(batch = "Multiome")
sp_ord <- c(size_bat1$samplename, size_bat2$samplename)
size_data <- rbind(size_bat1, size_bat2) %>%
			tidyr::gather(., "type", "value", c("ratio", "total_cells")) %>%
			mutate(samplename = factor(as.character(samplename), levels = sp_ord))


p <- ggplot(size_data, aes_string(x = "samplename", y = "value", fill = "batch")) +
			geom_bar(stat = "identity", size = 0.2) +
			theme_classic() +
			RotatedAxis() +
			facet_wrap(vars(type), nrow = 2, ncol = 1, scales = "free_y") +
			theme(strip.background = element_blank(),
				axis.line = element_line(size = 0.2),
				axis.ticks = element_line(size = 0.2),
				axis.text = element_text(size = rel(0.8)))
pdf(paste0(outputdir, "SF.rAstro.ratio.twobatches.pdf"), width = 4, height = 4)
print(p)
dev.off()


size_data2 <- rbind(size_bat1, size_bat2) %>%
			ungroup() %>%
			mutate(age = c(64, 36, 19, 50, 45, 60, 43, 68, 51))
p2 <- ggplot(size_data2, aes_string(x = "age", y = "ratio", color = "batch")) +
			geom_point(size = 2, shape = 16) +
			theme_classic() +
			RotatedAxis() +
			theme(axis.line = element_line(size = 0.2),
				axis.ticks = element_line(size = 0.2),
				axis.text = element_text(size = rel(0.8)))
pdf(paste0(outputdir, "SF.rAstro.ratio.across-age.pdf"), width = 4, height = 3.5)
print(p2)
dev.off()
























