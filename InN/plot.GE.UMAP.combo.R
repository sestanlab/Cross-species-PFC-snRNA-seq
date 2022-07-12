## Plot the integrated UMAP of MGE/CGE cells
args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")


## Subset to all a specific GE
ge <- args[1]
cls <- switch(ge, 
			CGE = c("ADARB2", "LAMP5_RELN", "VIP"),
			MGE = c("PVALB", "PVALB_CHC", "SST", "SST_NPY", "TH", "LAMP5_LHX6"))


## Run the UMAP for the subset of cells
pdataFile <- paste0(inputdir, "GE.UMAP", ge, ".pdata.rds")
if (!file.exists(pdataFile)){
	pfc <- readRDS(file = paste0(dataDir, "Final_PFC_HPRC.InN.07102021.rds"))
	subseu <- pfc %>%
				subset(mres %in% cls)
	subseu <- RunUMAP(subseu, dims = 1:30, seed.use = 10, min.dist = 0.4, n.neighbors = 35)


	## Plot the UMAP with cell type labeled.
	p_anno <- DimFig(subseu, file = "asdf", group.by = "fig1cluster", output.ggplot = TRUE)[[2]]
	pdf(paste0(outputdir, "InN.UMAPcombo.", ge, ".sup1.pdf"), 7, 7)
	print(p_anno)
	dev.off()



	## Do all the related plots
	plot_data <- subseu@meta.data[, c("fig1cluster", "mres", "species")] %>%
					cbind(., subseu$umap@cell.embeddings)
	layer_infor <- readRDS(file = paste0(inputdir, "PFC_transed_layers_by_HMTG.meta.rds"))
	plot_data$layer <- layer_infor[rownames(plot_data), "predicted.id"]
	saveRDS(plot_data, file = pdataFile)	
}



plot_data <- readRDS(file = pdataFile)



## randomalize the order
set.seed(42)
plot_data <- plot_data[sample(1:nrow(plot_data)), ]


cls_anno <- readRDS(paste0(dataDir, "cluster.color.rds"))
cls_cols <- setNames(cls_anno$color, rownames(cls_anno))
cls_cols <- cls_cols[names(cls_cols) %in% as.character(plot_data$fig1cluster)]


sp_cols <- setNames(c("#FF420E","#4CB5F5","#89DA59","#FFBB00"), c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
gp_cols <- gg_color_hue(length(cls)) %>% setNames(., cls)
layer_cols <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f") %>% setNames(., paste0("L", 1:6))



## Plot the subregion
cols <- list(fig1cluster = cls_cols,
			mres = gp_cols,
			species = sp_cols, 
			layer = layer_cols)
plist <- lapply(names(cols), function(x) {
	p <- ggplot(plot_data, aes_string(x = "UMAP_1", y = "UMAP_2", color = x)) +
			geom_point(size = 0.01) + 
			scale_color_manual(values = cols[[x]])+
			theme_classic() + 
			theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.position = "none", plot.margin = unit(c(0, -0.1, 0, -0.1), "in"))
	p
	})



plot.scale <- 1
jpeg(paste0(outputdir, "InN.UMAPcombo.", ge, ".jpeg"), width = plot.scale * 4 * 6, height = plot.scale * 6, units = "in", res = 300)
plot_grid(plotlist = plist, nrow = 1, ncol = 4) %>% print()
dev.off()






## Plot the legend 
lp_data <- data.frame(cluster = lapply(cols, names) %>% unlist(use.names = FALSE), 
					color = unlist(cols, use.names = FALSE),
					stringsAsFactors = FALSE) %>%
				mutate(value = 1) %>%
				mutate(cluster = factor(cluster, levels = lapply(cols, names) %>% unlist(use.names = FALSE)))
lp <- ggplot(lp_data, aes(y = cluster,  x = 1, color = cluster))+
		geom_point(size = 3) + 
		scale_color_manual(values = setNames(lp_data$color, lp_data$cluster)) +
		theme_classic() +
		theme(legend.position = "right")
pdf(paste0(outputdir, "InN.UMAPcombo.", ge, ".legend.pdf"), width = 4, height = 6)
print(lp)
dev.off()










