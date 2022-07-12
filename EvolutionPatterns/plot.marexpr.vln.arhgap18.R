source("../scripts/pfc.fun.R")


allmeta <- readRDS(file = paste0(dataDir, "PFC.filtered.meta.05082022.rds"))
all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
seu_list <- lapply(all_sps, function(sp) {
	seu <- readRDS(file = paste0(inputdir, "L23IT_insp_", sp, ".rds"))
	seu@meta.data$fig1cluster <- as.character(allmeta[colnames(seu), "fig1cluster"])
	seu@meta.data$mres <- as.character(allmeta[colnames(seu), "mres"])
	seu
	}) %>%
		setNames(., all_sps)


genes <- c("CUX2", "PRLR")
## Visualize the AUC scores
pdata <- lapply(all_sps, function(sp){
	obj <- seu_list[[sp]]
	data <- as.matrix(t(obj[["RNA"]]@data[genes, ])) %>%
				cbind(., obj@meta.data[, c("species", "fig1cluster")]) %>%
				cbind(., obj$umap@cell.embeddings) %>%
				tidyr::gather(., "gene", "expr", genes) %>%
	return(data)
	}) %>%
		do.call(rbind, .) %>%
		mutate(species = factor(species, levels = all_sps)) %>%
		mutate(gene = factor(gene, levels = genes))


p <- ggplot() + 
		geom_violin(data = pdata, aes_string(x = "fig1cluster", y = "expr", fill = "species"), scale = "width", size = 0.1, adjust = 2,trim =TRUE) + 
		theme_cowplot() + 
		RotatedAxis() + 
		scale_fill_manual(values = c("#FF420E", "#4CB5F5", "#89DA59", "#FFBB00", "#AE017E") %>% setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset", "Mouse"))) + 
		facet_wrap(vars(gene), nrow = length(genes), ncol = 1, strip.position = 'left', scales = "free_y") +
		theme(axis.text.x=element_text(size = 8), strip.background = element_blank(), strip.text = element_text(size = 10, angle = 90), panel.spacing = unit(0.05, "in"), strip.placement = 'outside', axis.line = element_line(size = 0.25), axis.ticks = element_line(size = 0.25), legend.position = "bottom", legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size = rel(0.5)), legend.title = element_blank())

pdf(paste0(outputdir, "Violin_RepresentativeMarkers_L2-3 IT ARHGAP18.pdf"), width = 8, height = 4)
print(p)
dev.off()





