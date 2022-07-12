source("../scripts/pfc.fun.R")


marlist <- readRDS(file = paste0("../SpecclsMars/load_files/", "LAMP5.speccls.putative.markers.rds"))


seu_list <- list(Motor_maromoset = readRDS(file = paste0(inputdir, "Motor_marmoset_LAMP5_RELN_inte_sub.rds")),
				Motor_human = readRDS(file = paste0(inputdir, "Motor_Human_LAMP5_RELN_inte_sub_fil.rds")))
aa <- readRDS(file = paste0("../MF6_InN/load_files/LAMP5_RELN_sepcls_all.rds"))[["Marmoset"]]
seu_list[["dlPFC_marmoset"]] <- aa


## Add AUC scores to the data
seu_nlist <- lapply(names(seu_list), function(sp) {
	print(paste0("Workig on : ", sp))
	obj <- seu_list[[sp]]
	aucsub <- GetModuleScore(assay.data = obj$RNA@data, features = marlist, method = "aucell", input_dir = inputdir, file_name = paste0("LAMP5MarkerAUC_", sp), output_dir = outputdir)$auc
	obj[["AUC"]] <- CreateAssayObject(data = t(as.matrix(aucsub)), min.cells = 0, min.features = 0)
	return(obj)
	}) %>%
		setNames(., names(seu_list))


##--------------------------------------------------------------------------------
## Visualize the AUC scores
pdata <- lapply(names(seu_nlist), function(sp){
	obj <- seu_nlist[[sp]]
	obj$cluster <- ifelse_check(sp == "dlPFC_marmoset", as.character(obj$fig1cluster), as.character(obj$cross_species_cluster_label))
	data <- as.matrix(t(obj[["AUC"]]@data)) %>%
				cbind(., obj@meta.data[, c("cluster"), drop = FALSE]) %>%
				cbind(., obj$umap@cell.embeddings) %>%
				tidyr::gather(., "type", "auc", rownames(obj[["AUC"]])) %>%
				mutate(species = sp)
	return(data)
	}) %>%
		do.call(rbind, .)
set.seed(42)
pdata <- pdata[sample(1:nrow(pdata)), ]



library(ggrastr)
set.seed(42)
anno_list <- lapply(names(seu_nlist), function(sp){
		sp_data <- filter(pdata, species == sp)
		all_cls <- levels(as.factor(sp_data$cluster))
		cls_cols <- setNames(gg_color_hue(length(all_cls)), all_cls)
		p1 <- ggplot(sp_data, aes_string(x = "UMAP_1", y = "UMAP_2", color = "cluster")) +
				rasterise(geom_point(size = 0.5, shape = 16), dpi = 300, scale = 1) +
	            theme_classic() + 
	            scale_color_manual(values = cls_cols) + 
	            theme(legend.position = "none",
	                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
	                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
	                    plot.title = element_blank())
	    
	    p2 <- lapply(names(marlist), function(tp) {
	    	max_auc <- pdata %>%
						filter(type == tp) %>% .$auc %>% max()
			sp_data <- filter(pdata, species == sp & type == tp)
			sp_data$auc <- MinMax(sp_data$auc, min = 0, max = max_auc)
			p <- ggplot(sp_data, aes_string(x = "UMAP_1", y = "UMAP_2", color = "auc")) +
				rasterise(geom_point(size = 0.6, shape = 16), dpi = 300, scale = 1) +
	            theme_classic() + 
	            scale_color_gradientn(colors = c("#BEBEBE", "#C4AAAA", "#E35050", "#FF0000"), limits = c(0, max_auc)) + 
	            theme(legend.position = "bottom",
	                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
	                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
	                    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
	        p
	    	})

	    return(c(list(p1), p2))
		}) %>%
		do.call(c, .)
pdf(paste0(outputdir, "SF.SYT10.spec-scores.motor.dlPFC.pdf"), width = 2.5 * 3, height = 2.5 * 3)
plot <- patchwork::wrap_plots(plotlist = anno_list, nrow = 3, ncol = 3, guides = "collect") & theme(legend.position = "none")
print(plot)
dev.off()












