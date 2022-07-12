source("../scripts/pfc.fun.R")


marlist <- readRDS(file = paste0(inputdir, "LAMP5.speccls.putative.markers.rds"))


allmeta <- readRDS(file = paste0(dataDir, "PFC.filtered.meta.05082022.rds"))
all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
seu_list <- lapply(all_sps, function(sp) {
	seu <- readRDS(file = paste0(inputdir, "CGE_insp_", sp, ".rds"))
	seu@meta.data$fig1cluster <- as.character(allmeta[colnames(seu), "fig1cluster"])
	seu@meta.data$mres <- as.character(allmeta[colnames(seu), "mres"])
	seu <- subset(seu, mres %in% c("LAMP5 RELN")) %>%
			RunUMAP(., dims = 1:30, metric = "correlation", umap.method = "umap-learn")
	seu
	}) %>%
		setNames(., all_sps)



## Add AUC scores to the data
seu_nlist <- lapply(all_sps, function(sp) {
	print(paste0("Workig on : ", sp))
	obj <- seu_list[[sp]]
	aucsub <- GetModuleScore(assay.data = obj$RNA@data, features = marlist, method = "aucell", input_dir = inputdir, file_name = paste0("LAMP5MarkerAUC_", sp), output_dir = outputdir)$auc
	obj[["AUC"]] <- CreateAssayObject(data = t(as.matrix(aucsub)), min.cells = 0, min.features = 0)
	return(obj)
	}) %>%
		setNames(., all_sps)


## Visualize the AUC scores
pdata <- lapply(all_sps, function(sp){
	obj <- seu_nlist[[sp]]
	data <- as.matrix(t(obj[["AUC"]]@data)) %>%
				cbind(., obj@meta.data[, c("species", "fig1cluster", "samplename")]) %>%
				cbind(., obj$umap@cell.embeddings) %>%
				tidyr::gather(., "type", "auc", rownames(obj[["AUC"]])) %>%
				mutate(samplename = paste0("rep_", as.numeric(as.factor(samplename))))
	return(data)
	}) %>%
		do.call(rbind, .)
set.seed(42)
pdata <- pdata[sample(1:nrow(pdata)), ]




all_cls <- levels(as.factor(pdata$fig1cluster))
cls_cols <- setNames(gg_color_hue(length(all_cls)), all_cls)
rep_cols <- setNames(c("#1b9e77", "#e41a1c", "#6a3d9a", "#a65628"), paste0("rep_", 1:4))
library(ggrastr)
set.seed(42)
anno_list <- lapply(all_sps, function(sp){
		sp_data <- filter(pdata, species == sp)
		p1 <- ggplot(sp_data, aes_string(x = "UMAP_1", y = "UMAP_2", color = "fig1cluster")) +
				rasterise(geom_point(size = 0.5, shape = 16), dpi = 300, scale = 1) +
	            theme_classic() + 
	            scale_color_manual(values = cls_cols) + 
	            theme(legend.position = "none",
	                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
	                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
	                    plot.title = element_blank())
	    p2 <- ggplot(sp_data, aes_string(x = "UMAP_1", y = "UMAP_2", color = "samplename")) +
				rasterise(geom_point(size = 0.5, shape = 16), dpi = 300, scale = 1) +
	            theme_classic() + 
	            scale_color_manual(values = rep_cols) + 
	            theme(legend.position = "none",
	                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
	                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
	                    plot.title = element_blank())
	    list(p1, p2)
		}) %>%
		do.call(c, .)





all_types <- "Marmoset SYT10" 
type_list <- list()
for (cur_type in all_types){
	max_auc <- pdata %>%
			filter(type == cur_type) %>% .$auc %>%
			max()
	max_auc <- max_auc * 0.95
	plist <- lapply(all_sps, function(sp){
		sp_data <- filter(pdata, species == sp & type == cur_type)
		sp_data$auc <- MinMax(sp_data$auc, min = 0, max = max_auc)
		p <- ggplot(sp_data, aes_string(x = "UMAP_1", y = "UMAP_2", color = "auc")) +
				rasterise(geom_point(size = 0.6, shape = 16), dpi = 300, scale = 1) +
	            theme_classic() + 
	            scale_color_gradientn(colors = c("#BEBEBE", "#C4AAAA", "#E35050", "#FF0000"), limits = c(0, max_auc)) + 
	            theme(legend.position = "bottom",
	                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
	                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
	                    plot.title = element_text(size = 12, hjust = 0.5, face = "bold")) ##+ 
		})
	type_list[[cur_type]] <- plist
}

total <- length(all_types) * 4

final_list <- c(anno_list[1:2], do.call(c, type_list)[seq(1, total, 4)], 
				anno_list[3:4], do.call(c, type_list)[seq(2, total, 4)], 
				anno_list[5:6], do.call(c, type_list)[seq(3, total, 4)], 
				anno_list[7:8], do.call(c, type_list)[seq(4, total, 4)])
pdf(paste0(outputdir, "SF.validate.spec-clusters.SYT10.AUC.cbn.pdf"), width = 3 * (2 + length(all_types)), height = 3 * 4)
plot <- patchwork::wrap_plots(plotlist = final_list, nrow = 4, ncol = 2 + length(all_types), guides = "collect") & theme(legend.position = "none")
print(plot)
dev.off()



