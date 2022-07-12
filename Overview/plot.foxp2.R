source("../scripts/pfc.fun.R")
source("../MF4_FOXP2/foxp2.fun.R") 



all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
pdata <- lapply(all_sps, function(sp) {
	seu <- readRDS(file = paste0(inputdir, "Cluster_insp_", sp, ".rds"))

	DefaultAssay(seu) <- "RNA"
	genes <- c("FOXP2", "TH")
	mapgenes <- unlist(mclapply(genes, function(i) ifelse(length(grep(paste0("-",i,"$"), rownames(seu), ignore.case= TRUE))==1, rownames(seu)[grep(paste0("-",i,"$"),rownames(seu), ignore.case= TRUE)], "empty")))
	mapgenes <- mapgenes[mapgenes != "empty"]
	print(mapgenes)

	exp_data <- t(as.matrix(seu$RNA@data[mapgenes, , drop = FALSE]))
	colnames(exp_data) <- extract_field(colnames(exp_data), -1, "-")

	data <- data.frame(cell = rownames(seu@meta.data),
					cluster = as.character(seu@meta.data$fig1cluster),
					mres = as.character(seu@meta.data$mres),
					species = sp,
					stringsAsFactors = FALSE) %>%
				cbind(., seu$umap@cell.embeddings)
	print(paste0("Finish species: ", sp))
	return(list(meta = data, exp = exp_data))
	}) %>%
		setNames(., all_sps)

plot_data <- lapply(pdata, function(x) {
	cbind(x$meta, x$exp)
	}) %>%
		do.call(rbind, .)



cls_ord <- c("L2-3 IT","L3-5 IT","L5 ET","L5-6 NP","L6 CT","L6 IT","L6B","LAMP5","ADARB2 KCNG1","VIP","SST","PVALB","Astro","OPC","Oligo","Micro","Immune","Vas")
plot_data$mres2 <- gsub("L3-5 IT-[0-9]", "L3-5 IT", plot_data$mres) %>%
					gsub("L6 IT-[0-9]", "L6 IT", .) %>%
					gsub(" HGF| NPY", "", .) %>%
					gsub(" ChC", "", .) %>%
					gsub(" LHX6| RELN", "", .) %>%
					gsub("Endo|SMC|VLMC", "Vas", .) %>%
					gsub("^PC$", "Vas", .) %>%
					gsub("^RB$", "Vas", .)
plot_data$mres2 <- factor(plot_data$mres2, levels = cls_ord)
plot_data$species <- factor(plot_data$species, levels = all_sps)	



## Plot Type
me <- max(plot_data$FOXP2)
p <- ggplot() + 
		geom_violin(data = plot_data, aes_string(x = "mres2", y = "FOXP2", fill = "species"), scale = "width", size = 0.1, adjust = 2,trim =TRUE) + 
		theme_cowplot() + 
					RotatedAxis() + 
					scale_fill_manual(values = c("#FF420E", "#4CB5F5", "#89DA59", "#FFBB00", "#AE017E") %>% setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset", "Mouse"))) + 
					facet_wrap(vars(species), nrow = 4, ncol = 1, strip.position = 'left') +
					scale_y_continuous(breaks = seq(0, me, by = round(me/2, digits = 1))) + 
					theme(axis.text.x=element_text(size = 8), strip.background = element_blank(), strip.text = element_text(size = 10, angle = 90), panel.spacing = unit(0.05, "in"), strip.placement = 'outside', axis.line = element_line(size = 0.25), axis.ticks = element_line(size = 0.25), legend.position = "bottom", legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size = rel(0.5)), legend.title = element_blank())
pdf(paste0(outputdir, "FOXP2_expression_nonLiftover.pdf"), width =6, height = 4)
print(p)
dev.off() 


saveRDS(plot_data, file = paste0(inputdir, "FOXP2_expression_nonLiftover.rds"))










