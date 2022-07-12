## Plot the expression of SST across layers
source("../scripts/pfc.fun.R")



###---------------------------------------------------------------------------------
## Visualize in primates (HPRC)
pfc <- readRDS(paste0(dataDir, file = "Final_PFC_HPRC.InN.07102021.rds"))
pfc <- subset(pfc, fig1cluster != "InN LAMP5 SYT10")



## Clean the laminar information [Remove some noisy layer information]
layer_infor <- readRDS(file = paste0(inputdir, "PFC_transed_layers_by_HMTG.meta.rds"))
pfc@meta.data$layer <- layer_infor[colnames(pfc), "predicted.id"]


all_cls <- levels(as.factor(as.character(pfc@meta.data$fig1cluster)))
all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
cells <- lapply(all_sps, function(sp) {
	sp_seu <- subset(pfc, species == sp)
	sp_cells <- lapply(all_cls, function(cls) {
		incls <- sp_seu@meta.data$fig1cluster == cls
		sl <- table(sp_seu@meta.data$layer[incls])/sum(incls)
		sl <- names(sl)[sl >= 0.025]
		cs <- colnames(sp_seu)[incls & sp_seu@meta.data$layer %in% sl]
		cs
		}) %>%
			unlist()
	sp_cells
	}) %>%
		unlist()



## Get the plot data
plot_data <- data.frame(cell = colnames(pfc),
						layer = pfc@meta.data$layer, 
						species = pfc@meta.data$species,
						exp = pfc$RNA@data["SST", ],
						stringsAsFactors = FALSE) %>%
					mutate(species = factor(as.character(species), levels = c("Human", "Chimpanzee", "Rhesus", "Marmoset")))



## Violin plot showing the laminar expression of SST across species.
p <- ggplot(plot_data, aes_string(x = "layer", y = "exp", fill = "species")) + 
			geom_violin(color = "black", scale = "width", size = 0.1, adjust = 2, trim =TRUE) + 
			theme_cowplot() + 
			RotatedAxis() + 
			scale_fill_manual(values = c("#FF420E", "#4CB5F5", "#89DA59", "#FFBB00") %>% setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset"))) + 
			theme(axis.text.x=element_text(size = rel(1)), strip.background = element_blank(), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), legend.position = "bottom")

pdf(paste0(outputdir, "SST_Laminar_exp.species_v2.pdf"), width = 7, height = 3)
print(p)
dev.off()



## Calculate the significance of SST laminar enrichment in human vs other 3 species
lapply(paste0("L", 1:6), function(x) {
	sapply(c("Chimpanzee", "Rhesus", "Marmoset"), function(sp) {
		sub_data <- filter(plot_data, layer == x & species %in% c("Human", sp)) %>%
						mutate(species = factor(as.character(species), levels = c("Human", sp)))
        means <- aggregate(exp ~ species, sub_data, mean) %>%
                    column_to_rownames("species") %>% as.matrix() %>% .[, 1]
        fc <- (means["Human"] + 0.01)/(means[sp] + 0.01)
		p_val <- wilcox.test(exp ~ species, sub_data, alternative = "greater")$p.value
        if (fc <= 2){
            p_val <- 1
        }
        return(p_val)
		}) %>% max()
	})


##------------------------------------------------------------
## Plot laminar scores
pfc <- readRDS(paste0(dataDir, file = "Final_PFC_HPRC.InN.07102021.rds"))
pfc <- subset(pfc, fig1cluster != "InN LAMP5 SYT10")


aa <- pfc@meta.data
aa$fig1cluster <- as.character(aa$fig1cluster)
aa$layer <- factor(as.character(aa$layer), levels = paste0("L", 1:6))
layer_data <- aa %>%
				group_by(fig1cluster, layer, .drop=FALSE) %>%
				summarize(ncells = n()) %>%
				mutate(ratio = ncells/sum(ncells)) %>%
				ungroup() %>%
				mutate(fig1cluster = factor(fig1cluster, levels = rev(cls_order)))


p.anno <- ggplot(layer_data) + 
				geom_tile(aes(x = layer, y = fig1cluster, fill = ratio), color = "lightgrey", height = 1, width = 1, size = 0.1) + 
				theme_classic() +
				scale_fill_gradientn(colors = c("#FFFFFF", "#b2abd2","#542788")) +
				RotatedAxis() + 
	            theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 8), axis.line = element_blank(), axis.title = element_blank(), axis.ticks.x = element_blank(), legend.position = "left")
pdf(paste0(outputdir, "MF6_InN_SST_layer1_laminarInfor.pdf"), width = 15, height = 7, useDingbats = FALSE)
patchwork::wrap_plots(p.anno, p1, p2, nrow = 1, ncol = 3)
dev.off()







