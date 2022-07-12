source("../scripts/pfc.fun.R") 
source("../scripts/exp.fun.R")


sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
seulist <- lapply(sps, function(sp) readRDS(paste0("../recluster_v2/load_files/Cls_bysp_detail_Immune_", sp, ".rds"))) %>%
			setNames(., sps)


## Load Meta file to update cluster identity
meta <- readRDS(file = paste0("../data/", "PFC.filtered.meta.05082022.rds"))
seulist <- lapply(names(seulist), function(sp) {
	xx <- seulist[[sp]]
	xx[["fig1cluster"]] <- as.character(meta[colnames(xx), "fig1cluster"])
	xx[["species"]] <- sp
	DefaultAssay(xx) <- "RNA"
	xx
	}) %>%
		setNames(., names(seulist))


seu <- Reduce(merge, seulist)
seu[["cluster"]] <- gsub("^B .*", "B", seu$fig1cluster) %>%
					gsub("^T .*", "T", .) %>%
					gsub("Macro.*", "Macro", .) %>%
					gsub("Micro P2RY12 APBB1IP", "Micro", .) %>%
					gsub("Micro P2RY12 CCL3", "huMicro", .) %>%
					gsub("Micro P2RY12 GLDN", "hoMicro", .) %>%
					gsub("Myeloid LSP1 LYZ", "Myeloid", .)

seu@meta.data$spcls <- paste0(seu@meta.data$species, "|", as.character(seu@meta.data$cluster))



## Calculate AVG & Exp ratio
Idents(seu) <- "spcls"
avgs <- log(AverageExpression(object = seu, assay = "RNA", slot = "data")$RNA + 1) %>% 
				as.matrix()

all_cls <- levels(as.factor(seu@meta.data[, "spcls"]))
ratios <- lapply(all_cls, function(x) {
			message(paste0("Calculate expression ratio for cluster: ", x))	
			sub_cdt <- seu$RNA@data[, seu@meta.data[, "spcls"] == x, drop = FALSE]
			er <- Matrix::rowMeans(sub_cdt != 0)
			er
			}) %>% setNames(., all_cls) %>%
				as.data.frame(., check.names = FALSE) %>% as.matrix()



## Plot data
genes <- c("P2RY12", "FOS", "CD83", "CCL3", "CCL4", "CH25H", "CD86", "BTG2", "CDKN1A", "GLDN", "MYO1E", "PADI2", "PPARG", "F13A1", "LYZ", "CD247", "EBF1")
sp_cols <- setNames(c("#FF420E","#4CB5F5","#89DA59","#FFBB00"), c("Human", "Chimpanzee", "Rhesus", "Marmoset"))

avg_use <- avgs[genes, ] %>%
			as.matrix() %>%
			t() %>% scale() %>% t() %>%
			MinMax(., min = -2.5, max = 2.5)



plot_data <- CombnMelt(avg = avg_use[genes, ], ratio = ratios[genes, ], stroke = NULL) %>%
				mutate(exp = as.numeric(cut(avg.exp.scaled, 20))) %>%
				mutate(color = sp_cols[species]) %>%
				mutate(spcls = paste0(species, "|", cluster))
plot_data$color <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("lightgrey", color))(20)[value])
        }, color = sp_cols[as.character(plot_data$species)], value = plot_data$exp)
plot_data$pct.exp <- plot_data$pct.exp * 100



cls_order <- data.frame(spcls = colnames(avg_use), stringsAsFactors = FALSE) %>%
				mutate(species = extract_field(spcls, 1, "|"), cluster = extract_field(spcls, 2, "|")) %>%
				mutate(species = factor(species, levels = c("Human", "Chimpanzee", "Rhesus", "Marmoset")), 
						cluster = factor(cluster, levels = c("Micro", "huMicro", "hoMicro", "Macro", "Myeloid", "T", "B"))) %>%
				arrange(species, cluster) %>%
				.$spcls



p <- ggplot(plot_data, aes(x = spcls, y = features.plot, color = color, size = pct.exp)) +
			geom_point(shape = 16) +
			scale_radius(range = c(0, 5.5), limits = c(0, 100)) +
			scale_color_identity() +
			theme_classic() +
			scale_x_discrete(limits = cls_order, labels = extract_field(cls_order, 2, "|")) +
			scale_y_discrete(limits = rev(genes)) +
			RotatedAxis() + 
			theme(legend.position = "bottom", axis.line = element_blank(), axis.text.x = element_text(size = rel(0.8)), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = rel(0.8)), panel.grid.major = element_line(size = 0.15, color = "#4D4D4D"))


pdf(paste0(outputdir, "MF3.Immune.micro.spec.genes.v2.pdf"), width = 6, height = 3.7, useDingbats = FALSE)
print(p)
dev.off()



