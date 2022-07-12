source("../scripts/pfc.fun.R")
library(ggpubr)


resolution <- "fig1cluster"
all_species <- c('Human', 'Chimpanzee', 'Rhesus', 'Marmoset')


## Similarity to CGE, MGE
pfc <- readRDS(paste0(dataDir, file = "Final_PFC_HPRC.InN.07102021.rds"))
sel.cls <- c("InN LAMP5 LHX6 PROX1", "InN LAMP5 LHX6 TAC1")
seu <- subset(pfc, fig1cluster %in% sel.cls)



## subtype markers
all.sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
mar_df <- lapply(all.sps, function(sp) {
	sub.seu <- subset(seu, species == sp)
	Idents(sub.seu) <- "fig1cluster"
	sp_df <- lapply(sel.cls, function(cls){
					m1 <- FindMarkers(sub.seu, ident.1 = cls,pct.1 = 0.2, only.pos = TRUE, logfc.threshold = 0.2) %>%
						rownames_to_column("gene") %>%
						mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
						filter(ratio_fc >= 1.5) %>%
						mutate(species = sp) %>%
						mutate(cluster = cls)
					m1
				}) %>%
					do.call(rbind, .)
	sp_df
	}) %>%
		setNames(., all.sps)

mars <- lapply(mar_df, function(xx) {
	split(xx$gene, xx$cluster) %>%
		.[sel.cls]
	})


## Calculate similarity
avgs <- readRDS(file = paste0(dataDir, "Expr_avg_ratio_by_species.rds"))$fig1cluster$avg

plot_data <- lapply(all.sps, function(sp) {
	g1 <- mars[[sp]][["InN LAMP5 LHX6 PROX1"]]
	g2 <- mars[[sp]][["InN LAMP5 LHX6 TAC1"]]

	avg_use <- avgs[, grep(paste0(sp, "\\|InN"), colnames(avgs), value = TRUE)]
	mge_cls <- grep("SST|PVALB|LHX6", colnames(avg_use), value = TRUE) %>%
				grep("PROX1|TAC1|SYT10", ., value = TRUE, invert = TRUE)

	p1data <- as.matrix(avg_use[g1, ]) %>%
					t() %>% scale() %>% t() %>%
					MinMax(., min = -3, max = 3) %>%
					reshape2::melt() %>%
					setNames(., c("gene", "cluster", "exp")) %>%
					mutate(gtype = "PROX1 markers")
	p2data <- as.matrix(avg_use[g2, ]) %>%
					t() %>% scale() %>% t() %>%
					MinMax(., min = -3, max = 3) %>%
					reshape2::melt() %>%
					setNames(., c("gene", "cluster", "exp")) %>%
					mutate(gtype = "TAC1 markers")
	pdata <- rbind(p1data, p2data) %>%
					mutate(GE = ifelse(as.character(cluster) %in% mge_cls, "MGE", "CGE")) %>%
					mutate(GE = factor(GE, levels = c("MGE", "CGE"))) %>%
					mutate(species = sp) %>%
					mutate(newcls = paste0(species, "|", gtype))
	return(pdata)
	}) %>%
		do.call(rbind, .)


p <- ggplot(plot_data, aes_string(x = "newcls", y = "exp", fill = "GE")) +
				geom_boxplot(outlier.shape = NA, size = 0.2, notch = TRUE, alpha = 0.7) +
                theme_cowplot() + 
                RotatedAxis() + 
                stat_compare_means(label = "p.signif", method = "wilcox.test", size = 5, symnum.args = list(cutpoints = c(0, 0.01, 0.05, 1), symbols = c("**", "*", "ns"))) +
                scale_fill_manual(values = c(MGE = "#f46d43", CGE = "#542788")) + 
                scale_x_discrete(limits = paste0(rep(all.sps, times = 2), "|", rep(c("TAC1 markers", "PROX1 markers"), each = 4))) +
                theme(axis.text.x=element_text(size = 8), axis.text.y = element_text(size = 8), axis.line = element_line(size = 0.15), axis.ticks = element_line(size = 0.15), axis.title.y = element_blank(), legend.position = "right")

pdf(paste0(outputdir, "LAMP5_LHX6_MCGE_similarity.violin.pdf"), width = 8, height = 5, useDingbats = FALSE)
print(p)
dev.off()




## Find the key markers that drive the differences
sp <- "Human"

g1 <- mars[[sp]][["InN LAMP5 LHX6 PROX1"]]
g2 <- mars[[sp]][["InN LAMP5 LHX6 TAC1"]]

avg_use <- avgs[, grep(paste0(sp, "\\|InN"), colnames(avgs), value = TRUE)]
mge_cls <- grep("SST|PVALB|LHX6", colnames(avg_use), value = TRUE) %>%
				grep("PROX1|TAC1|SYT10", ., value = TRUE, invert = TRUE)
cge_cls <- grep("SST|PVALB|LHX6", colnames(avg_use), value = TRUE, invert = TRUE) %>%
				grep("PROX1|TAC1|SYT10", ., value = TRUE, invert = TRUE)


exp_data <- as.matrix(avg_use[g1, ]) %>%
					t() %>%
					as.data.frame(., check.names = FALSE) %>%
					aggregate(., by = list(GE = ifelse(colnames(avg_use) %in% mge_cls, "MGE", "CGE")), FUN = mean) %>%
					column_to_rownames("GE") %>%
					as.matrix()

(exp_data["CGE", ]/exp_data["MGE", ]) %>% sort(., decreasing = TRUE)


exp_data <- as.matrix(avg_use[g2, ]) %>%
					t() %>%
					as.data.frame(., check.names = FALSE) %>%
					aggregate(., by = list(GE = ifelse(colnames(avg_use) %in% mge_cls, "MGE", "CGE")), FUN = mean) %>%
					column_to_rownames("GE") %>%
					as.matrix()

(exp_data["MGE", ]/exp_data["CGE", ]) %>% sort(., decreasing = TRUE)



pfc@meta.data$newcls <- pfc@meta.data$mres
idx <- which(pfc@meta.data$mres == "LAMP5_LHX6")
pfc@meta.data$newcls[idx] <- as.character(pfc@meta.data$fig1cluster[idx])
pfc@meta.data$newcls <- gsub("SST_NPY|TH", "SST", pfc@meta.data$newcls)


cls_order <- c("SST","PVALB","PVALB_CHC","InN LAMP5 LHX6 TAC1","InN LAMP5 LHX6 PROX1","LAMP5_RELN","ADARB2","VIP")
genes <- c("LAMP5", "LHX6","TAC1", "RSPO2", "CUX2", "SULF1", "SAMD5", "GHR", "GPC6", "PROX1", "THSD7B", "BMP7", "JAM2")

p <- PointPlot(pfc, features = genes, dot.scale = 5, cols = sp_cols, dot.min = 0, group.by = "newcls", split.by = "species", cluster_order = cls_order, shape = 16, species_order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), scale.by = "radius") +
				coord_flip() +
		        RotatedAxis() +
		        facet_wrap(vars(species), nrow = 1, ncol = 4) +
		        theme(axis.text.x=element_text(size = 8), strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(0.1, "in"), legend.position = "bottom")
pdf(paste0(outputdir, "LAMP5_LHX6_MCGE_exp_scaleby-radius.pdf"), width = 8, height = 7)
print(p)
dev.off()




















