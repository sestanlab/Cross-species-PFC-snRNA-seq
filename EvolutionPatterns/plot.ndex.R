## Plot the number of DEX genes for each clusters 
source("../scripts/pfc.fun.R")
library(ggpubr)


## Calculate the median nGenes in each cluster in each species
expr_mat <- readRDS(file = paste0(dataDir, "Expr_avg_ratio_by_species.rds"))$fig1cluster$ratio

all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
homo_cls <- readRDS(file = paste0(dataDir, "homolog.cluster.rds"))

ng_mat <- lapply(all_sps, function(sp) {
	sp_cls <- grep(sp, colnames(expr_mat), value = TRUE)
	cls <- extract_field(sp_cls, 2, "|") %>% setNames(., NULL)
	out <- colSums(expr_mat[, sp_cls] >= 0.3) %>% setNames(., cls)
	return(out[homo_cls])
	}) %>%
		setNames(., all_sps) %>%
		as.data.frame() %>%
		as.matrix()




## Get the nDEX genes
all.dex <- readRDS(file = paste0(dataDir, "DEX.hprc.", "fig1cluster", ".pre.rds"))
sp_pairs <- c("Human|Chimpanzee", "Human|Rhesus", "Chimpanzee|Rhesus", "Rhesus|Marmoset", "Human|Marmoset", "Chimpanzee|Marmoset")
dex_abs <- lapply(sp_pairs, function(pp) {
	subdata <- filter(all.dex, speciespair == gsub("\\|", "_", pp))

	pos_idx <- subdata$pct.1 >= 0.3 & subdata$ratio_fc >= 1.25 & subdata$avg_logFC >= log(2.5) & subdata$p_val_adj <= 0.01
	neg_idx <- subdata$pct.2 >= 0.3 & subdata$ratio_fc <= 0.8 & subdata$avg_logFC <= log(0.4) & subdata$p_val_adj <= 0.01

	fdata <- subdata[pos_idx | neg_idx, ]
	mars <- split(fdata$gene, fdata$cluster) %>%
				.[homo_cls]
	ndex <- sapply(mars, length) %>%
				setNames(., homo_cls)
	return(ndex)
	}) %>%
			setNames(., sp_pairs) %>%
			as.data.frame(., check.names = FALSE)
dex_abs$group <- 'NNC'
dex_abs$group[grep("^L[0-9]", rownames(dex_abs))] <- "ExN"
dex_abs$group[grep("^InN", rownames(dex_abs))] <- "InN"
dex_abs$group[grep("Astro|OPC|COP|Oligo|Micro", rownames(dex_abs))] <- "Glia"






## Get the relative number of DEX (DEX percentage = n DEX normalized by the # expressed genes)
dex_rel <- dex_abs
for (pp in sp_pairs){
	sps <- strsplit(pp, "|", fixed = TRUE)[[1]]
	cls_umi <- ng_mat[, sps] %>%
				rowMeans()
	dex_rel[, pp] <- dex_rel[, pp]/cls_umi[rownames(dex_abs)]
}



## Incorporate the 1-cor
div_mat <- readRDS(file = paste0(inputdir, "DIV_1minusCor_raw.rds"))


data_list <- lapply(colnames(div_mat), function(pair) {
	plot_data <- data.frame(cluster = rownames(div_mat), 
							div = div_mat[, pair], 
							ndex = dex_abs[rownames(div_mat), pair], 
							rdex = dex_rel[rownames(div_mat), pair], 
							group = dex_abs[rownames(div_mat), "group"],
							stringsAsFactors = FALSE)
	return(plot_data)
	}) %>% 
		setNames(., colnames(div_mat))


gp_cols <- setNames(c("#C51B7D", "#4D9221", "#08519c", "#999999"), c("ExN", "InN", "Glia", "NNC"))
plist <- lapply(c("ndex", "rdex"), function(type) {
	min_x <- min(div_mat)
	max_x <- max(div_mat)
	max_y <- lapply(data_list, function(mm) max(mm[, type])) %>% unlist() %>% max()
	subp <- lapply(names(data_list), function(x) {
		p <- ggplot(data_list[[x]], aes_string(x = "div", y = type)) + 
				geom_point(aes_string(color = "group", fill = "group"), size = 2.5, alpha = 0.7, shape = 21) +
				theme_bw() +
				scale_color_manual(values = gp_cols) +
				scale_fill_manual(values = gp_cols) +
				scale_x_continuous(limits = c(min_x, max_x)) +
				scale_y_continuous(limits = c(0, max_y)) +
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), legend.position = "right")
		p
		})
	subp
	}) %>%
		do.call(c, .)

plist[[1]] <- plist[[1]] + theme(axis.text.y = element_text(size = rel(0.9)))
plist[[7]] <- plist[[7]] + theme(axis.text.y = element_text(size = rel(0.9)))
for (i in 7:12){
	plist[[i]] <- plist[[i]] + theme(axis.text.x = element_text(size = rel(0.9)))
}


pdf(paste0(outputdir, "NDEX_vs_div.dot.pdf"), width = 19, height = 6, useDingbats = FALSE)
patchwork::wrap_plots(plist, nrow = 2, ncol = 6, guides = "collect") %>% print()
dev.off()

















