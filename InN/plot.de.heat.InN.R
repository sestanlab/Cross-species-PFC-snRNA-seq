## Visualize species-exclusive markers on Heatmap
source("../scripts/pfc.fun.R")


##------------------------------------------------------------------
# load data
resolution <- "fig1cluster"
all.markers <- readRDS(file = paste0(dataDir, "Markers.hprc.", resolution, ".rds"))
all.dex <- readRDS(file = paste0(dataDir, "DEX.hprc.", resolution, ".rds")) 



# cluster order
cls_order <- readRDS(file = paste0(dataDir, "PFC.cluster.order.rds"))
homocls <- readRDS(file = paste0(dataDir, "homolog.cluster.rds"))
cls_order <- cls_order[match(homocls, cls_order)]
cls_order <- cls_order[c(1:14, 16:20, 15, 21:length(cls_order))]


ctp <- c("ExN", "InN", "NNC")[2]
homo.order <- switch(ctp, ExN = cls_order[1:36],
						InN = cls_order[37:81],
						NNC = cls_order[82:103])


# cluster color
colss <- readRDS(file = paste0(dataDir, "cluster.color.rds"))
cls_col <- setNames(colss$color, rownames(colss))

all_species <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")



##------------------------------------------------------------------
## Get the dex & markers (DEM) genes
slim.dex <- all.dex %>%
					filter(pct.1 >= 0.25 & pct.2 <= 0.35 & ratio_fc >= 1.5 & avg_logFC >= log(2.5)) %>%
					group_by(cluster, speciespair) %>% 
					mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
					filter(p_val_adj <= 0.01) %>%
					ungroup()
slim.mar <- all.markers %>%
					subset(pct.1 >= 0.25 & ratio_fc >= 1.25 & avg_logFC >= 0.25) %>%
					group_by(cluster, species) %>% 
					mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
					filter(p_val_adj <= 0.01) %>%
					ungroup()
demarkers_idx <- lapply(homo.order, function(cls) {
	cls_idx <- lapply(all_species, function(sp) {
		mar <- slim.mar$gene[slim.mar$species == sp & slim.mar$cluster == cls]
		idx <- which(slim.dex$species1 == sp & slim.dex$gene %in% mar & slim.dex$cluster == cls)
		idx
		}) %>% unlist() %>% unique()
	}) %>% unlist() %>% sort() %>% unique()
demarkers <- slim.dex[demarkers_idx, ]
saveRDS(demarkers, file = paste0(dataDir, "HPRC_DEMarkers.", ctp, ".rds"))





dex.up.base <- lapply(all_species, function(sp) {
	sub.dem <- filter(demarkers, species1 == sp) %>%
				group_by(cluster, gene) %>%
				summarize(nhits = n(), mfc = min(ratio_fc)) %>%
				filter(nhits == 3) %>%
				ungroup() %>%
				group_by(cluster) %>%
				arrange(desc(mfc)) %>%
				mutate(color = cls_col[cluster]) %>%
				ungroup() %>%
				mutate(species = sp)
	dem <- split(sub.dem, sub.dem$cluster) %>%
			.[homo.order] %>%
			do.call(rbind, .)
	dem
	}) %>% 
		do.call(rbind, .)
up.idx <- which(!duplicated(dex.up.base$gene))
up.genes <- dex.up.base$gene[up.idx]
up.cols <- dex.up.base$color[up.idx]



dex.down.base <- lapply(all_species, function(sp) {
	sub.dem <- filter(demarkers, species2 == sp) %>%
				group_by(cluster, gene) %>%
				summarize(nhits = n(), mfc = min(ratio_fc)) %>%
				filter(nhits == 3) %>%
				ungroup() %>%
				group_by(cluster) %>%
				arrange(desc(mfc)) %>%
				mutate(color = cls_col[cluster]) %>%
				ungroup() %>%
				mutate(species = sp)
	dem <- split(sub.dem, sub.dem$cluster) %>%
			.[homo.order] %>%
			do.call(rbind, .)
	dem
	}) %>% 
		do.call(rbind, .)
down.idx <- which(!duplicated(dex.down.base$gene))
down.genes <- dex.down.base$gene[down.idx]
down.cols <- dex.down.base$color[down.idx]



##------------------------------------------------------------------
## Get the label genes
set.seed(0)
dex.up.top <- split(dex.up.base[up.idx, ], dex.up.base$species[up.idx]) %>%
				.[all_species] %>%
				lapply(., function(x) {
					mm <- filter(x, mfc >= 2) %>%
								group_by(cluster) %>%
								top_n(., 2, mfc) %>% 
								.$gene 
					mm <- grep("^CTD-", mm, invert = TRUE, value = TRUE) %>%
							grep("^CTC-", ., invert = TRUE, value = TRUE) %>%
							grep("^RP[0-9]", ., invert = TRUE, value = TRUE) %>%
							grep("^AC[0-9]", ., invert = TRUE, value = TRUE) %>%
							grep("-AS1$", ., invert = TRUE, value = TRUE) %>%
							sample(., 10)
					mm
					}) 
dex.down.top <- split(dex.down.base[down.idx, ], dex.down.base$species[down.idx]) %>%
				.[all_species] %>%
				lapply(., function(x) {
					mm <- filter(x, mfc >= 2) %>%
								group_by(cluster) %>%
								top_n(., 2, mfc) %>% 
								.$gene 
					mm <- grep("^CTD-", mm, invert = TRUE, value = TRUE) %>%
							grep("^CTC-", ., invert = TRUE, value = TRUE) %>%
							grep("^RP[0-9]", ., invert = TRUE, value = TRUE) %>%
							grep("^AC[0-9]", ., invert = TRUE, value = TRUE) %>%
							grep("-AS1$", ., invert = TRUE, value = TRUE) %>%
							sample(., 6)
					mm
					}) 
dex.up.top$Human <- c(dex.up.top$Human, c("NMU", "ADGRG6", "MELTF"))
##dex.down.top$Human <- c(dex.down.top$Human, "TWIST1")


lbs <- list(c(dex.up.top$Human, dex.up.top$Rhesus),
			c(dex.down.top$Human, dex.down.top$Rhesus), 
			c(dex.up.top$Chimpanzee, dex.up.top$Marmoset),
			c(dex.down.top$Chimpanzee, dex.down.top$Marmoset))
label_genes <- list(setNames(c(lbs[[1]], lbs[[2]]), c(paste0("up_", lbs[[1]]), paste0("down_", lbs[[2]]))),
			setNames(c(lbs[[3]], lbs[[4]]), c(paste0("up_", lbs[[3]]), paste0("down_", lbs[[4]])))) 





##------------------------------------------------------------------
## Plot 
sp.cls <- paste0(rep(all_species, each = length(homo.order)), "|",
			rep(homo.order, times = length(all_species)))
avg_data <- readRDS(file = paste0(dataDir, "Expr_avg_ratio_by_species.rds"))[["fig1cluster"]]$avg

avg_data <- avg_data[, sp.cls] %>% 
				as.matrix() %>% 
				t() %>% scale() %>% t() %>% 
				MinMax(., min = -2.5, max = 2.5)


up_data <- avg_data[up.genes, ]; rownames(up_data) <- paste0("up_", rownames(up_data))
down_data <- avg_data[down.genes, ]; rownames(down_data) <- paste0("down_", rownames(down_data))
data_use <- rbind(up_data, down_data)
gene_split <- setNames(c(rep("Up", length(up.genes)), rep("Down", length(down.genes))) %>% factor(., levels = c("Up", "Down")), rownames(data_use))


## Set label genes
label_cols <- setNames(c(up.cols, down.cols), rownames(data_use)) %>%
				.[unlist(lapply(label_genes, names))]
source("../scripts/pfc.fun.R")
plot_dex_heatmap(data = data_use, split.by = "species", group.by = "cluster", label_genes = label_genes, label_col = label_cols, font_scale = 4, height_unit = 0.01, file_name = paste0("HPRC_DEX_", ctp), output_dir = outputdir, show_rownames = FALSE, color_limits = seq(-1, 2.5, 0.5), pdf_width = 18, pdf_height = 10, row_split = gene_split, heat_width = 6)












