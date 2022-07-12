## Plot the markers that differentiate L6IT clusters in heatmap
source("../scripts/pfc.fun.R")


ctp <- "L6_IT_2"
seu <- readRDS(file = paste0(inputdir, "Recls_MNN_", ctp, ".mnn.rds"))


##-----------------------------------------------------------------
## Plot subtypes across species
cls_cols <- c("#009E73", "#CC79A7", "#0072B2", "#D55E00") %>%
				setNames(., c("L6 OPRK1 SMYD1 ADAMTS17", "L6 OPRK1 SMYD1 SNTB1", "L6 OPRK1 SMYD1 KCND2", "L6 OPRK1 SMYD1 TSHZ2"))
seu$fig1cluster <- as.character(seu$fig1cluster)
p1 <- DimFig(seu, group.by = "fig1cluster", pt.size = 0.8, output.ggplot = TRUE, file_name = "test", split.by = "species", split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), cols = cls_cols)[[3]]+
					theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), plot.title = element_blank())
jpeg(paste0(outputdir, "SF5_ExN_Heter_cluster_", ctp, "_UMAP.jpeg"), width = 12, height = 3, unit = "in", res = 300)
print(p1)
dev.off()		



##-----------------------------------------------------------------
## Extract the markers within each species
seu.list <- SplitObject(seu, split.by = "species")
all.cls <- c("L6 OPRK1 SMYD1 ADAMTS17", "L6 OPRK1 SMYD1 SNTB1", "L6 OPRK1 SMYD1 KCND2", "L6 OPRK1 SMYD1 TSHZ2")
res_df <- lapply(names(seu.list), function(sp) {
	Idents(seu.list[[sp]]) <- "fig1cluster"
	sp_res <- lapply(all.cls, function(cls) {
		res <- FindMarkers(seu.list[[sp]], ident.1 = cls, min.pct = 0.2, max.cells.per.ident = 500, only.pos = TRUE, logfc.threshold = 0.3) %>%
				rownames_to_column("gene") %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(species = sp, cluster = cls)
		res
		}) %>%
		do.call(rbind, .)
	sp_res
	}) %>%
		do.call(rbind, .)


mar_df <- res_df %>%
			group_by(cluster, species) %>%
			filter(p_val_adj <= 0.01) %>%
			ungroup() 
mar_df <- mar_df[!(mar_df$species == "Marmoset" & mar_df$cluster == "L6 OPRK1 SMYD1 ADAMTS17"), ]

mars <- split(mar_df$gene, mar_df$species) %>%
			lapply(., unique)



##-----------------------------------------------------------------
## Plot the Heatmap that visualize the distinctions of gene expression
seu@meta.data$sp_cls <- paste0(seu@meta.data$species, "|", seu@meta.data$fig1cluster)
cls_ord <- rep(c("Human", "Chimpanzee", "Rhesus", "Marmoset"), each = 4) %>%
			paste0(., "|", rep(c("L6 OPRK1 SMYD1 ADAMTS17", "L6 OPRK1 SMYD1 SNTB1", "L6 OPRK1 SMYD1 KCND2", "L6 OPRK1 SMYD1 TSHZ2"), times = 4))
cells <- lapply(cls_ord, function(x) colnames(seu)[seu@meta.data$sp_cls == x]) %>% unlist()



plot_meta <- seu@meta.data[cells, c("species", "fig1cluster")] %>%
				mutate(species = factor(as.character(species), levels = c("Human", "Chimpanzee", "Rhesus", "Marmoset")))
rna_data <- seu$RNA@data[, cells]


data <- lapply(c("Human", "Chimpanzee", "Rhesus", "Marmoset"), function(sp) {
			genes <- mars[[sp]]
			mat <- rna_data[genes, ] %>% 
				as.matrix() %>%
				t() %>% scale() %>% t()
			rownames(mat) <- paste0(sp, "_", rownames(mat))
			mat
		}) %>%
			do.call(rbind, .)

data <- MinMax(data, min = -1.5, max = 2)
color_breaks <- seq(-1.5, 2, 0.25)

anno_cols <- list(species = c("#FF420E","#4CB5F5","#89DA59","#FFBB00") %>%
			setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset")),
				fig1cluster = gg_color_hue(4) %>%
			setNames(., c("L6 OPRK1 SMYD1 ADAMTS17", "L6 OPRK1 SMYD1 KCND2", "L6 OPRK1 SMYD1 SNTB1", "L6 OPRK1 SMYD1 TSHZ2")))

column_ha <- HeatmapAnnotation(df = plot_meta, col = anno_cols, annotation_height = unit(c(0.01, 0.01), "in"))
htlist <- Heatmap(data, name = "scaled_expr", 
                    col = colorRamp2(color_breaks, viridis(20)[c(1:11, 14, 16, 18, 20)]), na_col = "white",
                    row_split = factor(extract_field(rownames(data), 1, "_"), levels = c("Human", "Chimpanzee", "Rhesus", "Marmoset")), 
                    column_split = plot_meta[, "species"], 
                    cluster_rows = FALSE, cluster_columns = FALSE, 
                    row_title_rot = 90, column_title = NULL, row_title_gp = gpar(fontsize = 12), 
                    show_column_names = FALSE, show_row_names = FALSE, column_names_gp = gpar(fontsize = 10),column_names_rot = 45, 
                    width = unit(8, "in"),
                    top_annotation = column_ha, 
                    heatmap_legend_param = list(title = "Scaled expr", color_bar = "continuous", at = color_breaks),
                    use_raster = TRUE, raster_quality = 3)


pdf(paste0(outputdir, "Heter_markers_", ctp, "_heat.pdf"), width = 12, height = 8.5)
draw(htlist)
dev.off()


















