## Plot the expression of SST across layers
source("../scripts/pfc.fun.R")
source("../MF6_InN/inn.fun.R")


###---------------------------------------------------------------------------------
## Visualize in primates (HPRC)
pfc <- readRDS(paste0(dataDir, file = "Final_PFC_HPRC.InN.rm1680.07102021.rds"))
sel.cls <- c("InN PVALB GRIN2C", "InN PVALB PDE3A")
seu <- subset(pfc, fig1cluster %in% sel.cls)


exp.thre <- 0.25

all.sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
mar_df <- lapply(all.sps, function(sp) {
	sub.seu <- subset(seu, species == sp)
	Idents(sub.seu) <- "fig1cluster"
	sp_df <- lapply(sel.cls, function(cls){
					m1 <- FindMarkers(sub.seu, ident.1 = cls, pct.1 = exp.thre, only.pos = TRUE, logfc.threshold = 0.4) %>%
						rownames_to_column("gene") %>%
						mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
						filter(ratio_fc >= 1.1) %>%
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
		.[sel.cls] %>%
		unlist() %>%
		unique()
	})


## Check the enrichment of gene families in DEGs.
all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
sp_cls <- paste0(rep(all_sps, each = 2), "|", sel.cls) 
ratios <- readRDS(file = paste0(dataDir, "Expr_avg_ratio_by_species.rds"))[["fig1cluster"]]$ratio
exp_genes <- rownames(ratios)[apply(ratios[, sp_cls], 1, max) >= exp.thre]


load(file = "~/project/public_data/HGNC/HGNC_gset.09152021.Rdata") ##glist, gtable
fgset <- lapply(glist, function(x) intersect(x, exp_genes))
fgset <- fgset[sapply(fgset, length) >= 5]


enr.sp.list <- lapply(mars, function(xx) {
	henr <- gset_enrich(universe = exp_genes, interesting_genes = xx, gset = fgset, min_nodeSize = 6)
	henr
	})

gfs <- lapply(enr.sp.list, function(x) x$GFname[x$FDR <= 0.1])
gfs


## Barplot showing the significance
all_gfs <- unlist(gfs, use.names = FALSE) %>% unique()
pdata <- lapply(all_sps, function(sp) {
	df <- data.frame(row.names = all_gfs,
			gset = all_gfs,
			species = sp,
			FDR = 1,
			stringsAsFactors = FALSE)
	henr <- enr.sp.list[[sp]]
	sh_gf <- intersect(rownames(df), rownames(henr))
	df[sh_gf, "FDR"] <- henr[sh_gf, "FDR"]
	df
	}) %>%
		do.call(rbind, .) %>%
		mutate(logFDR = -log10(FDR))

pdata$gset <- gsub("Glutamate ionotropic receptor NMDA type subunits", "GRIN", pdata$gset) %>%
					gsub("Non-clustered protocadherins", "ncProtocadherins", .)


sp_cols <- setNames(c("#FF420E","#4CB5F5","#89DA59","#FFBB00"), c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
p <- ggplot(pdata, aes_string(x = "gset", y = "logFDR", fill = "species")) +
		geom_bar(stat = "identity", size = 0.2, position = position_dodge()) +
		scale_fill_manual(values = sp_cols) +
		theme_cowplot() + 
		RotatedAxis() +
		scale_x_discrete(limits = c("GRIN", "ncProtocadherins")) +
		geom_hline(yintercept = 1, size = 0.4, linetype = "dashed") +
		theme(axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.text.x = element_text(size = rel(0.6)), axis.text.y = element_text(size = rel(1)), axis.title = element_blank(), legend.position = "none")


pdf(paste0(outputdir, "PVALB_CHC_gset_enriched_in_markers.pdf"), width = 5, height = 5)
print(p)
dev.off()


###---------------------------------------------------------------------------------
## Visualize in primates (HPRC)
vis_gfs <- unlist(gfs) %>% unique() ##"Glutamate ionotropic receptor NMDA type subunits" & "Non-clustered protocadherins"

enriched_gfgs <- lapply(vis_gfs, function(gf) {
	gf_genes <- lapply(sel.cls, function(cls) {
		mm <- mar_df %>%
				do.call(rbind, .) %>%
				filter(cluster == cls) %>%
				group_by(gene) %>%
				summarize(nhits = n()) %>%
				ungroup() %>%
				filter(gene %in% fgset[[gf]]) %>%
				arrange(desc(nhits)) %>%
				.$gene
		mm
		})
	unlist(gf_genes) %>% unique()
	}) %>% 
		setNames(., vis_gfs)

all_gfgs <- lapply(vis_gfs, function(gf) c(enriched_gfgs[[gf]], setdiff(fgset[[gf]], enriched_gfgs[[gf]]))) %>% setNames(., vis_gfs)


gene_plot <- unlist(all_gfgs) %>% unique()


sp_cls <- paste0(rep(all_sps, each = 2), "|", sel.cls) 
pdata <- seu@meta.data[, c("fig1cluster", "species")] %>%
			cbind(., t(seu$RNA@data[gene_plot, ])) %>%
			mutate(newcls = paste0(species, "|", fig1cluster)) %>%
			mutate(newcls = factor(newcls, levels = sp_cls))

plist <- lapply(gene_plot, function(gg) {
	p <- ggplot(pdata, aes_string(x = "newcls", y = gg, fill = "species")) + 
			geom_violin(scale = "width", size = 0.1, adjust = 2,trim =TRUE) + 
			theme_cowplot() + 
			labs(y = gg) +
			scale_fill_manual(values = c("#FF420E", "#4CB5F5", "#89DA59", "#FFBB00") %>% setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset"))) + 
			theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.9), angle = 0, hjust = 0.5, vjust = 1), axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = rel(0.8)), axis.line.y = element_line(size = 0.2), axis.ticks.y = element_line(size = 0.2), legend.position = "none", plot.margin = unit(c(-0.05, 0, -0.05, 0), "in"))
	p
	})

plist[[length(gene_plot)]] <- plist[[length(gene_plot)]] +
			theme(axis.text.x=element_text(size = rel(0.9), angle = 45, hjust = 1, vjust = 1), axis.line.x = element_line(size = 0.2), axis.ticks.x = element_line(size = 0.2)) 
out_name <- gsub("Glutamate ionotropic receptor NMDA type subunits", "GRIN", gf) %>%
					gsub("Non-clustered protocadherins", "ncProtocadherins", .)
pdf(paste0(outputdir, "PVALB_CHC_marker_sep.pdf"), width = 5, height = 10)
patchwork::wrap_plots(plist, nrow = length(gene_plot), ncol = 1)
dev.off()




##-----------------------------------------------------------
## Heatmaps showing all the markers
seu@meta.data$sp_cls <- paste0(seu@meta.data$species, "|", seu@meta.data$fig1cluster)
sp_cls <- paste0(rep(all_sps, each = 2), "|", sel.cls) 
cell_ord <- lapply(sp_cls, function(x) colnames(seu)[seu@meta.data$sp_cls == x]) %>% unlist()


plot_meta <- seu@meta.data[cell_ord, c("species", "fig1cluster")] %>%
				mutate(species = factor(as.character(species), levels = c("Human", "Chimpanzee", "Rhesus", "Marmoset")))
rna_data <- seu$RNA@data[, cell_ord]


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
##-2.5, -2, 
## Get the range of the heatmap legends
color_breaks <- seq(-1.5, 2, 0.25)

anno_cols <- list(species = c("#FF420E","#4CB5F5","#89DA59","#FFBB00") %>%
			setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset")),
				fig1cluster = c("#8c2905", "#b51217") %>% setNames(., sel.cls))

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


pdf(paste0(outputdir, "PVALB_CHC_subtype_markers_heat.pdf"), width = 12, height = 8.5)
draw(htlist)
dev.off() 













