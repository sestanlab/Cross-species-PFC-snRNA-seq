source("../scripts/pfc.fun.R")


seu <- readRDS(file = paste0(inputdir, "Inte_Geirsdottir_Micro_v2.slim.rds"))



##------------------------------------------------------------------------------------
## Visualize the Identity
pdata <- seu@meta.data[, c("stage", "fig1cluster")] %>%
			cbind(., seu$umap@cell.embeddings)
pdata$cluster <- ifelse(is.na(pdata$fig1cluster), "unknown", pdata$fig1cluster)
set.seed(42)
pdata <- pdata[sample(1:nrow(pdata)), ]



all_cls <- levels(as.factor(pdata$cluster)) %>%
			setdiff(., "unknown")
cls_cols <- c("#00C094", "#C49A00", "#53B400", "#F8766D", "#FB61D7", "#A58AFF", "#00B6EB")[2:5] %>%
				setNames(., all_cls)
library(ggrastr)
p1 <- ggplot(subset(pdata, cluster != "unknown"), aes_string(x = "UMAP_1", y = "UMAP_2", color = "cluster")) +
			rasterise(geom_point(size = 0.4, shape = 16), dpi = 300, scale = 1) +
	        coord_equal(ratio = 1) + 
	        theme_classic() + 
	        scale_color_manual(values = cls_cols) + 
	        theme(legend.position = "none",
	                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
	                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
	                    plot.title = element_blank())
p2 <- ggplot(pdata, aes_string(x = "UMAP_1", y = "UMAP_2", color = "stage")) +
			rasterise(geom_point(size = 0.4, shape = 16), dpi = 300, scale = 1) +
	        coord_equal(ratio = 1) + 
	        theme_classic() + 
	        scale_color_manual(values = c(adult = "#FF420E", Geirsdottir = "#762a83")) + 
	        theme(legend.position = "none",
	                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
	                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
	                    plot.title = element_blank())
pdf(paste0(outputdir, "SF.inte.geir.vis.pdf"), width = 5 * 2, height = 5)
plot <- patchwork::wrap_plots(plotlist = list(p2, p1), nrow = 1, ncol = 2, guides = "collect") & theme(legend.position = "none")
print(plot)
dev.off()




##------------------------------------------------------------------------------------
## Visualize gene expression
plot_hex_byregion_v2 <- function(sce, feature, action = "mean", cols = colorRampPalette(c("lightgrey", "red"))(3), split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), assay = "RNA", qq = 1){
    if (!feature %in% rownames(sce[[1]][[assay]])) {
        stop("Specify a valid assay.")
    }

    ##First prepare the data list
    all_regions <- ifelse_check(is.null(split.order), names(sce), split.order)
    exp_list <- lapply(all_regions, function(reg) as.numeric(GetAssayData(sce[[reg]], slot = "data")[feature, ])) %>% setNames(., all_regions)
    plot.data <- lapply(all_regions, function(reg) {
        x <- exp_list[[reg]]
        out <- sce[[reg]]@misc$hexbin[[2]]
        cID <- sce[[reg]]@misc$hexbin[[1]]
        df <- data.frame(out, 
                exp = schex:::.make_hexbin_function(x, action, cID), 
                region = reg, 
                stringsAsFactors = FALSE, check.names = FALSE)
        return(df)
        }) %>% do.call(rbind, .)
    

    maxexp <- quantile(plot.data$exp, qq) %>% 
    			max(., 1)
    plot.data$exp[plot.data$exp > maxexp] <- maxexp

    regp_list <- list()
    for (reg in all_regions){
        sub.data <- subset(plot.data, region == reg)
        ##cur.cols <- colorRampPalette(colors = cols)(nbks)[1:max(sub.data$scale.exp)]
        regp_list[[reg]] <- ggplot(data = sub.data) +
                    rasterise(geom_hex(mapping = aes_string(x = "x", y = "y", fill = "exp"), stat = "identity"), dpi = 300, scale = 1) + 
                    theme_classic() + 
                    labs(title = reg) + 
                    coord_fixed() + 
                    scale_fill_gradientn(colors = cols, limits = c(0, maxexp)) + 
                    theme(legend.position = "bottom", axis.line = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), legend.title = element_blank())
    }    
    return(regp_list)
}

seu_hexlist <- make_hexbin_byregion(seu, nbins = 50, dimension_reduction = "UMAP", split.by = "stage")



sel_genes <- c("CCL3", "CCL4", "EGR3", "SLC11A1", "SCIN", "GPNMB")
plist <- lapply(sel_genes, function(gg) {
	pp <- plot_hex_byregion_v2(sce = seu_hexlist, feature = gg, action = "mean", cols = colorRampPalette(c("lightgrey", "red"))(3), split.order = c("Geirsdottir", "adult"), assay = "RNA", qq = 0.98)
	pp[[1]] <- pp[[1]] +
			theme(plot.title = element_blank())
	pp[[2]] <- pp[[2]] +
			theme(plot.title = element_blank())
	return(pp)
	}) %>%
		do.call(c, .)


pdf(paste0(outputdir, "SF.inte.geir.vis.exp.CCL3_GLDN.pdf"), width = 3.5 * 2 + 0.5, height = 3.5*length(sel_genes))
plot <- patchwork::wrap_plots(plotlist = plist, nrow = length(sel_genes), ncol = 2, guides = "collect") & theme(legend.position = "none")
print(plot)
dev.off()


## Median nGenes
#adult          826
#Geirsdottir          326
## Mean nGenes
#adult          854
#Geirsdottir          361





