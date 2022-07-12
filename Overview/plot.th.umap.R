source("../scripts/pfc.fun.R")
library(schex)



pfc <- readRDS(paste0(dataDir, file = "Final_PFC_HPRC.InN.rm1680.07102021.rds"))
pfc_hexlist <- make_hexbin_byregion(pfc, nbins = 80, dimension_reduction = "UMAP", split.by = "species")



plot.data <- lapply(all_sps, function(sp) {
        genes <- c("TH")
        sce <- seu_list[[sp]]
		mapgenes <- unlist(mclapply(genes, function(i) ifelse(length(grep(paste0("-",i,"$"), rownames(sce$RNA), ignore.case= TRUE))==1, rownames(sce$RNA)[grep(paste0("-",i,"$"),rownames(sce$RNA), ignore.case= TRUE)], "empty")))
		exp <- as.numeric(GetAssayData(sce, slot = "data")[mapgenes, ]) %>% 
			setNames(., colnames(sce)) %>%
			.[colnames(pfc_hexlist[[sp]])]
        out <- pfc_hexlist[[sp]]@misc$hexbin[[2]]
        cID <- pfc_hexlist[[sp]]@misc$hexbin[[1]]
        df <- data.frame(out, 
                exp = schex:::.make_hexbin_function(exp, "mean", cID), 
                species = sp, 
                stringsAsFactors = FALSE, check.names = FALSE)
        return(df)
        }) %>% 
		do.call(rbind, .)
me <- max(plot.data$exp)

plist <- lapply(all_sps, function(sp) {
	p <- ggplot(data = subset(plot.data, species == sp)) +
            ggrast::rasterise(geom_hex(mapping = aes_string(x = "x", y = "y", fill = "exp"), stat = "identity"), dpi = 300, scale = 1) +
            theme_classic() + 
            labs(title = sp) + 
            coord_fixed() + 
            scale_fill_gradientn(colors = colorRampPalette(c("lightgrey", "red"))(3), limits = c(0, me)) + 
            theme(legend.position = "bottom", axis.line = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), legend.title = element_blank())
    p
	})
pdf(paste0(outputdir, "TH_expression_nonLiftover_UMAP.pdf"), width =24, height = 6.5)
plot <- patchwork::wrap_plots(plist, nrow = 1, ncol = 4, guides = "collect") & theme(legend.position = "bottom")
print(plot)
dev.off() 




















