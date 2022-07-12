library(Seurat)
library(tibble)
library(dplyr)
library(ggplot2)


all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
seu_list <- lapply(all_sps, function(sp) {
	readRDS(file = paste0("./load_files/", "Cluster_NNC_insp_", sp, ".rds"))
	}) %>%
		setNames(., all_sps)



allcls <- lapply(seu_list, function(x) unique(x$fig1cluster)) %>% unlist() %>% unique()



cluster.color <- readRDS(file = paste0("../data/", "cluster.color.rds")) %>%
						rownames_to_column("fig1cluster")
updata_name <- function(x) {
    xx <- gsub("L5 BCL11B SYT2 POU3F1", "L5 FEZF2 BCL11B POU3F1", x) %>%
            gsub("InN SST TH GABRQ", "InN SST HGF GABRQ", .) %>%
            gsub("InN LHX6 TH STON2", "InN LHX6 HGF STON2", .) %>%
            gsub("iAstro", "Astro", .) %>%
            gsub("fAstro", "Astro", .) %>%
            gsub("pAstro", "Astro", .) %>%
            gsub("rAstro", "Astro", .) %>%
            gsub("cOPC", "OPC", .) %>%
            gsub("nOligo", "Oligo", .) %>%
            gsub("eOligo", "Oligo", .) %>%
            gsub("mOligo", "Oligo", .) %>%
            gsub("lOligo", "Oligo", .) %>%
            gsub("aEndo", "Endo", .) %>%
            gsub("vEndo", "Endo", .) %>%
            gsub("cEndo", "Endo", .) %>%
            gsub("aaSMC", "SMC", .) %>%
            gsub("aSMC", "SMC", .) %>%
            gsub("vSMC", "SMC", .) %>%
            gsub("ABC", "VLMC",.) %>%
            gsub("L6b", "L6B",.)
    return(xx)
}
cluster.color$fig1cluster <- updata_name(cluster.color$fig1cluster)
cluster.color <- cluster.color %>%
					filter(fig1cluster %in% allcls)

cls_cols <- setNames(cluster.color$color, cluster.color$fig1cluster)

plist <- lapply(all_sps, function(sp){
	Idents(seu_list[[sp]]) <- "fig1cluster"
	p <- DimPlot(seu_list[[sp]], group.by = "fig1cluster", label = TRUE, cols = cls_cols, raster = TRUE, raster.dpi = c(512, 512)) +
			theme_void() +
			labs(title = sp) +
			theme(legend.position = "none")
	return(p)
	})

pdf(paste0("./report/Cluster_NNC_insp_UMAP.pdf"), width = 24, height = 6)
plot <- patchwork::wrap_plots(plist, nrow = 1, ncol = 4)
print(plot)
dev.off()
















