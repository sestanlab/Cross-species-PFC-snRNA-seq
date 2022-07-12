args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")

gp <- args[1]


all_sps <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
seu_list <- lapply(all_sps, function(sp) {
	seu <- readRDS(file = paste0(inputdir, "Subcluster_insp_", gp, "_", sp, ".rds"))
	if (gp %in% c("CGE")){
		seu <- subset(seu, mres == "LAMP5 RELN") %>%
					RunUMAP(., dims = 1:35)
	} else if (gp %in% "Oligo") {
		seu <- subset(seu, mres == "OPC") %>%
					RunUMAP(., dims = 1:35)
	}
	return(seu)
	}) %>%
		setNames(., all_sps)


pdata <- lapply(all_sps, function(sp) {
	data <- data.frame(cell = rownames(seu_list[[sp]]@meta.data),
					cluster = as.character(seu_list[[sp]]@meta.data$fig1cluster),
					species = sp,
					stringsAsFactors = FALSE) %>%
				cbind(., seu_list[[sp]]$umap@cell.embeddings)
	data
	}) %>%
		do.call(rbind, .)


library(ggrastr)
set.seed(42)
all_cls <- levels(as.factor(pdata$cluster))
cls_cols <- setNames(gg_color_hue(length(all_cls)), all_cls)
if (gp == "Immune"){
	cls_cols <- c("#00C094", "#C49A00", "#53B400", "#F8766D", "#FB61D7", "#A58AFF", "#00B6EB") %>%
				setNames(., all_cls)
}
plist <- lapply(all_sps, function(sp) {
	subdata <- filter(pdata, species == sp)
	p <- ggplot(subdata, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
				rasterise(geom_point(size = ifelse(gp == "CGE", 0.5, 0.2), shape = 16), dpi = 300, scale = 1) +
				theme_void()+
				scale_color_manual(values = cls_cols)
	p
	})

pdf(paste0(outputdir, "Species-spec-clusters_vis_", gp, ".pdf"), width = 4, height = 4 * 4, useDingbats = FALSE)
plot <- patchwork::wrap_plots(plist, nrow = 4, ncol = 1, guides = "collect") & theme(legend.position = "bottom")
print(plot)
dev.off()



