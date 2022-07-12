source("../scripts/pfc.fun.R")

ctype <- "InN"
hprc <- readRDS(file = paste0(dataDir, "Final_PFC_HPRC.", ctype, ".07092021.rds"))
meta <- hprc@meta.data 
cell_coordinates <- hprc$umap@cell.embeddings



xlim <- c(floor(min(cell_coordinates[, 1])), ceiling(max(cell_coordinates[, 1])))
ylim <- c(floor(min(cell_coordinates[, 2])), ceiling(max(cell_coordinates[, 2])))




jpeg(paste0(outputdir, "HPRC.", ctype, ".umap.bycluster.UMAP.jpeg"), width = 6, height = 6, units = "in", res = 300)
group.by <- "fig1cluster"
allcls <- levels(as.factor(as.character(meta[, group.by])))
cluster.color <- readRDS(file = paste0(dataDir, "cluster.color.rds")) %>%
						rownames_to_column("fig1cluster") %>%
						.[.$fig1cluster %in% allcls, ]
color_vec <- setNames(cluster.color$color, cluster.color$fig1cluster)

plot(cell_coordinates, asp = 1, cex = 0.1, pch = 16, col= color_vec[as.character(meta[, group.by])], frame.plot= FALSE, axes = FALSE, xlim = xlim, ylim = ylim)
dev.off()






jpeg(paste0(outputdir, "HPRC.", ctype, ".umap.byspecies.UMAP.jpeg"), width = 6, height = 6, units = "in", res = 300)
split.by <- "species"
sp.cols <- setNames(c("#FF420E","#4CB5F5","#89DA59","#FFBB00"), c("Human", "Chimpanzee", "Rhesus", "Marmoset"))

set.seed(0)
new_ord <- sample(1:nrow(cell_coordinates), nrow(cell_coordinates)) 


plot(cell_coordinates[new_ord, ], asp = 1, cex = 0.1, pch = 16, col= sp.cols[as.character(meta[, split.by])][new_ord], frame.plot= FALSE, axes = FALSE, xlim = xlim, ylim = ylim)
dev.off()






