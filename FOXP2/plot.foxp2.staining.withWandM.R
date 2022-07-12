library(dplyr)
library(ggplot2)
library(ape)
library(phytools)
library(tibble)
library(Seurat)

inputdir <- paste0("./load_files/staining/")
outputdir <- paste0("./report/staining/")


tree.txt <- readLines(paste0(inputdir, "latin.sp.withWandM.nwk"))
newTree <- ape::read.tree(text=tree.txt)


rt_idx <- c(52, 58, 94, 86, 87, 88, 89, 81, 82, 83, 78, 75, 68, 66, 63, 99)
for (i in rt_idx){
	newTree <- rotate(newTree, i)
}

newTree$tip.label <- gsub("_", " ", newTree$tip.label)
pdf(paste0(outputdir, "species.tree.label.pdf"), width = 20, height = 20)
plotTree(newTree,node.numbers=T)
axisPhylo()
dev.off()


pdf(paste0(outputdir, "species.tree.pdf"), width = 6, height = 8)
plot(newTree, cex = .75)
axisPhylo()
dev.off()



## Plot the expression
is_tip <- newTree$edge[,2] <= length(newTree$tip.label)
tips <- newTree$tip.label[newTree$edge[is_tip, 2]] %>% rev()
new_tips <- tips %>%
				gsub("Cuniculus paca", "Agouti paca", .)

exp_data <- read.table(file = paste0(inputdir, "foxp2_staining.update.withWandM.1104.2020.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) %>%
				select(-c(Order, CommonName)) %>%
				column_to_rownames("Species") %>%
				as.matrix()



plot_data <- reshape2::melt(exp_data) %>%
				setNames(., c("species", "layer", "exp_ratio")) %>%
				mutate(species = as.character(species)) %>%
				mutate(species = factor(species, levels = rev(new_tips)))

pp <- ggplot(data = plot_data, mapping = aes_string(x = "layer", y = "species", size = "exp_ratio")) +
		geom_point(color = "black", stroke = 0.25, shape = 21, fill = "violet") +
		theme_classic() + 
		scale_radius(range = c(0, 5), limits = c(0, 100)) + 
		RotatedAxis() + 
		theme(legend.position = "bottom")

pdf(paste0(outputdir, "FOXP2_staing_results.pdf"), width = 3.5, height = 7)
print(pp) 
dev.off() 













