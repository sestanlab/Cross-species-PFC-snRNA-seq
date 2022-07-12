### R-4.1.0
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)



setwd("./adult_human_immune")
Convert("global.h5ad", dest = "h5seurat", overwrite = TRUE)
immune <- LoadH5Seurat("global.h5seurat")


## Save as RDS file
saveRDS(immune, file = paste0("Human_adult_immune_cross-organs.rds"))


## Subset to macrophages across tissues
macro <- subset(immune, Majority_voting_CellTypist_high == "Macrophages")


## Only keep the organs with sufficient cells
org_keep <- table(macro$Organ) %>% .[.>= 20] %>% names()
macro <- subset(macro, Organ %in% org_keep)
macro@meta.data$Organ <- as.character(macro@meta.data$Organ)
saveRDS(macro, file = paste0("Human_adult_immune_macrophage_cross-organs.rds"))



## Do the violin plots
features <- c("C1QC", "FOXP2")
org_ord <- c("BLD","LIV","LNG","OME","DUO","ILE","JEJEPI", "JEJLP", "BMA", "MLN", "SPL", "LLN")
pdata <- macro@meta.data[, "Organ", drop = FALSE] %>%
            cbind(., t(macro$RNA@data[features, ])) %>%
			mutate(Organ = factor(as.character(Organ), levels = org_ord))

plist <- lapply(features, function(gg) {
    p <- ggplot(pdata, aes_string(x = "Organ", y = gg, fill = "Organ")) + 
            geom_violin(scale = "width", size = 0.1, adjust = 2,trim =TRUE) + 
            theme_classic() + 
            labs(y = gg) +
            theme(axis.title.x = element_blank(), axis.title.y = element_text(size = rel(0.9), angle = 0, hjust = 1, vjust = 0.5), axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size = rel(0.8)), axis.line.y = element_line(size = 0.2), axis.ticks.y = element_line(size = 0.2), legend.position = "none", plot.margin = unit(c(-0.05, 0, -0.05, 0), "in"))
    p
    })

plist[[length(features)]] <- plist[[length(features)]] +
            theme(axis.text.x=element_text(size = rel(0.9), angle = 45, hjust = 1, vjust = 1), axis.line.x = element_line(size = 0.2), axis.ticks.x = element_line(size = 0.2)) 
pdf(paste0("../report/", "Human_adult_macro_FOXP2_expr_violin.pdf"), width = 4, height = 2)
patchwork::wrap_plots(plist, nrow = length(features), ncol = 1)
dev.off()


















