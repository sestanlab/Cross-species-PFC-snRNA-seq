source("../scripts/hip.fun.R") 
library(ggpubr)


hrpm <- readRDS(file = paste0(inputdir, "HRPM.all.seurat.expHVG.1500.slim.rds"))
hrpm@meta.data$cluster <- hrpm@meta.data$fig2cluster
hrpm@meta.data$cluster[hrpm@meta.data$fig1cluster %in% c("Micro C1QB CD83", "Micro C1QB P2RY12")] <- "Micro"

DimFig(hrpm, file_name = "FOXP2_UMAP_species", group.by = "cluster", plot.scale = 1, pt.size = 0.01)



plot_data <- data.frame(exp = hrpm$RNA@data["FOXP2", ], 
                    species = hrpm@meta.data$species,
                    cluster = hrpm@meta.data$fig2cluster,
                    xaxis = hrpm$umap@cell.embeddings[, 1], 
                    yaxis = hrpm$umap@cell.embeddings[, 2], 
                    stringsAsFactors = FALSE)
me <- max(plot_data$exp)
bb <- 0
xrange <- c(min(plot_data$xaxis)- bb, 
            max(plot_data$xaxis) + bb)
yrange <- c(min(plot_data$yaxis)- bb, 
            max(plot_data$yaxis) + bb)


all_sps <- c("Human", "Rhesus", "Pig", "Mouse")
umap_list <- lapply(all_sps, function(sp) {
            sub_data <- subset(plot_data, species == sp)
            p <- ggplot(sub_data, aes(x = xaxis, y = yaxis, color = exp)) +
                geom_point(size = 0.01) +
                scale_color_gradientn(colors = c("lightgrey", "darkred"), limits = c(0, me)) +
                theme_classic() + 
                theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "bottom")
            p <- p + 
                        xlim(limits = xrange) +
                        ylim(limits = yrange)
            return(p)
        })

umap_list <- umap_list %>%
          lapply(., function(x) x + theme(plot.margin = unit(c(-0.1, -0.1, -0.1, -0.1), "inch")))

jpeg(paste0(outputdir, "FOXP2_UMAP_exp.jpeg"), width = 32, height = 8, units = "in", res = 300)
plot_grid(plotlist = umap_list, nrow = 1, ncol = 4) %>% print()
dev.off() 








