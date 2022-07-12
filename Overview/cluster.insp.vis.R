args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")

sp <- args[1]

seu <- readRDS(file = paste0(inputdir, "Cluster_insp_", sp, ".rds"))


plist <- DimFig(seu, group.by = c("fig1cluster", "mres", "samplename"), file_name = "AA", plot.scale = 1, output.ggplot = TRUE)
jpeg(paste0(outputdir, "Cluster_insp_", sp, "_Idents.jpeg"), width = 2 * 15, height = 2 * 15, units = "in", res = 300)
plot_grid(plotlist = plist[c(1, 2, 5, 7)], nrow = 2, ncol = 2) %>% print()
dev.off()


