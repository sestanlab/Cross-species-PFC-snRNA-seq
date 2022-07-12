## Visualize gene expression on HEX plots
source("../scripts/pfc.fun.R")
library(schex)


pfc <- readRDS(paste0(dataDir, file = "Final_PFC_HPRC.InN.07102021.rds"))
pfc_hexlist <- make_hexbin_byregion(pfc, nbins = 80, dimension_reduction = "UMAP", split.by = "species")


cols <- colorRampPalette(c("lightgrey", "darkred"))(10)[c(1, 5, 7, 8, 9, 10)]
genes <- c("TH", "GCH1", "DRD2", "SLC18A2", "SLC47A1", "SLC6A3")

for (gene in genes){
	plist <- plot_hex_byregion(sce = pfc_hexlist, feature = gene, action = "mean", cols = cols, split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
	jpeg(paste0(outputdir, "HPRC_InN_hex_", gene, ".jpeg"), width = 24, height = 6, res = 300, units = "in")
	patchwork::wrap_plots(plist, nrow = 1, ncol = 4) %>% print()
	dev.off()
}







