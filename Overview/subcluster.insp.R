args <- commandArgs(trailingOnly = TRUE)
source("../scripts/pfc.fun.R")
source("./inte.fun.R")


library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 28*1000*1024^2)


gp <- args[1]
sp <- args[2]
nfeatures <- ifelse(gp == "CGE", 1200, 1000)
sel_gps <- list(`L23IT` = "L2-3 IT",
				CGE = c("VIP", "ADARB2 KCNG1", "LAMP5 RELN", "LAMP5 LHX6"),
				Astro = "Astro",
				Immune = c("Micro", "Immune"), 
				Oligo = c("OPC", "Oligo")) %>%
				.[[gp]]
seu <- readRDS(file = paste0(inputdir, "RawSeurat_", sp, ".rds"))
seu <- subset(seu, mres %in% sel_gps)



## Prepare HVGs
seu_list <- SplitObject(seu, split.by = "samplename")
hvg <- lapply(seu_list[setdiff(names(seu_list), "CJB1680")], function(x) FindVariableFeatures(x, nfeatures = (nfeatures + 500))) %>%
				SelectIntegrationFeatures(., nfeatures = nfeatures)
if (gp == "Oligo"){
	mars <- readRDS(file = "../SpecclsMars/load_files/Oligo.speccls.putative.markers.rds") %>% unlist() %>% unique()
    idsymbol <- unlist(mclapply(mars, function(i) ifelse(length(grep(paste0("-",i,"$"), rownames(seu), ignore.case= TRUE))==1, rownames(seu)[grep(paste0("-",i,"$"),rownames(seu), ignore.case= TRUE)], "empty")))
    idsymbol <- idsymbol[idsymbol != "empty"]
    hvg <- union(hvg, idsymbol)
}
rm(seu_list)


cbn <- SeuInte(object = seu, hvg = hvg, file_name = paste0("Subcluster_insp_", gp, "_", sp), input_dir = inputdir, do_cluster = FALSE, merge.order = NULL)


plist <- DimFig(cbn, group.by = c("fig1cluster", "mres", "samplename"), file_name = "AA", plot.scale = 1, output.ggplot = TRUE)
jpeg(paste0(outputdir, "Subcluster_insp_", gp, "_", sp, "_Idents.jpeg"), width = 2 * 9, height = 2 * 9, units = "in", res = 300)
plot_grid(plotlist = plist[c(1, 2, 5, 7)], nrow = 2, ncol = 2) %>% print()
dev.off()	



