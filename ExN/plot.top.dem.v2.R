### Cluster the disease genes
source("../scripts/pfc.fun.R")


load(file = paste0(inputdir, "HPRC.ExN.topDEM.mres.Rdata"))

##exl_up, exl_down


## Set object
resolution <- "mres"
object <- readRDS(file = paste0(dataDir, "DotHeatPlot.pre.rds"))
avgs <- object[[resolution]]$avg
ratios <- object[[resolution]]$ratio


cls_transfer <- c("L2-3 IT","L3-5 IT-1","L3-5 IT-2","L3-5 IT-3","L6 IT-1","L6 IT-2","L5 ET","L5-6 NP","L6 CT","L6B") %>%
            setNames(., c("L2_3_IT","L3_5_IT_1","L3_5_IT_2","L3_5_IT_3","L6_IT_1","L6_IT_2","L5_PT","L5_6_NP","L6_CT","L6b"))
cls_order <- c("L2-3 IT","L3-5 IT-1","L3-5 IT-2","L3-5 IT-3","L6 IT-1","L6 IT-2","L5 ET","L5-6 NP","L6 CT","L6B")
nametrans <- function(x) {
    for (ii in 1:length(cls_transfer)){
        x <- gsub(names(cls_transfer)[ii], setNames(cls_transfer[ii], NULL), x)
    }
    return(x)
}


colnames(avgs) <- nametrans(colnames(avgs))
colnames(ratios) <- nametrans(colnames(ratios))


## Filter genes
sp <- "Human"
source("./exp.fun.R")
#for (sp in c("Human", "Chimpanzee", "Rhesus", "Marmoset")){
	features <- ValidateFeatures(features = rownames(exl_up[[sp]]), all.genes = rownames(avgs))

	file_name <- "ExN_top_dem"
	CirclePlot.horizontal(avg = avgs, ratio = ratios, features = features, dot.min = 0.05, file_name = paste0(file_name, "_byradius_", sp), dot.scale = 6, scale.min = 5, scale.max = 100, scale.by = "radius", cluster.order = cls_order, stroke.size = 0, stroke.matrix = NULL, width.scale = 1 * 0.4, height.base = 2, font.scale = 1 * 1, height.unit = switch(sp, Human = 0.2, Chimpanzee = 0.28, Rhesus = 0.15, Marmoset = 0.1)) 
#}



file_name <- "ExN_calcium_ion_ECM"
gene_list <- list(Human = c("PKD2L1", "ANXA1", "PLS1", "SMOC1", "VCAN", "CALB1", "S100A10", "EFEMP1", "ANXA3", "SMOC2", "HMCN1", "CASQ2", "COL16A1", "ADTRP", "ADAM12", "COL9A1", "SERPINF1", "SPP1", "ITGA6", "LAMA3", "COL13A1", "TNC", "LAMC2", "COL5A1", "ADAMTS20"))
sp <- "Human"
CirclePlot.horizontal(avg = avgs, ratio = ratios, features = gene_list[[sp]], dot.min = 0.05, file_name = paste0(file_name, "_", sp, "_byradius"), dot.scale = 6, scale.min = 5, scale.max = 100, scale.by = "radius", cluster.order = cls_order, stroke.size = 0, stroke.matrix = NULL, width.scale = 1 * 0.4, height.base = 1.5, font.scale = 1 * 1, height.unit = 0.15)


















