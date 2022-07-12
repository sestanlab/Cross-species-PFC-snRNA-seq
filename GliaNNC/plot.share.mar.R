## This script will Plot the species-conserved markers across clusters.
source("../scripts/pfc.fun.R")
source("../scripts/exp.fun.R")



update_name <- function(x){
    x <- gsub("VLMC COL1A2 ABCA8", "VLMC1 COL1A2 ABCA8", x) %>%
            gsub("VLMC COL1A2 SLC13A3", "VLMC2 COL1A2 SLC13A3", x) %>%
            extract_field(., 1, " ") %>%
            setNames(., NULL) %>%
            gsub("ABC", "VLMC3", .) 
    return(xx)
}



object <- readRDS(file = paste0("../data/", "DotHeatPlot.pre.rds"))


cls_ord <- readRDS(file = paste0(dataDir, "PFC.cluster.order.rds"))[c(87:115)]
cls_ord <- cls_ord[c(1:3, 6:7,9:12,15:17,19:26, 28, 29, 27)]

genes <- c("GFAP", "AQP4", "PDGFRA", "GPR17", "OPALIN", "MOG", "PTPRC", "P2RY12", "F13A1", "LSP1", "SKAP1", "CLDN5", "DKK2", "SLC7A5", "IL1R1", "HBA1", "GRM8", "ACTA2", "CNN1", "COL1A2")

CirclePlot.horizontal(avg = object$hres$avg, ratio = object$hres$ratio, features = genes, dot.min = 0.05, file_name = paste0("NNC_share_mars", "_radius"), dot.scale = 3.5, scale.by = "radius", cluster.order = cls_ord, stroke.size = 0, stroke.matrix = NULL, width.scale = 0.5, height.base = 2.5, height.unit = 0.075)



















