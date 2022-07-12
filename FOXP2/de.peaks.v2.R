library(Seurat)
library(Signac)
library(DT)
library(tibble)
library(ggplot2)
library(dplyr)
source("./de.peak.fun.R")



cls_transfer <- c("L2-3 IT","L3-5 IT-1","L3-5 IT-2","L3-5 IT-3","L6 IT-1","L6 IT-2","L5 ET","L5-6 NP","L6 CT","L6B", "LAMP5 LHX6", "LAMP5 RELN", "VIP", "ADARB2 KCNG1", "SST NPY", "SST", "SST HGF", "PVALB","PVALB ChC", "Astro","OPC","Oligo","Micro","Immune","Endo","PC","SMC","VLMC") %>%
            setNames(., c("L2_3_IT","L3_5_IT_1","L3_5_IT_2","L3_5_IT_3","L6_IT_1","L6_IT_2","L5_PT","L5_6_NP","L6_CT","L6b", "LAMP5_LHX6", "LAMP5_RELN", "VIP", "ADARB2", "SST_NPY", "SST", "TH", "PVALB","PVALB_CHC", "Astro","OPC","Oligo","Micro","immune","Endo","PC","SMC","VLMC"))

## Load data
object <- readRDS("/home/sm2726/myshare/AdultHumanMultiome/CombinedSeurat/Four.smps.adultPFC.addUMAP.v05092022.rds")
object <- subset(object, ATAC_QC == 1)  ## Subset to only high quality cells


object@meta.data$subclass <- factor(setNames(cls_transfer[object@meta.data$mres], NULL), levels = setNames(cls_transfer, NULL))
Idents(object) <- object@meta.data$subclass
object@assays$RNA<-NULL
object@assays$ARCRNA<-NULL
DefaultAssay(object)<-"peaks"



## Find Differentially Accessible Peaks
## Overlaps with FOXP2 coordinate
foxp2_coords <- "chr7-114087747-114690035"
foxp2_coords_x1.5 <- 'chr7-113937174-114840607'
foxp2_chr <- 'chr7'
foxp2_start <- 113937174
foxp2_end <- 114840607



#########################
## Define groups
object$peak_region_fragments <- GetAssayData(object = object, assay = "peaks", slot = "counts") %>%
				colSums()


## Find all peaks within the FOXP2 region
input_peaks <- overlap_foxp2_peaks(rownames(object$peaks))
peaks_neuron_res <- lapply(c("L3-5 IT-1", "L3-5 IT-2", "L3-5 IT-3", "L5-6 NP", "L6 CT", "L6B", "Micro"), function(cls){
	res <- FindMarkers(object = object, features = input_peaks, ident.1 = cls, ident.2 = NULL, min.pct = 0.025, test.use = 'LR', latent.vars = 'peak_region_fragments', logfc.threshold = 0.01, only.pos = TRUE) %>%
			rownames_to_column("peak") %>%
			mutate(cluster = cls)
	return(res) 
	}) %>%
	do.call(rbind, .)


peaks_neuron_fil <- peaks_neuron_res %>%
			filter(p_val_adj <= 0.01) %>%
			mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01))
peaks_neuron_fil <- peaks_neuron_fil[(peaks_neuron_fil$pct.1 >= 0.1 & peaks_neuron_fil$ratio_fc >= 1.5) |
								(peaks_neuron_fil$pct.1 < 0.1 & peaks_neuron_fil$ratio_fc >= 2), ]
peaks_neuron_fil <- peaks_neuron_fil %>%
			##filter(peak != "chr7-114081562-114088841") %>%
			mutate(start = as.numeric(extract_field(peak, 2, "-"))) %>%
			mutate(end = as.numeric(extract_field(peak, 3, "-")))


uniq_peaks <- unique(peaks_neuron_fil$peak)
uniq_peaks <- uniq_peaks[order(as.numeric(extract_field(uniq_peaks, 2, "-")))]
peak_ord <- setdiff(uniq_peaks, c("chr7-114051725-114052632", "chr7-114081562-114088841")) %>%
				c(., c("chr7-114051725-114052632", "chr7-114081562-114088841"))

peaks_neuron_fil$peakorder <- factor(peaks_neuron_fil$peak, levels = peak_ord) %>% as.numeric() %>% paste0("DAP-", .)
write.table(peaks_neuron_fil, file = "./load_files/Differentially_accessible_peaks_aroundFOXP2.xls", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



##----------------------------------------------------------------------
#coverage plot to check the peaks
cov.plt <- CoveragePlot(
  object = object,
  assay.scale = "common",
  region = "chr7-113937174-114840607",
  annotation = TRUE,
  peaks = TRUE
)

uniq_peaks <- unique(peaks_neuron_fil$peak)
for (ii in 1:nrow(peaks_neuron_fil)){
	cov.plt <- cov.plt +
			annotate("rect", xmin = as.numeric(extract_field(uniq_peaks[ii], 2, "-")), xmax = as.numeric(extract_field(uniq_peaks[ii], 3, "-")), ymin = 1, ymax = 3.6, alpha = .2, color = NA, fill = "red")
}


pdf('report/FOXP2.CoveragePlot.pdf', width =10, height = 5)
print(cov.plt)
dev.off()



##----------------------------------------------------------------------
## Highlight Differntially-accessible peaks
high_region <- sapply(uniq_peaks, function(pk) {
	start <- as.numeric(extract_field(pk, 2, "-"))
	end <- as.numeric(extract_field(pk, 3, "-"))

	nstart <- start - 3500 * 5
	nend <- end + 3500 * 5

	paste0("chr7-", nstart, "-", nend)
	}) %>%
	setNames(., NULL)
high_peaks <- uniq_peaks[order(as.numeric(extract_field(uniq_peaks, 2, "-")))]


high_reg1 <- paste0("chr7-", 
			as.numeric(extract_field(high_peaks[1], 2, "-")) - 6000, 
			"-", 
			as.numeric(extract_field(high_peaks[1], 3, "-")) + 6000)
len <- as.numeric(extract_field(high_reg1, 3, "-")) - as.numeric(extract_field(high_reg1, 2, "-"))
p1 <- CoveragePlot(object = object, assay.scale = "common", region = high_reg1, annotation = TRUE, peaks = TRUE, ncol = 3)
for (ii in 1:1){
	p1 <- p1 +
			annotate("rect", xmin = as.numeric(extract_field(high_peaks[ii], 2, "-")), xmax = as.numeric(extract_field(high_peaks[ii], 3, "-")), ymin = 1, ymax = 3.6, alpha = .2, color = NA, fill = "red")
}
pdf('report/FOXP2.Fil.CoveragePlot.zoom-in-1.pdf', width = 3, height = 5)
print(p1)
dev.off()

