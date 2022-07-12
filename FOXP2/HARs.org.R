source("./de.peak.fun.R")
library(dplyr)
library(ggplot2)


foxp2_coords_x1.5 <- 'chr7-113937174-114840607'
foxp2_chr <- 'chr7'
foxp2_start <- 113937174
foxp2_end <- 114840607



## Organize HARs peaks
pklist <- list(HARs = c("hg38_prabhakar", "hg38_bush", "ANC_Bird2007.hg38", "hg38_pollard", "2xHARs_Lindblad2013.hg18.bed.hg38", "2xPARs_Lindblad2013.hg18.bed.hg38"), 
				SelSites = c("hg38_yu_selected", "hg38_atkinson_snps", "hg38_maricic_snp"),
				HgainCREs = c("hg38_kanton_hum_vs_chip_organoids_differential", "hg38_vermunt_vs_chimp_prom_pfc", "hg38_vermunt_vs_macaque_prom_pfc", "haDHS.gittelman2015.hg38", 
					"hg38_vermunt_vs_chimp_enh_pfc", "hg38_vermunt_vs_macaque_enh_pfc", 
					"vermunt_human_gain_h3k27ac_cortex_hg38", "hg38_kozlenkov_glut_hum_enhancer"))


rawpks <- lapply(unlist(pklist), function(x) {
		print(x)
		tb <- read.table(paste0("./HARs/", x, ".bed"), sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>%
			.[, 1:3] %>%
			setNames(., c("chr", "start", "end")) %>%
			mutate(peaks = paste0(chr, "-", start, "-", end))
		fpks <- overlap_foxp2_peaks(unique(tb$peaks))
		return(fpks)
		}) %>%
		setNames(., unlist(pklist))


## Plot the rectangle plots 
pdata <- lapply(names(pklist), function(type) {
	base_df <- data.frame(xmin = foxp2_start, xmax = foxp2_end,
					ymin = switch(type, HgainCREs = 0, SelSites = 1, HARs = 2),
					ymax = switch(type, HgainCREs = 1, SelSites = 2, HARs = 3),
					color = "lightgrey",
					stringsAsFactors = FALSE)

	cur_pks <- rawpks[pklist[[type]]] %>% unlist() %>% unique()

	pk_df <- data.frame(peaks = cur_pks, stringsAsFactors = FALSE) %>%
					mutate(chr = extract_field(peaks, 1, "-"), xmin = extract_field(peaks, 2, "-"), 
						xmax = extract_field(peaks, 3, "-")) %>%
					select(xmin, xmax) %>%
					mutate(ymin = switch(type, HgainCREs = 0, SelSites = 1, HARs = 2)) %>%
					mutate(ymax = switch(type, HgainCREs = 1, SelSites = 2, HARs = 3)) %>%
					mutate(color = "black") %>%
					.[, colnames(base_df)] %>%
					rbind(base_df, .) %>%
					mutate(type = type)
	pk_df
	}) %>%
	do.call(rbind, .) %>%
	mutate(xmin = as.numeric(xmin), xmax = as.numeric(xmax))


p <- ggplot(pdata, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))+
		geom_rect(aes(fill = color), color = NA, size = 1) +
		scale_fill_identity() + 
		theme_classic() +
		scale_x_continuous(breaks = seq(113900000, 114700000, 200000), limits = c(foxp2_start, foxp2_end)) +
		theme(axis.line = element_line(size = 0.2), 
			axis.ticks = element_line(size = 0.2))
pdf("./report/DE_peaks_intersect_HARs.pdf", 6, 3)
print(p)
dev.off()




