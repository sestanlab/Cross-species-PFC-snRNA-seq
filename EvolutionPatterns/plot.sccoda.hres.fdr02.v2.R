## Plot sccoda results
library(dplyr)
library(ggplot2)




pair_ord <- c("Human_Chimpanzee", "Human_Rhesus", "Human_Marmoset",
			"Chimpanzee_Rhesus", "Chimpanzee_Marmoset", 
			"Rhesus_Marmoset")
all_ctps <- c("ExN", "InN", "Astro", "Oligo", "Immune", "Vas")

flist <- rep(pair_ord, each = length(all_ctps)) %>%
		paste0("HRES.new.", ., "_", rep(all_ctps, length(pair_ord)), ".meta.sccoda.effects.tsv")


pdata <- lapply(flist, function(fname) {
	pair <- strsplit(fname, ".", fixed = TRUE)[[1]][3] %>%
				strsplit(., "_", fixed = TRUE) %>%
				.[[1]] %>% .[1:2] %>% paste(., collapse = "_")
	sp1 <- strsplit(pair, "_", fixed = TRUE)[[1]][1]
	sp2 <- strsplit(pair, "_", fixed = TRUE)[[1]][2]


	tb <- read.table(paste0("./load_files/SCCODA_v2_defref/", fname), sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
			mutate(treatment = gsub("species\\[T\\.", "", Covariate)) %>%
			mutate(treatment = gsub("\\]", "", treatment)) %>%
			select(Cell.Type, Final.Parameter, log2.fold.change, treatment)
	tb$treatment <- gsub("HSB", "Human", tb$treatment) %>%
					gsub("PTB", "Chimpanzee", .) %>%
					gsub("RMB", "Rhesus", .) %>%
					gsub("CJB", "Marmoset", .)
	tb$reference <- setdiff(c(sp1, sp2), unique(tb$treatment))
	tb$pair <- paste0(tb$treatment, "_", tb$reference)
	if (unique(tb$reference) == sp1){
		tb2 <- tb
		tb2$Final.Parameter <- 0 - tb2$Final.Parameter
		tb2$log2.fold.change <- 0 - tb2$log2.fold.change
		tb2$reference <- sp2
		tb2$pair <- paste0(sp1, "_", sp2)
	} else {
		tb2 <- tb
		tb2$Final.Parameter <- 0 - tb2$Final.Parameter
		tb2$log2.fold.change <- 0 - tb2$log2.fold.change
		tb2$reference <- sp2
		tb2$pair <- paste0(sp2, "_", sp1)
	}
	tb_final <- rbind(tb, tb2)
	return(tb_final)
	}) %>%
		do.call(rbind, .)
pdata$Final.Parameter[pdata$Final.Parameter > 1] <- 1
pdata$Final.Parameter[pdata$Final.Parameter < -1] <- -1

pdata$log2.fold.change[pdata$Final.Parameter == 0] <- 0
pdata$log2.fold.change[pdata$log2.fold.change > 3] <- 3
pdata$log2.fold.change[pdata$log2.fold.change < -3] <- -3


all_pairs <- c("Human_Chimpanzee", "Human_Rhesus", "Human_Marmoset",
			"Chimpanzee_Human", "Chimpanzee_Rhesus", "Chimpanzee_Marmoset", 
			"Rhesus_Human", "Rhesus_Chimpanzee", "Rhesus_Marmoset",
			"Marmoset_Human", "Marmoset_Chimpanzee", "Marmoset_Rhesus")
## Rename cluster names
updata_name <- function(x) {
    xx <- gsub("L5 BCL11B SYT2 POU3F1", "L5 FEZF2 BCL11B POU3F1", x) %>%
            gsub("InN SST TH GABRQ", "InN SST HGF GABRQ", .) %>%
            gsub("InN LHX6 TH STON2", "InN LHX6 HGF STON2", .) %>%
            gsub("iAstro", "Astro", .) %>%
            gsub("fAstro", "Astro", .) %>%
            gsub("pAstro", "Astro", .) %>%
            gsub("rAstro", "Astro", .) %>%
            gsub("cOPC", "OPC", .) %>%
            gsub("nOligo", "Oligo", .) %>%
            gsub("eOligo", "Oligo", .) %>%
            gsub("mOligo", "Oligo", .) %>%
            gsub("lOligo", "Oligo", .) %>%
            gsub("aEndo", "Endo", .) %>%
            gsub("vEndo", "Endo", .) %>%
            gsub("cEndo", "Endo", .) %>%
            gsub("aaSMC", "SMC", .) %>%
            gsub("aSMC", "SMC", .) %>%
            gsub("vSMC", "SMC", .) %>%
            gsub("ABC", "VLMC",.) %>%
            gsub("L6b", "L6B",.)
    return(xx)
}
dend <- readRDS(file = paste0("../MF1_overview/load_files/MF1_plot_dend_v3.rds"))
cls_ord <- rev(labels(dend))


p <- ggplot(pdata, aes_string(x = "pair", y = "Cell.Type", fill = "log2.fold.change")) +
                geom_tile(width = 1, height = 1, size = 0.1, color = "black") +
                scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
                theme_classic() +
                scale_x_discrete(limits = all_pairs) +
                scale_y_discrete(limits = rev(cls_ord)) +
                coord_fixed() +
                theme(legend.position = "right", axis.text.x = element_text(size = rel(0.8), angle = 45, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = rel(0.8)))
pdf(paste0("./report/", "Prop.difs.hres.custom_ref.FDR-0.2.pdf"), width = 5, height = 10)
print(p)
dev.off()

saveRDS(pdata, file = paste0(inputdir, "Prop.difs.hres.custom_ref.FDR-02.rds"))




