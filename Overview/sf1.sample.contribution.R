## Plot SF1 quality overview
source("../scripts/pfc.fun.R")
library(lemon)




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




meta_use <- readRDS(file = paste0(dataDir, "PFC.filtered.meta.05082022.rds")) %>%
                mutate(cluster = as.character(fig1cluster))


group.by <- "cluster"
split.by <- "species"
rep.by <- "repname"
ind.by <- "samplename"
sp_colors <- c("#FF420E", "#4CB5F5", "#89DA59", "#FFBB00") %>% setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset"))
dend <- readRDS(file = paste0("../MF1_overview/load_files/MF1_plot_dend_v2.rds"))
cls_order <- rev(labels(dend))
cluster_colors <- rep(c("#C51B7D", "#4D9221", "#984ea3", "#999999"), c(40, 46, 14, 15)) %>% setNames(., cls_order)


###----------------------------------------------------------------------------------------
## Plot sample contribution 
plot_data <- meta_use %>%
				rownames_to_column("cell") %>%
				group_by(!!sym(ind.by), !!sym(group.by)) %>%
				summarize(size = n(), !!split.by := unique(!!sym(split.by))) %>%
				ungroup() %>% 
				group_by(!!sym(group.by)) %>%
				mutate(size = size * 100/sum(size)) %>%
				ungroup() %>%
				mutate(!!group.by := factor(!!sym(group.by), levels = cls_order)) %>%
				mutate(!!split.by := factor(!!sym(split.by), levels = rev(names(sp_colors))))


mm <- ggplot(plot_data, aes_string(x = group.by, y = "size")) +
		geom_bar(aes_string(x = group.by, y = "size", fill = split.by), color = "black", position = position_stack(reverse = FALSE), stat = "identity", lwd = 0.1) +
		coord_capped_cart(bottom='both', left='both') +
		scale_fill_manual(values = sp_colors) +
		theme_classic() + 
		labs(y = "Cluster", x = "Sample Ratio") +
		theme(axis.line=element_line(size = 0.25), axis.ticks=element_line(size = 0.25), axis.text.x = element_text(size = 7, color = cluster_colors, angle = 90, hjust = 1))


pdf(paste0(outputdir, "Sample_contribution_v4.pdf"), width = 14, height = 5)
print(mm)
dev.off()



## Show the subclass information for the subtypes
allmres <- c("L2-3 IT", "L6 IT-2", "L6 IT-1", "L3-5 IT-2", "L3-5 IT-1", "L3-5 IT-3", "L5 ET", "L5-6 NP", "L6 CT", "L6B", 
	"LAMP5 LHX6", "LAMP5 RELN", "VIP", "ADARB2 KCNG1", "SST NPY", "SST", "SST HGF", "PVALB", "PVALB ChC", 
	"Astro", "OPC", "Oligo", "Micro", "Immune", "Endo", "RB", "PC", "SMC", "VLMC")
mres_data <- meta_use %>%
            group_by(fig1cluster) %>%
            summarize(mres = unique(mres)) %>%
            mutate(fig1cluster = factor(as.character(fig1cluster), levels = cls_order))
mres_sum <- mres_data[match(cls_order, mres_data$fig1cluster), ] %>%
            group_by(mres) %>%
            summarize(ncls = n()) %>%
            ungroup()
mres_sum <- mres_sum[match(allmres, mres_sum$mres), ]
mres_cusum <- cumsum(mres_sum$ncls)

pdata2 <- data.frame(cluster = allmres,
					color = rep(c("#C51B7D", "#4D9221", "#984ea3", "#999999"), c(10, 9, 4, 6)),
					xleft = c(0, mres_cusum[1:(length(allmres) - 1)]) + 0.1,
					xright = as.numeric(mres_cusum) - 0.1,
					ybottom = 0,
					ytop = 1,
					stringsAsFactors = FALSE)
p2 <- ggplot(pdata2) +
			geom_rect(aes(xmin=xleft, xmax=xright, ymin=ybottom, ymax=ytop, fill = color), color = NA)+
			scale_fill_identity()+
			theme_classic()
pdf(paste0(outputdir, "Sample_contribution_v4_mres.pdf"), width = 14, height = 5)
print(p2) 
dev.off()			












