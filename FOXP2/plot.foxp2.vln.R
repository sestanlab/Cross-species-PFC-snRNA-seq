source("../scripts/pfc.fun.R")
source("./foxp2.fun.R") 



plot.pre <- paste0(dataDir, "plot.HPRC.Mouse.ExonIntron.plot.Rdata")
if (!file.exists(plot.pre)){
	hprc <- readRDS(file = paste0(dataDir, "PFC_filtered_seu.11052020.rds"))
	mctx <- readRDS(file = paste0(dataDir, "public.MouseVISp-ALM.org.rds"))


	hprc_data <- hprc$RNA@data
	## Use the exon + intron data for visualiz expression
	mvisp.eidata <- readRDS("/home/sm2726/project/public_data/Allen_mouse_datasets_2019_09_26/load_files/1a_mouse_VISp_both_both_seu.rds")$RNA@data
	malm.eidata <- readRDS("/home/sm2726/project/public_data/Allen_mouse_datasets_2019_09_26/load_files/1a_mouse_ALM_both_both_seu.rds")$RNA@data
	colnames(mvisp.eidata) <- paste0("VISP_", colnames(mvisp.eidata))
	colnames(malm.eidata) <- paste0("ALM_", colnames(malm.eidata))
	cbn.data <- cbind(mvisp.eidata, malm.eidata)
	rownames(cbn.data) <- toupper(rownames(cbn.data))
	mctx_data <- cbn.data[, colnames(mctx)]

	cls_order <- c("L2-3 IT","L3-5 IT-1","L3-5 IT-2","L3-5 IT-3","L5 PT","L5-6 NP","L6 CT","L6 IT-1","L6 IT-2","L6b","LAMP5 LHX6","LAMP5 RELN","ADARB2","VIP","SST NPY","SST","TH","PVALB","PVALB CHC","Astro","OPC","Oligo","Micro","immune","blood","Endo","PC","SMC","VLMC")
	mrestransfer <- unique(hprc@meta.data$mres) %>%
						gsub("^(L[0-9])_([0-9].*)$", "\\1-\\2", .) %>%
						gsub("^(.*IT)_([0-9])$", "\\1-\\2", .) %>%
						gsub("_", " ", .) %>%
						setNames(., unique(hprc@meta.data$mres))

	meta <- data.frame(cluster = c(mrestransfer[hprc@meta.data$mres], mctx@meta.data$mres), 
		species = c(hprc@meta.data$species, rep("Mouse", ncol(mctx))), stringsAsFactors = FALSE) %>%
					mutate(cluster = factor(cluster, levels = cls_order)) %>%
					mutate(species = factor(species, levels = c("Human", "Chimpanzee", "Rhesus", "Marmoset", "Mouse")))
	save(hprc_data, mctx_data, meta, file = plot.pre)
} else {
	load(plot.pre)
}




plotHMvln <- function(genes, output_dir, ptype = c("all", "ExN", "InN", "NNC")[2], file_name = "HPRC-Mouse.mres.exp") {
	for (gene in genes){
		message(paste0("Plotting gene: ", gene))
		if ((!gene %in% rownames(hprc_data)) && (!gene %in% rownames(mctx_data))){
			stop(paste0("gene ", gene, " are not present in hprc"))
		} else {
			## Whether the gene is present in mouse
			if (!gene %in% rownames(mctx_data)){
				sp_order <- c("Human", "Chimpanzee", "Rhesus", "Marmoset")
				plot.data <- meta %>%
							subset(species %in% sp_order) %>%
							mutate(exp = hprc_data[gene, ]) %>% 
							mutate(species = factor(as.character(species), levels = c("Human", "Chimpanzee", "Rhesus", "Marmoset")))
				pdf_height <- 4
			} else {
				plot.data <- meta %>%
							mutate(exp = c(hprc_data[gene, ], mctx_data[gene, ]))
				pdf_height <- 4
			}


			## Plot Type
			cls_use <- switch(ptype, all = cls_order, ExN = cls_order[1:7], InN = cls_order[8:12], NNC = cls_order[13:18])
			plot.data <- plot.data %>%
								subset(cluster %in% cls_use) %>%
								mutate(cluster = factor(as.character(cluster), levels = cls_use))



			## Plot Type
			me <- max(plot.data$exp)
			p <- ggplot() + 
					geom_violin(data = plot.data, aes_string(x = "cluster", y = "exp", fill = "species"), scale = "width", size = 0.1, adjust = 2,trim =TRUE) + 
					theme_cowplot() + 
					RotatedAxis() + 
					scale_fill_manual(values = c("#FF420E", "#4CB5F5", "#89DA59", "#FFBB00", "#AE017E") %>% setNames(., c("Human", "Chimpanzee", "Rhesus", "Marmoset", "Mouse"))) + 
					facet_wrap(vars(species), nrow = 5, ncol = 1, strip.position = 'left') +
					scale_y_continuous(breaks = seq(0, me, by = round(me/2, digits = 1))) + 
					theme(axis.text.x=element_text(size = 8), strip.background = element_blank(), strip.text = element_text(size = 10, angle = 90), panel.spacing = unit(0.05, "in"), strip.placement = 'outside', axis.line = element_line(size = 0.25), axis.ticks = element_line(size = 0.25), legend.position = "bottom", legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size = rel(0.5)), legend.title = element_blank())
			pdf(paste0(output_dir, file_name, ".", ptype, ".", gene, ".pdf"), width =ceiling(length(levels(plot.data$cluster)) * 0.1) + 1.5, height = pdf_height)
			print(p)
			dev.off() 
		}
	}
}




meta_bc <- meta
meta <- meta %>%
		mutate(cluster = gsub("^(L3-5 IT)-[0-9]", "\\1", cluster)) %>%
		mutate(cluster = gsub("^(L6 IT)-[0-9]", "\\1", cluster)) %>%
		mutate(cluster = gsub("^(LAMP5).*", "\\1", cluster)) %>%
		mutate(cluster = gsub("^(SST).*", "\\1", cluster)) %>%
		mutate(cluster = gsub("^(PVALB).*", "\\1", cluster)) %>%
		mutate(cluster = gsub("^TH$", "SST", cluster)) %>%
		mutate(cluster = gsub("blood|Endo|VLMC", "Vas", cluster)) %>%
		mutate(cluster = gsub("^PC$", "Vas", cluster)) %>%
		mutate(cluster = gsub("^SMC$", "Vas", cluster))


source("./foxp2.fun.R") 

cls_order <- c("L2-3 IT","L3-5 IT","L5 PT","L5-6 NP","L6 CT","L6 IT","L6b","LAMP5","ADARB2","VIP","SST","PVALB","Astro","OPC","Oligo","Micro","immune","Vas")
plotHMvln(genes = "FOXP2", output_dir = paste0(outputdir), ptype = c("all", "ExN", "InN", "NNC")[1], cls_order = cls_order,file_name = "HPRC-Mouse.lres.exp.ExIntron")







