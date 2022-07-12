## Plot the transcriptomic divergences change along the OPC-Oligo lineage.
source("../scripts/pfc.fun.R") 


div_mat <- readRDS(file = paste0("../MF2_globalevo/load_files/", "DIV_1minusCor_raw.rds"))

cls_use <- c("OPC PDGFRA PCDH15", "COP GPR17 SOX4", "eOligo MOG FRY", "mOligo MOG OPALIN", "lOligo MOG GSN")

plot_data <- div_mat[cls_use, ] %>%
				as.matrix() %>%
				reshape2::melt() %>%
				setNames(., c("cluster", "pair", "div")) %>%
				mutate(cluster = as.character(cluster), pair = as.character(pair)) %>%
				mutate(pair = factor(pair, levels = c("Human|Chimpanzee", "Human|Rhesus", "Chimpanzee|Rhesus", "Rhesus|Marmoset", "Human|Marmoset", "Chimpanzee|Marmoset"))) %>%
				mutate(cluster = factor(cluster, levels = cls_use))


p <- ggplot(plot_data, aes(x = cluster, y = div, color = pair, group = pair)) + 
		geom_line(size = 0.5) + 
		geom_point(size = 2, shape = 16) + 
		theme_bw() +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.text.x = element_blank(), axis.title = element_blank()) 

pdf(paste0(outputdir, "OPC-Oligo_div_lineage.pdf"), width = 6, height = 2.5, useDingbats = FALSE)
print(p)
dev.off()



