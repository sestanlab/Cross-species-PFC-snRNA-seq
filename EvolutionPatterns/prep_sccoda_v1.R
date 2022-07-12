library(dplyr)
library(ggplot2)
library(limma)


meta <- readRDS(file = paste0("../data/", "PFC.filtered.meta.05082022.rds")) %>%
			mutate(fig1cluster = as.character(fig1cluster)) %>%
			mutate(fig1cluster = gsub("OPC PDGFRA TOP2A", "OPC PDGFRA PCDH15", fig1cluster))
meta$group2 <- meta$group
meta$group2[meta$mres %in% "Astro"] <- "Astro"
meta$group2[meta$mres %in% c("Immune", "Micro")] <- "Immune"
meta$group2[meta$mres %in% c("OPC", "Oligo")] <- "Oligo"
meta$group2[meta$mres %in% c("RB", "Endo", "SMC", "PC", "VLMC")] <- "Vas"



## Write prop for each group
all_gps <- c("ExN", "InN", "Astro", "Oligo", "Immune", "Vas")
for (gp in all_gps){
	submeta <- meta[, c("samplename", "fig1cluster", "group2", "species")] %>%
					filter(group2 == gp)
	size_df <- table(submeta$samplename, as.character(submeta$fig1cluster)) %>%
				as.data.frame() %>%
				reshape2::dcast(., Var1 ~ Var2) %>%
				rename(samplename = Var1)
	write.table(size_df, file = paste0("./load_files/", "HRES.allsps.", gp, ".meta.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}






















